"""
Unified Boltz-2 launcher for Vast.ai -- zero-setup predictions.

Uses custom Docker image with pre-baked Boltz-2, CCD, and model weights.
Supports persistent instances (stop/start) to avoid repeated setup.

Usage:
    python scripts/boltz_runner.py search                        # Show available GPUs
    python scripts/boltz_runner.py launch [--tier small|large]   # Launch instance
    python scripts/boltz_runner.py warm                          # Reuse/restart existing
    python scripts/boltz_runner.py run --yamls <dir> [--samples 3] [--seed 42]
    python scripts/boltz_runner.py download [--dest <dir>]       # Download CIF + JSON
    python scripts/boltz_runner.py stop                          # Pause billing (keep disk)
    python scripts/boltz_runner.py start                         # Resume stopped instance
    python scripts/boltz_runner.py status                        # Show instance info + cost
    python scripts/boltz_runner.py ssh                           # Print SSH command
    python scripts/boltz_runner.py destroy                       # Full teardown
    python scripts/boltz_runner.py full --yamls <dir> [--keep]   # End-to-end

Requires: pip install vastai
API key: set via `vastai set api-key <KEY>` or VAST_API_KEY env var
"""
import subprocess
import sys
import json
import os
import time
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

# IMPORTANT: Update this to your Docker Hub image after building
IMAGE = "jerempire/boltz2-ancient:latest"

GPU_TIERS = {
    "small": {
        "min_vram_gb": 11,      # RTX 2080 Ti minimum
        "max_price_per_hr": 0.80,
        "disk_gb": 30,
        "order": "dph",        # cheapest first
    },
    "large": {
        "min_vram_gb": 40,      # A100 40GB minimum
        "max_price_per_hr": 3.00,
        "disk_gb": 50,
        "order": "gpu_ram-",   # most VRAM first
    },
}

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCE_FILE = os.path.join(PROJECT_DIR, ".vast_instance_boltz.json")


# ── Helpers ──────────────────────────────────────────────────────────────────

def run_vastai(*args, parse_json=True):
    """Run a vastai CLI command and return output."""
    cmd = ["vastai"] + list(args)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"vastai error: {result.stderr.strip()}")
        sys.exit(1)
    if parse_json:
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError:
            return result.stdout.strip()
    return result.stdout.strip()


def save_instance(info):
    with open(INSTANCE_FILE, "w") as f:
        json.dump(info, f, indent=2)


def load_instance():
    if not os.path.exists(INSTANCE_FILE):
        print(f"No saved instance ({INSTANCE_FILE}). Run 'launch' or 'warm' first.")
        sys.exit(1)
    with open(INSTANCE_FILE) as f:
        return json.load(f)


def ssh_base(info):
    return f"ssh -o StrictHostKeyChecking=no -p {info['ssh_port']} root@{info['ssh_host']}"


def scp_to(info, local_src, remote_dst):
    return (f"scp -o StrictHostKeyChecking=no -r "
            f"-P {info['ssh_port']} {local_src} root@{info['ssh_host']}:{remote_dst}")


def scp_from(info, remote_src, local_dst):
    return (f"scp -o StrictHostKeyChecking=no -r "
            f"-P {info['ssh_port']} root@{info['ssh_host']}:{remote_src} {local_dst}")


def elapsed_hours(info):
    try:
        launched = time.mktime(time.strptime(info["launched_at"], "%Y-%m-%d %H:%M:%S"))
        return (time.time() - launched) / 3600
    except (KeyError, ValueError):
        return 0


def choose_tier(yaml_dir, diffusion_samples):
    n = len(list(Path(yaml_dir).glob("*.yaml")))
    if n > 30 or diffusion_samples >= 5:
        return "large"
    return "small"


def wait_for_running(instance_id, timeout_sec=300):
    """Poll until instance is running. Returns updated info dict."""
    for i in range(timeout_sec // 5):
        time.sleep(5)
        instances = run_vastai("show", "instances", "--raw")
        if isinstance(instances, list):
            for inst in instances:
                if str(inst.get("id")) == str(instance_id):
                    status = inst.get("actual_status", "unknown")
                    if status == "running":
                        return inst
                    print(f"  [{(i+1)*5}s] {status}...")
                    break
    print(f"Timeout after {timeout_sec}s waiting for instance {instance_id}")
    sys.exit(1)


# ── Commands ─────────────────────────────────────────────────────────────────

def cmd_search(tier_name="small"):
    tier = GPU_TIERS[tier_name]
    print(f"Searching: >= {tier['min_vram_gb']}GB VRAM, <= ${tier['max_price_per_hr']}/hr...")
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {tier['min_vram_gb']}",
        f"dph <= {tier['max_price_per_hr']}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", tier["order"],
        "--type", "on-demand",
        "--limit", "15",
        "--raw",
    )
    if not results:
        print("No instances found. Try --tier large or wait.")
        return None

    print(f"\n{'ID':>8}  {'GPU':>20}  {'VRAM':>6}  {'$/hr':>6}  {'CUDA':>5}  {'DL Mbps':>8}")
    print("-" * 70)
    for r in results:
        vram_gb = r.get('gpu_ram', 0) / 1024  # API returns MB
        print(f"{r.get('id','?'):>8}  {r.get('gpu_name','?'):>20}  "
              f"{vram_gb:>5.0f}G  ${r.get('dph_total',0):>5.2f}  "
              f"{r.get('cuda_max_good','?'):>5}  {r.get('inet_down',0):>7.0f}")
    return results


def cmd_launch(tier_name="small"):
    tier = GPU_TIERS[tier_name]
    print(f"Finding best GPU (tier={tier_name})...")
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {tier['min_vram_gb']}",
        f"dph <= {tier['max_price_per_hr']}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", tier["order"],
        "--type", "on-demand",
        "--limit", "1",
        "--raw",
    )
    if not results:
        print("No suitable instances. Try --tier large or wait.")
        sys.exit(1)

    offer = results[0]
    gpu = offer.get("gpu_name", "?")
    vram = offer.get("gpu_ram", 0) / 1024  # API returns MB
    price = offer.get("dph_total", 0)
    print(f"Launching: {gpu} ({vram:.0f}GB) at ${price:.2f}/hr")
    print(f"Image: {IMAGE}")

    result = run_vastai(
        "create", "instance", str(offer["id"]),
        "--image", IMAGE,
        "--disk", str(tier["disk_gb"]),
        "--raw",
    )

    instance_id = result.get("new_contract") if isinstance(result, dict) else str(result)
    print(f"Instance {instance_id} created. Waiting for startup...")

    inst = wait_for_running(instance_id)
    info = {
        "id": instance_id,
        "gpu": gpu,
        "vram_gb": vram,
        "price_per_hr": price,
        "ssh_host": inst.get("ssh_host", ""),
        "ssh_port": inst.get("ssh_port", ""),
        "image": IMAGE,
        "disk_gb": tier["disk_gb"],
        "tier": tier_name,
        "launched_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "weights_ready": True,  # Custom image has weights pre-baked
    }
    save_instance(info)
    print(f"\nRUNNING: {gpu} ({vram:.0f}GB) at ${price:.2f}/hr")
    print(f"SSH: ssh -p {info['ssh_port']} root@{info['ssh_host']}")
    return info


def cmd_warm(tier_name="small"):
    """Reuse existing instance if possible, otherwise launch new."""
    if not os.path.exists(INSTANCE_FILE):
        print("No saved instance. Launching fresh...")
        return cmd_launch(tier_name)

    info = load_instance()
    instances = run_vastai("show", "instances", "--raw")
    if isinstance(instances, list):
        for inst in instances:
            if str(inst.get("id")) == str(info["id"]):
                status = inst.get("actual_status", "unknown")
                if status == "running":
                    # Update SSH info in case it changed
                    info["ssh_host"] = inst.get("ssh_host", info["ssh_host"])
                    info["ssh_port"] = inst.get("ssh_port", info["ssh_port"])
                    save_instance(info)
                    print(f"Reusing instance {info['id']} ({info['gpu']}, "
                          f"running {elapsed_hours(info):.1f}h)")
                    return info
                elif status == "stopped":
                    print(f"Instance {info['id']} stopped. Restarting...")
                    run_vastai("start", "instance", str(info["id"]), parse_json=False)
                    inst_data = wait_for_running(info["id"])
                    info["ssh_host"] = inst_data.get("ssh_host", info["ssh_host"])
                    info["ssh_port"] = inst_data.get("ssh_port", info["ssh_port"])
                    save_instance(info)
                    print(f"Resumed: {info['gpu']} at ${info['price_per_hr']:.2f}/hr")
                    return info
                elif status in ("loading", "provisioning"):
                    print(f"Instance {info['id']} is {status}, waiting...")
                    inst_data = wait_for_running(info["id"])
                    info["ssh_host"] = inst_data.get("ssh_host", info["ssh_host"])
                    info["ssh_port"] = inst_data.get("ssh_port", info["ssh_port"])
                    save_instance(info)
                    return info
                else:
                    print(f"Instance {info['id']} is {status}. Launching new...")
                    break

    return cmd_launch(tier_name)


def cmd_stop():
    info = load_instance()
    hrs = elapsed_hours(info)
    est_cost = hrs * info.get("price_per_hr", 0)
    print(f"Stopping {info['gpu']} (ran {hrs:.1f}h, ~${est_cost:.2f})...")
    print("Disk preserved. Resume with: python scripts/boltz_runner.py start")
    run_vastai("stop", "instance", str(info["id"]), parse_json=False)
    print("Instance stopped. Billing paused (disk-only charges apply).")


def cmd_start():
    info = load_instance()
    print(f"Starting instance {info['id']} ({info['gpu']})...")
    run_vastai("start", "instance", str(info["id"]), parse_json=False)
    inst = wait_for_running(info["id"])
    info["ssh_host"] = inst.get("ssh_host", info["ssh_host"])
    info["ssh_port"] = inst.get("ssh_port", info["ssh_port"])
    save_instance(info)
    print(f"Running: ssh -p {info['ssh_port']} root@{info['ssh_host']}")


def cmd_status():
    instances = run_vastai("show", "instances", "--raw")
    if not instances:
        print("No running instances.")
        return

    # Check if we have a saved instance
    saved_id = None
    if os.path.exists(INSTANCE_FILE):
        saved_id = str(load_instance().get("id"))

    for inst in instances:
        iid = str(inst.get("id", "?"))
        gpu = inst.get("gpu_name", "?")
        status = inst.get("actual_status", "?")
        dph = inst.get("dph_total", 0)
        marker = " <- saved" if iid == saved_id else ""
        print(f"  {iid}: {gpu} | {status} | ${dph:.2f}/hr{marker}")


def cmd_ssh():
    info = load_instance()
    print(f"ssh -o StrictHostKeyChecking=no -p {info['ssh_port']} root@{info['ssh_host']}")


def cmd_destroy():
    info = load_instance()
    hrs = elapsed_hours(info)
    est_cost = hrs * info.get("price_per_hr", 0)
    print(f"Destroying {info['gpu']} (ran {hrs:.1f}h, ~${est_cost:.2f})...")
    run_vastai("destroy", "instance", str(info["id"]), parse_json=False)
    os.remove(INSTANCE_FILE)
    print("Instance destroyed.")


def cmd_run(yaml_dir, out_dir="/workspace/results", diffusion_samples=3, seed=42):
    """Upload YAMLs and run Boltz-2 predictions on the remote instance."""
    info = load_instance()
    ssh = ssh_base(info)
    yaml_path = Path(yaml_dir)

    if not yaml_path.exists():
        print(f"YAML directory not found: {yaml_dir}")
        sys.exit(1)

    yamls = sorted(yaml_path.glob("*.yaml"))
    if not yamls:
        print(f"No .yaml files in {yaml_dir}")
        sys.exit(1)

    n = len(yamls)
    est_min = n * diffusion_samples * 0.5  # ~30sec per sample on RTX 4090

    print(f"Uploading {n} YAMLs to instance...")
    remote_yaml = "/workspace/boltz_inputs"
    subprocess.run(f'{ssh} "mkdir -p {remote_yaml} {out_dir}"', shell=True, check=True)

    # Upload all YAMLs
    subprocess.run(
        scp_to(info, f'"{yaml_path}/"*.yaml', f"{remote_yaml}/"),
        shell=True, check=True,
    )

    print(f"Running {n} complexes x {diffusion_samples} samples (~{est_min:.0f} min est.)")
    print(f"GPU: {info['gpu']} ({info.get('vram_gb', '?')}GB)")
    print()

    # Build the remote run script
    # Uses BOLTZ_CACHE_DIR to find pre-baked weights/CCD
    # Skips already-completed predictions for resume capability
    run_script = f'''set -e
REMOTE_YAML="{remote_yaml}"
OUT_DIR="{out_dir}"
SAMPLES={diffusion_samples}
SEED={seed}

DONE=0
FAILED=0
SKIPPED=0
TOTAL=$(ls $REMOTE_YAML/*.yaml 2>/dev/null | wc -l)

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No YAML files found in $REMOTE_YAML"
    exit 1
fi

echo "=== Boltz-2 batch run: $TOTAL complexes x $SAMPLES samples ==="
echo "Started: $(date)"
echo ""

for yaml in $REMOTE_YAML/*.yaml; do
    name=$(basename "$yaml" .yaml)
    out="$OUT_DIR/$name"
    DONE=$((DONE + 1))

    # Skip if already completed
    if [ -d "$out" ] && find "$out" -name "confidence_*.json" -print -quit 2>/dev/null | grep -q .; then
        echo "[$DONE/$TOTAL] $name - SKIP (done)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo -n "[$DONE/$TOTAL] $name ... "
    START_T=$(date +%s)

    if BOLTZ_CACHE=/opt/boltz boltz predict "$yaml" \\
        --diffusion_samples $SAMPLES \\
        --seed $SEED \\
        --out_dir "$out" \\
        --accelerator gpu \\
        2>&1 | tail -2; then
        END_T=$(date +%s)
        echo "  OK ($((END_T - START_T))s)"
    else
        echo "  FAILED"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Summary ==="
echo "Total: $TOTAL | Completed: $((DONE - SKIPPED - FAILED)) | Skipped: $SKIPPED | Failed: $FAILED"
echo "Finished: $(date)"
'''

    # Run it
    subprocess.run(f"""{ssh} 'bash -s' << 'REMOTE_SCRIPT'\n{run_script}\nREMOTE_SCRIPT""",
                   shell=True)


def cmd_download(dest=None, remote_dir="/workspace/results"):
    """Download results including CIF files (mandatory)."""
    info = load_instance()

    if dest is None:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        dest = os.path.join(PROJECT_DIR, "data", "results", "boltz2_complexes", timestamp)
    os.makedirs(dest, exist_ok=True)

    print(f"Downloading results from {info['gpu']}...")
    subprocess.run(
        scp_from(info, f"{remote_dir}/*", f'"{dest}/"'),
        shell=True,
    )

    # Verify downloads -- CIF files are mandatory per CLAUDE.md
    dest_path = Path(dest)
    cif_files = list(dest_path.rglob("*.cif"))
    json_files = list(dest_path.rglob("confidence_*.json"))
    print(f"\nDownloaded: {len(cif_files)} CIF files, {len(json_files)} confidence JSON files")

    if not cif_files:
        print("WARNING: No CIF files found! PyRosetta analysis will fail.")
        print("Check remote: ssh into instance and run: find /workspace/results -name '*.cif'")

    if json_files:
        # Quick summary of scores
        print(f"\nResults saved to: {dest}")
        _print_score_summary(dest_path)
    else:
        print("WARNING: No confidence JSON files found. Predictions may have failed.")

    return dest


def _print_score_summary(results_dir):
    """Print a quick ipTM summary from downloaded results."""
    import json as json_mod
    scores = []
    for jf in results_dir.rglob("confidence_*.json"):
        try:
            data = json_mod.loads(jf.read_text())
            if isinstance(data, list):
                data = data[0]
            if isinstance(data, dict):
                iptm = data.get("iptm") or data.get("ipTM") or data.get("interface_ptm")
                if iptm is not None:
                    # Get complex name from directory structure
                    name = jf.parent.parent.parent.name
                    scores.append((name, float(iptm)))
        except Exception:
            pass

    if scores:
        scores.sort(key=lambda x: x[0])
        print(f"\n{'Complex':>35s}  {'ipTM':>8s}")
        print("-" * 48)
        for name, iptm in scores:
            print(f"{name:>35s}  {iptm:>8.4f}")


def cmd_full(yaml_dir, diffusion_samples=3, seed=42, keep=False, tier_name=None):
    """End-to-end: warm/launch -> run -> download -> stop/destroy."""
    if tier_name is None:
        tier_name = choose_tier(yaml_dir, diffusion_samples)
    print(f"Tier: {tier_name}")

    info = cmd_warm(tier_name)
    print("\nWaiting 10s for SSH to stabilize...")
    time.sleep(10)

    print()
    cmd_run(yaml_dir, diffusion_samples=diffusion_samples, seed=seed)

    print()
    local_dir = cmd_download()

    # Cost report
    hrs = elapsed_hours(info)
    cost = hrs * info.get("price_per_hr", 0)
    print(f"\nCost: ~${cost:.2f} ({hrs:.1f}h at ${info.get('price_per_hr', 0):.2f}/hr)")

    if keep:
        print("\nStopping instance (--keep). Resume later with: boltz_runner.py start")
        cmd_stop()
    else:
        print("\nDestroying instance.")
        cmd_destroy()

    print(f"Results: {local_dir}")
    return local_dir


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Boltz-2 on Vast.ai -- zero-setup predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command")

    # search
    p = sub.add_parser("search", help="Show available GPUs")
    p.add_argument("--tier", choices=["small", "large"], default="small")

    # launch
    p = sub.add_parser("launch", help="Launch new instance")
    p.add_argument("--tier", choices=["small", "large"], default="small")

    # warm
    p = sub.add_parser("warm", help="Reuse existing or launch new")
    p.add_argument("--tier", choices=["small", "large"], default="small")

    # run
    p = sub.add_parser("run", help="Upload YAMLs and run predictions")
    p.add_argument("--yamls", required=True, help="Directory of .yaml input files")
    p.add_argument("--samples", type=int, default=3, help="Diffusion samples (default: 3)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

    # download
    p = sub.add_parser("download", help="Download results (CIF + JSON)")
    p.add_argument("--dest", help="Local destination directory")

    # stop / start / status / ssh / destroy
    sub.add_parser("stop", help="Pause billing (keep disk)")
    sub.add_parser("start", help="Resume stopped instance")
    sub.add_parser("status", help="Show instance info")
    sub.add_parser("ssh", help="Print SSH command")
    sub.add_parser("destroy", help="Destroy instance")

    # full
    p = sub.add_parser("full", help="End-to-end: launch -> run -> download -> cleanup")
    p.add_argument("--yamls", required=True, help="Directory of .yaml input files")
    p.add_argument("--samples", type=int, default=3, help="Diffusion samples (default: 3)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--keep", action="store_true", help="Stop (not destroy) after run")
    p.add_argument("--tier", choices=["small", "large"], default=None)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    dispatch = {
        "search": lambda: cmd_search(args.tier),
        "launch": lambda: cmd_launch(args.tier),
        "warm": lambda: cmd_warm(args.tier),
        "run": lambda: cmd_run(args.yamls, diffusion_samples=args.samples, seed=args.seed),
        "download": lambda: cmd_download(dest=args.dest),
        "stop": cmd_stop,
        "start": cmd_start,
        "status": cmd_status,
        "ssh": cmd_ssh,
        "destroy": cmd_destroy,
        "full": lambda: cmd_full(
            args.yamls,
            diffusion_samples=args.samples,
            seed=args.seed,
            keep=args.keep,
            tier_name=args.tier,
        ),
    }

    dispatch[args.command]()


if __name__ == "__main__":
    main()
