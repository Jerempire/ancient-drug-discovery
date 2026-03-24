"""Launch Vast.ai instance and run validation experiment end-to-end.

Targets H100 > A100 > L40S for fastest throughput.

Usage:
    python v4/validation/launch_validation.py           # Search, launch, upload, run, download, destroy
    python v4/validation/launch_validation.py search    # Just search for available GPUs
    python v4/validation/launch_validation.py upload    # Upload to existing instance
    python v4/validation/launch_validation.py run       # SSH and run experiment
    python v4/validation/launch_validation.py download  # Download results
    python v4/validation/launch_validation.py destroy   # Tear down instance
"""
import subprocess
import sys
import json
import os
import time

# --- Config: target fastest GPUs ---
MIN_VRAM_GB = 40              # H100=80, A100=40/80, L40S=48
MAX_PRICE_PER_HR = 4.00       # H100 ~$2-3/hr, A100 ~$1-2/hr on Vast.ai
DISK_GB = 50
IMAGE = "nvidia/cuda:12.4.1-devel-ubuntu22.04"  # Generic CUDA image (Boltz-2 pip-installed)
PREFERRED_ORDER = "gpu_ram-"  # Prefer most VRAM (H100 > A100 80GB > A100 40GB)

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
VALIDATION_DIR = os.path.dirname(os.path.abspath(__file__))
INSTANCE_FILE = os.path.join(VALIDATION_DIR, ".vast_instance.json")


def run_vastai(*args, parse_json=True):
    cmd = ["vastai"] + list(args)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr.strip()}")
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
        print("No instance found. Run 'launch_validation.py' without args first.")
        sys.exit(1)
    with open(INSTANCE_FILE) as f:
        return json.load(f)


def ssh_cmd(info):
    return f"ssh -o StrictHostKeyChecking=no -p {info['ssh_port']} root@{info['ssh_host']}"


def scp_cmd(info, src, dst):
    return f"scp -o StrictHostKeyChecking=no -r -P {info['ssh_port']} {src} root@{info['ssh_host']}:{dst}"


def cmd_search():
    print(f"Searching: >= {MIN_VRAM_GB}GB VRAM, <= ${MAX_PRICE_PER_HR}/hr, fastest GPUs...")
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {MIN_VRAM_GB}",
        f"dph <= {MAX_PRICE_PER_HR}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", "gpu_ram-",
        "--type", "on-demand",
        "--limit", "15",
        "--raw",
    )
    if not results:
        print("No instances found. Try increasing MAX_PRICE_PER_HR.")
        return None

    print(f"{'ID':>8}  {'GPU':>20}  {'VRAM':>6}  {'$/hr':>6}  {'CUDA':>5}  {'DL Mbps':>8}")
    print("-" * 70)
    for r in results:
        print(f"{r.get('id','?'):>8}  {r.get('gpu_name','?'):>20}  "
              f"{r.get('gpu_ram',0):>5.0f}G  ${r.get('dph_total',0):>5.2f}  "
              f"{r.get('cuda_max_good','?'):>5}  {r.get('inet_down',0):>7.0f}")
    return results


def cmd_launch():
    print("Finding fastest available GPU...")
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {MIN_VRAM_GB}",
        f"dph <= {MAX_PRICE_PER_HR}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", "gpu_ram-",
        "--type", "on-demand",
        "--limit", "1",
        "--raw",
    )
    if not results:
        print("No suitable instances. Try 'search' to see what's available.")
        sys.exit(1)

    offer = results[0]
    gpu = offer.get("gpu_name", "?")
    vram = offer.get("gpu_ram", 0)
    price = offer.get("dph_total", 0)
    print(f"Launching: {gpu} ({vram:.0f}GB) at ${price:.2f}/hr")

    result = run_vastai(
        "create", "instance", str(offer["id"]),
        "--image", IMAGE,
        "--disk", str(DISK_GB),
        "--onstart-cmd", "apt-get update -qq && apt-get install -y -qq wget git > /dev/null 2>&1",
        "--raw",
    )

    instance_id = (result.get("new_contract") if isinstance(result, dict) else str(result))
    print(f"Instance {instance_id} created. Waiting for startup...")

    for i in range(90):
        time.sleep(5)
        instances = run_vastai("show", "instances", "--raw")
        if isinstance(instances, list):
            for inst in instances:
                if str(inst.get("id")) == str(instance_id):
                    status = inst.get("actual_status", "unknown")
                    if status == "running":
                        info = {
                            "id": instance_id,
                            "gpu": gpu,
                            "vram_gb": vram,
                            "price_per_hr": price,
                            "ssh_host": inst.get("ssh_host", ""),
                            "ssh_port": inst.get("ssh_port", ""),
                            "launched_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                        }
                        save_instance(info)
                        print(f"\nRUNNING: {gpu} ({vram:.0f}GB) at ${price:.2f}/hr")
                        print(f"SSH: ssh -p {info['ssh_port']} root@{info['ssh_host']}")
                        return info
                    print(f"  [{(i+1)*5}s] {status}...")
                    break

    print("Timeout. Check 'vastai show instances'.")
    sys.exit(1)


def cmd_upload():
    info = load_instance()
    print(f"Uploading validation experiment to {info['gpu']}...")

    # Create remote dirs
    subprocess.run(f'{ssh_cmd(info)} "mkdir -p /workspace/validation"', shell=True, check=True)

    # Upload validation scripts + YAMLs
    subprocess.run(
        scp_cmd(info, f'"{VALIDATION_DIR}/validation_experiment.py"', "/workspace/validation/"),
        shell=True, check=True,
    )
    subprocess.run(
        scp_cmd(info, f'"{VALIDATION_DIR}/run_validation_vastai.sh"', "/workspace/validation/"),
        shell=True, check=True,
    )
    subprocess.run(
        scp_cmd(info, f'"{VALIDATION_DIR}/boltz2_inputs"', "/workspace/validation/"),
        shell=True, check=True,
    )
    print("Upload complete.")


def cmd_run():
    info = load_instance()
    print(f"Running validation on {info['gpu']} ({info['vram_gb']:.0f}GB)...")
    print("This will take ~1-2 hours on H100, ~2-3 hours on A100.")
    print()

    # Run the experiment script
    subprocess.run(
        f'{ssh_cmd(info)} "cd /workspace/validation && bash run_validation_vastai.sh"',
        shell=True,
    )


def cmd_download():
    info = load_instance()
    local_out = os.path.join(VALIDATION_DIR, "boltz2_out")
    os.makedirs(local_out, exist_ok=True)

    print(f"Downloading results from {info['gpu']}...")
    subprocess.run(
        f'scp -o StrictHostKeyChecking=no -r -P {info["ssh_port"]} '
        f'root@{info["ssh_host"]}:/workspace/validation/boltz2_out/* "{local_out}/"',
        shell=True,
    )

    # Also grab the analysis JSON
    subprocess.run(
        f'scp -o StrictHostKeyChecking=no -P {info["ssh_port"]} '
        f'root@{info["ssh_host"]}:/workspace/validation/boltz2_out/validation_analysis.json '
        f'"{local_out}/"',
        shell=True,
    )
    print(f"Results saved to {local_out}/")


def cmd_destroy():
    info = load_instance()
    runtime_min = (time.time() - time.mktime(time.strptime(info["launched_at"], "%Y-%m-%d %H:%M:%S"))) / 60
    est_cost = (runtime_min / 60) * info["price_per_hr"]
    print(f"Destroying {info['gpu']} (ran {runtime_min:.0f} min, ~${est_cost:.2f})...")

    run_vastai("destroy", "instance", str(info["id"]), parse_json=False)
    os.remove(INSTANCE_FILE)
    print("Instance destroyed.")


def cmd_full():
    """Full pipeline: search → launch → upload → run → download → destroy."""
    results = cmd_search()
    if not results:
        return
    print()

    info = cmd_launch()
    print()

    # Brief pause for SSH to be ready
    print("Waiting 15s for SSH to stabilize...")
    time.sleep(15)

    cmd_upload()
    print()
    cmd_run()
    print()
    cmd_download()
    print()

    print("Experiment complete. Destroy instance? [y/N] ", end="", flush=True)
    if input().strip().lower() == "y":
        cmd_destroy()
    else:
        print(f"Instance still running at ${info['price_per_hr']:.2f}/hr.")
        print(f"Destroy later: python {__file__} destroy")


def main():
    commands = {
        "search": cmd_search,
        "launch": cmd_launch,
        "upload": cmd_upload,
        "run": cmd_run,
        "download": cmd_download,
        "destroy": cmd_destroy,
    }

    if len(sys.argv) < 2:
        cmd_full()
    elif sys.argv[1] in commands:
        commands[sys.argv[1]]()
    else:
        print(f"Unknown: {sys.argv[1]}")
        print(f"Available: {', '.join(commands)} (or no args for full pipeline)")
        sys.exit(1)


if __name__ == "__main__":
    main()
