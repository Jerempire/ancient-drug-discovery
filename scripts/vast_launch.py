"""
Vast.ai instance launcher for drug discovery GPU pipeline.

Usage:
    python scripts/vast_launch.py search       # Find available RTX 4090 / A100 instances
    python scripts/vast_launch.py launch       # Launch cheapest suitable instance
    python scripts/vast_launch.py status       # Check running instance
    python scripts/vast_launch.py ssh          # Print SSH command for running instance
    python scripts/vast_launch.py destroy      # Tear down running instance
    python scripts/vast_launch.py upload       # SCP data files to instance
    python scripts/vast_launch.py download     # SCP results from instance

Requires: pip install vastai
API key: set via `vastai set api-key <KEY>` or VAST_API_KEY env var
"""
import subprocess
import sys
import json
import os
import time

# --- Config ---
MIN_VRAM_GB = 20          # RTX 4090 = 24 GB, A100 = 40/80 GB
MAX_PRICE_PER_HR = 1.50   # USD
DISK_GB = 40              # Enough for model weights + results
IMAGE = "rosettacommons/rfdiffusion:latest"  # Pre-built RFdiffusion + DGL + weights
PREFERRED_GPUS = ["RTX_4090", "A100_SXM4", "A100_PCIE", "RTX_A6000", "L40S"]

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCE_FILE = os.path.join(PROJECT_DIR, ".vast_instance.json")


def run_vastai(*args, parse_json=True):
    """Run a vastai CLI command and return output."""
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
    """Save instance info to local file."""
    with open(INSTANCE_FILE, "w") as f:
        json.dump(info, f, indent=2)
    print(f"Instance info saved to {INSTANCE_FILE}")


def load_instance():
    """Load saved instance info."""
    if not os.path.exists(INSTANCE_FILE):
        print("No running instance found. Run 'vast_launch.py launch' first.")
        sys.exit(1)
    with open(INSTANCE_FILE) as f:
        return json.load(f)


def cmd_search():
    """Search for suitable GPU instances."""
    print(f"Searching for instances: >= {MIN_VRAM_GB} GB VRAM, <= ${MAX_PRICE_PER_HR}/hr...")
    print()

    # vastai search offers with filters
    # Query: GPU RAM >= MIN_VRAM_GB, price <= MAX_PRICE_PER_HR, CUDA >= 12.0
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {MIN_VRAM_GB}",
        f"dph <= {MAX_PRICE_PER_HR}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", "dph",
        "--type", "on-demand",
        "--limit", "15",
        "--raw",
    )

    if not results:
        print("No instances found matching criteria. Try increasing MAX_PRICE_PER_HR.")
        return

    print(f"{'ID':>8}  {'GPU':>16}  {'VRAM':>6}  {'$/hr':>6}  {'CUDA':>5}  {'DL Mbps':>8}  {'Reliability':>6}")
    print("-" * 75)
    for r in results:
        gpu_name = r.get("gpu_name", "?")
        gpu_ram = r.get("gpu_ram", 0)
        dph = r.get("dph_total", 0)
        cuda = r.get("cuda_max_good", "?")
        inet = r.get("inet_down", 0)
        rel = r.get("reliability2", 0)
        oid = r.get("id", "?")
        print(f"{oid:>8}  {gpu_name:>16}  {gpu_ram:>5.0f}G  ${dph:>5.3f}  {cuda:>5}  {inet:>7.0f}  {rel:>5.1%}")

    print(f"\nFound {len(results)} offers. Use 'vast_launch.py launch' to start cheapest.")


def cmd_launch():
    """Launch the cheapest suitable instance."""
    print("Finding best instance...")
    results = run_vastai(
        "search", "offers",
        f"gpu_ram >= {MIN_VRAM_GB}",
        f"dph <= {MAX_PRICE_PER_HR}",
        "cuda_vers >= 12.0",
        "reliability >= 0.95",
        "inet_down >= 200",
        "--order", "dph",
        "--type", "on-demand",
        "--limit", "1",
        "--raw",
    )

    if not results:
        print("No suitable instances available. Try again later or increase MAX_PRICE_PER_HR.")
        sys.exit(1)

    offer = results[0]
    offer_id = offer["id"]
    gpu = offer.get("gpu_name", "?")
    price = offer.get("dph_total", 0)

    print(f"Launching: {gpu} ({offer.get('gpu_ram', 0):.0f} GB) at ${price:.3f}/hr")
    print(f"Image: {IMAGE}")
    print(f"Disk: {DISK_GB} GB")

    # Create instance
    result = run_vastai(
        "create", "instance", str(offer_id),
        "--image", IMAGE,
        "--disk", str(DISK_GB),
        "--onstart-cmd", "apt-get update -qq && apt-get install -y -qq wget git > /dev/null 2>&1",
        "--raw",
    )

    if isinstance(result, dict) and "new_contract" in result:
        instance_id = result["new_contract"]
    elif isinstance(result, str) and "started" in result.lower():
        # Parse instance ID from text response
        instance_id = result.split()[-1] if result.split() else "unknown"
    else:
        instance_id = str(result)

    print(f"Instance created: {instance_id}")
    print("Waiting for instance to start...")

    # Poll until running
    for i in range(60):
        time.sleep(5)
        instances = run_vastai("show", "instances", "--raw")
        if isinstance(instances, list):
            for inst in instances:
                if str(inst.get("id")) == str(instance_id):
                    status = inst.get("actual_status", "unknown")
                    if status == "running":
                        ssh_host = inst.get("ssh_host", "")
                        ssh_port = inst.get("ssh_port", "")
                        print(f"\nInstance {instance_id} is RUNNING!")
                        print(f"SSH: ssh -p {ssh_port} root@{ssh_host}")

                        save_instance({
                            "id": instance_id,
                            "gpu": gpu,
                            "price_per_hr": price,
                            "ssh_host": ssh_host,
                            "ssh_port": ssh_port,
                            "launched_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                        })
                        return
                    else:
                        print(f"  Status: {status}... ({(i+1)*5}s)")
                        break

    print("Timeout waiting for instance. Check 'vastai show instances'.")


def cmd_status():
    """Show status of running instances."""
    instances = run_vastai("show", "instances", "--raw")
    if not instances:
        print("No running instances.")
        return

    for inst in instances:
        iid = inst.get("id", "?")
        gpu = inst.get("gpu_name", "?")
        status = inst.get("actual_status", "?")
        dph = inst.get("dph_total", 0)
        ssh_host = inst.get("ssh_host", "")
        ssh_port = inst.get("ssh_port", "")
        print(f"Instance {iid}: {gpu} | {status} | ${dph:.3f}/hr | ssh -p {ssh_port} root@{ssh_host}")


def cmd_ssh():
    """Print SSH command for running instance."""
    info = load_instance()
    print(f"ssh -p {info['ssh_port']} root@{info['ssh_host']}")


def cmd_destroy():
    """Destroy running instance."""
    info = load_instance()
    iid = info["id"]
    print(f"Destroying instance {iid} ({info['gpu']}, ${info['price_per_hr']:.3f}/hr)...")
    run_vastai("destroy", "instance", str(iid), parse_json=False)
    print("Instance destroyed.")

    if os.path.exists(INSTANCE_FILE):
        os.remove(INSTANCE_FILE)
        print("Cleaned up local instance file.")


def cmd_upload():
    """Upload structure + ligand files to the instance."""
    info = load_instance()
    host = info["ssh_host"]
    port = info["ssh_port"]

    data_dir = os.path.join(PROJECT_DIR, "data")
    structures = os.path.join(data_dir, "structures")
    ligands = os.path.join(data_dir, "ligands")

    print(f"Uploading data to {host}:{port}...")

    # Create workspace dirs on remote
    ssh_cmd = f"ssh -p {port} root@{host}"
    subprocess.run(f'{ssh_cmd} "mkdir -p /workspace/data/structures /workspace/data/ligands /workspace/scripts /workspace/results"', shell=True)

    # Upload structures
    for f in os.listdir(structures):
        if f.endswith(".pdb"):
            src = os.path.join(structures, f)
            subprocess.run(f'scp -P {port} "{src}" root@{host}:/workspace/data/structures/', shell=True)
            print(f"  Uploaded {f}")

    # Upload ligands
    for f in os.listdir(ligands):
        src = os.path.join(ligands, f)
        subprocess.run(f'scp -P {port} "{src}" root@{host}:/workspace/data/ligands/', shell=True)
        print(f"  Uploaded {f}")

    # Upload pipeline scripts
    scripts_dir = os.path.join(PROJECT_DIR, "scripts")
    for script in ["gpu_setup.sh", "run_pipeline_gpu.py"]:
        src = os.path.join(scripts_dir, script)
        if os.path.exists(src):
            subprocess.run(f'scp -P {port} "{src}" root@{host}:/workspace/scripts/', shell=True)
            print(f"  Uploaded {script}")

    print("Upload complete.")


def cmd_download():
    """Download results from instance."""
    info = load_instance()
    host = info["ssh_host"]
    port = info["ssh_port"]

    results_dir = os.path.join(PROJECT_DIR, "data", "results")
    os.makedirs(results_dir, exist_ok=True)

    print(f"Downloading results from {host}:{port}...")
    subprocess.run(
        f'scp -r -P {port} root@{host}:/workspace/results/* "{results_dir}/"',
        shell=True,
    )
    print(f"Results saved to {results_dir}/")

    # List what we got
    for f in os.listdir(results_dir):
        fpath = os.path.join(results_dir, f)
        size = os.path.getsize(fpath) if os.path.isfile(fpath) else "dir"
        print(f"  {f}: {size}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    commands = {
        "search": cmd_search,
        "launch": cmd_launch,
        "status": cmd_status,
        "ssh": cmd_ssh,
        "destroy": cmd_destroy,
        "upload": cmd_upload,
        "download": cmd_download,
    }

    cmd = sys.argv[1]
    if cmd not in commands:
        print(f"Unknown command: {cmd}")
        print(f"Available: {', '.join(commands)}")
        sys.exit(1)

    commands[cmd]()


if __name__ == "__main__":
    main()
