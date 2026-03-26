"""Generate MSA for EvoBind using ColabFold's MMseqs2 API.

Replaces the 74GB Uniclust30 local database requirement.
ColabFold API is free and returns a3m-format MSAs.

Usage: python evobind_msa.py <input.fasta> <output.a3m>
"""
import sys
import os
import time
import requests
import json

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, Exception):
    pass

COLABFOLD_API = "https://api.colabfold.com"


def read_fasta(fasta_path):
    """Read first sequence from FASTA file."""
    header, seq_lines = "", []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    break
                header = line[1:]
            else:
                seq_lines.append(line)
    return header, "".join(seq_lines)


def submit_msa_job(sequence):
    """Submit sequence to ColabFold MMseqs2 API."""
    print("  Submitting to ColabFold MMseqs2 API...")
    resp = requests.post(
        f"{COLABFOLD_API}/ticket/msa",
        data={"q": f">query\n{sequence}", "mode": "all"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    ticket_id = data.get("id")
    if not ticket_id:
        raise RuntimeError(f"No ticket ID returned: {data}")
    print(f"  Ticket: {ticket_id}")
    return ticket_id


def poll_result(ticket_id, max_wait=600):
    """Poll for MSA result."""
    print("  Waiting for MSA generation...")
    start = time.time()
    while time.time() - start < max_wait:
        resp = requests.get(f"{COLABFOLD_API}/ticket/{ticket_id}", timeout=30)
        data = resp.json()
        status = data.get("status")

        if status == "COMPLETE":
            print("  MSA generation complete!")
            return data
        elif status == "ERROR":
            raise RuntimeError(f"MSA generation failed: {data}")
        elif status in ("PENDING", "RUNNING"):
            elapsed = int(time.time() - start)
            print(f"    Status: {status} ({elapsed}s elapsed)")
            time.sleep(10)
        else:
            print(f"    Unknown status: {status}")
            time.sleep(10)

    raise TimeoutError(f"MSA generation timed out after {max_wait}s")


def download_a3m(ticket_id, output_path):
    """Download the a3m MSA result."""
    resp = requests.get(f"{COLABFOLD_API}/result/download/{ticket_id}", timeout=60)
    resp.raise_for_status()

    # ColabFold returns a tar.gz with multiple files
    import tarfile
    import io

    tar_data = io.BytesIO(resp.content)
    with tarfile.open(fileobj=tar_data, mode="r:gz") as tar:
        # Find the a3m file
        a3m_content = None
        for member in tar.getmembers():
            if member.name.endswith(".a3m"):
                f = tar.extractfile(member)
                if f:
                    a3m_content = f.read().decode("utf-8")
                    break

        if not a3m_content:
            # Try uniref.a3m or paired.a3m
            for member in tar.getmembers():
                print(f"    Archive contains: {member.name}")
            # Fall back to first text file
            for member in tar.getmembers():
                f = tar.extractfile(member)
                if f:
                    content = f.read()
                    try:
                        a3m_content = content.decode("utf-8")
                        if a3m_content.startswith(">"):
                            break
                    except UnicodeDecodeError:
                        continue

    if not a3m_content:
        raise RuntimeError("No a3m content found in API response")

    with open(output_path, "w") as f:
        f.write(a3m_content)

    n_seqs = a3m_content.count("\n>")
    print(f"  Saved: {output_path} ({n_seqs} sequences)")
    return output_path


def main():
    if len(sys.argv) != 3:
        print("Usage: python evobind_msa.py <input.fasta> <output.a3m>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    output_path = sys.argv[2]

    header, sequence = read_fasta(fasta_path)
    print(f"  Sequence: {header}")
    print(f"  Length: {len(sequence)} aa")

    ticket_id = submit_msa_job(sequence)
    poll_result(ticket_id)
    download_a3m(ticket_id, output_path)


if __name__ == "__main__":
    main()
