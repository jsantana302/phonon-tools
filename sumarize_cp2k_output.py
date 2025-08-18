#!/user/j.santanaandreo/u12658/miniconda3/envs/aim2dat/bin/python
import sys
import subprocess
import re

try:
    import yaml
except ImportError:
    sys.stderr.write("ERROR: PyYAML is required (pip install pyyaml)\n")
    sys.exit(1)

def cp2kparse_yaml(fn):
    """Run cp2kparse in YAML mode and return a Python dict."""
    cmd = ["cp2kparse", "--format", "yaml", fn]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"ERROR: cp2kparse failed:\n{proc.stderr}")
        sys.exit(proc.returncode)

    raw = proc.stdout
    try:
        return yaml.safe_load(raw)
    except Exception:
        # Strip out non-printable control chars, retry
        cleaned = re.sub(r'[\x00-\x08\x0B-\x0C\x0E-\x1F]+', '', raw)
        return yaml.safe_load(cleaned)

def count_warnings(fn):
    cnt = 0
    with open(fn, errors="ignore") as f:
        for L in f:
            if "WARNING" in L:
                cnt += 1
    return cnt

def extract_energy(data):
    # cp2kparse YAML puts energies under top-level 'energies'
    e = data.get("energies", {}) or {}
    # Prefer the "total force_eval" key if present
    if "total force_eval" in e:
        return e["total force_eval"]
    # else first available
    for v in e.values():
        return v
    return None

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <CP2K_output_file>")
        sys.exit(1)
    fn = sys.argv[1]
    data = cp2kparse_yaml(fn)

    summary = {
        "final_energy":   extract_energy(data),
        "run_type":       data.get("global", {}).get("run type"),
        "version":        data.get("cp2k", {}).get("version string"),
        "warnings_count": count_warnings(fn),
    }

    import json
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
