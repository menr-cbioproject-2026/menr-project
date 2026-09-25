import re
import glob

CANONICAL = {
    "PAAD_E0": 3.5, "PAAD_k": 0.8,
    "LIHC_E0": 5.5, "LIHC_k": 0.6,
    "TNBC_E0": 9.0, "TNBC_k": 0.35,
    "LUAD_E0": 5.5, "LUAD_k": 0.6,
    "kme": 0.15, "kde": 0.40,
    "time_window_hours": 48,
    "recovery_fraction": 0.40,
}

print("CANONICAL VALUES (from digital_twin.py):")
for k, v in CANONICAL.items():
    print(f"  {k}: {v}")

print("\n" + "="*60)
print("SCANNING ALL .py FILES")
print("="*60)

py_files = [f for f in glob.glob("**/*.py", recursive=True) if ".ipynb_checkpoints" not in f]

patterns = {
    "E0/k pairs": r'"E0":\s*([\d.]+),\s*"k":\s*([\d.]+)',
    "time window (0, N)": r't_span\s*=\s*\(0,\s*(\d+)\)|linspace\(0,\s*(\d+),|set_xlim\(0,\s*(\d+)\)',
    "recovery threshold": r'RECOVERY_FRACTION\s*=\s*([\d.]+)',
}

for f in sorted(py_files):
    with open(f, errors='ignore') as fh:
        content = fh.read()
    findings = []
    for label, pat in patterns.items():
        matches = re.findall(pat, content)
        if matches:
            findings.append((label, matches))
    if findings:
        print(f"\n{f}:")
        for label, matches in findings:
            print(f"  {label}: {matches}")
