path = "../data-acquisition/GSE54503_series_matrix.txt"

with open(path, 'r', errors='ignore') as f:
    lines = f.readlines()

print(f"Total lines: {len(lines)}")
print("\n--- First 15 lines (metadata) ---")
for l in lines[:15]:
    print(l.strip()[:150])

data_start = None
for i, l in enumerate(lines):
    if l.startswith('!series_matrix_table_begin') or l.startswith('"ID_REF"'):
        data_start = i
        break

print(f"\nData table appears to start around line: {data_start}")
if data_start is not None:
    print("\n--- Header + first few data rows ---")
    for l in lines[data_start:data_start+5]:
        print(l.strip()[:200])
