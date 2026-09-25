import re

with open('../model-b-kinetics/model_b.py', 'r') as f:
    content = f.read()

print("Current ODE section:")
for i, line in enumerate(content.split('\n')):
    if 'dnmt' in line.lower() or 'tet1' in line.lower() or 'dmdt' in line.lower():
        print(f"  Line {i}: {line}")
