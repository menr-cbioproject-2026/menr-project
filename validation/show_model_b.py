with open('../model-b-kinetics/model_b.py', 'r') as f:
    lines = f.readlines()

for i, line in enumerate(lines[23:50], start=23):
    print(f"{i}: {line}", end='')
