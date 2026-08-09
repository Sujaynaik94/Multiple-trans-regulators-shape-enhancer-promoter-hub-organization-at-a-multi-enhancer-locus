import os
import glob
import subprocess

# ✅ Use Python from your virtual environment
PY = r"C:\RSFISH\.venv\Scripts\python.exe"

# Your segmentation script name
SCRIPT = "nuclei_segmentation.py"

# Folder with your CZI files
DATA_DIR = r"C:\Users\pregerlab\Documents\Srijani\AP"

# Where to save the masks
OUT_DIR  = r"C:\Users\pregerlab\Documents\Srijani\AP\MASKS"

os.makedirs(OUT_DIR, exist_ok=True)

files = glob.glob(os.path.join(DATA_DIR, "*.czi"))

print(f"Found {len(files)} files.")

for f in files:
    print("\n---------------------------------------")
    print(f"Processing: {f}")
    print("---------------------------------------")

    # 👇 Force the venv python explicitly
    subprocess.run([
        PY,
        SCRIPT,
        f,
        "--outdir", OUT_DIR
    ], check=True)
