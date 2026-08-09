#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
nuclei_seg_gpu2.py (Batch-Ready Version with --outdir)
------------------------------------------------------

- AUTO-selects DAPI channel
- Napari disabled
- Batch-safe (no interactive input)
- Supports optional:  --outdir PATH
"""

import sys, os, re, time, argparse, xml.etree.ElementTree as ET
import numpy as np
import tifffile as tiff
import czifile
from scipy.ndimage import gaussian_filter, zoom
from skimage.transform import resize
from cellpose import models

HAVE_NAPARI = False  # Disabled for batch

# ============================================================
# Metadata helpers (unchanged)
# ============================================================

def _strip_ns(xml_text: str) -> str:
    return re.sub(r'xmlns="[^"]+"', '', xml_text)

def _unit_to_um(u: str | None) -> float:
    if not u:
        return 1e6
    u = u.strip().lower()
    if u in ("m", "meter", "meters"): return 1e6
    if u in ("um", "µm", "micrometer", "micrometers", "micron", "microns"): return 1.0
    if u in ("nm", "nanometer", "nanometers"): return 1e-3
    return 1e6

def vox_from_czi_xml(path):
    try:
        with czifile.CziFile(path) as czi:
            md = czi.metadata()
        if not md:
            return None
        root = ET.fromstring(_strip_ns(md))

        def get_axis_um(axis_id):
            for dist in root.findall(".//Distance"):
                aid = (dist.findtext("Id") or dist.get("Id") or "").strip().upper()
                if aid != axis_id:
                    continue
                val  = (dist.findtext("Value") or dist.get("Value"))
                unit = (dist.findtext("DefaultUnit") or dist.get("DefaultUnit"))
                if val:
                    return float(val) * _unit_to_um(unit)
            node = root.find(f".//Scaling/Items/Distance[@Id='{axis_id}']")
            if node is not None and node.get("Value"):
                return float(node.get("Value")) * _unit_to_um(node.get("DefaultUnit"))
            return None

        sx = get_axis_um("X")
        sy = get_axis_um("Y")
        sz = get_axis_um("Z")
        if sx and sy and sz:
            return (sz, sy, sx)
    except Exception:
        pass
    return None

def vox_from_ome_tiff(path):
    try:
        with tiff.TiffFile(path) as tf:
            ome = tf.ome_metadata
        if not ome:
            return None
        root = ET.fromstring(_strip_ns(ome))
        pix = root.find(".//Pixels")
        if pix is None:
            return None
        def getf(attr):
            v = pix.get(attr)
            return float(v) if v is not None else None
        sx = getf("PhysicalSizeX")
        sy = getf("PhysicalSizeY")
        sz = getf("PhysicalSizeZ")
        if sx and sy and sz:
            return (sz, sy, sx)
    except Exception:
        pass
    return None

def get_voxel_sizes_um(path):
    if path.lower().endswith(".czi"):
        vox = vox_from_czi_xml(path)
        if vox: return vox
    vox = vox_from_ome_tiff(path)
    if vox: return vox
    raise RuntimeError("Could not find voxel sizes in metadata.")

def chnames_from_czi_xml(path):
    try:
        with czifile.CziFile(path) as czi:
            md = czi.metadata()
        if not md:
            return None
        root = ET.fromstring(_strip_ns(md))
        names = []
        for ch in root.findall(".//Channels/Channel"):
            nm = (ch.findtext("Name") or ch.findtext("ShortName") or "").strip()
            if nm: names.append(nm)
        return names or None
    except Exception:
        return None

def chnames_from_ome_tiff(path):
    try:
        with tiff.TiffFile(path) as tf:
            ome = tf.ome_metadata
        if not ome:
            return None
        root = ET.fromstring(_strip_ns(ome))
        names = []
        for ch in root.findall(".//Channel"):
            nm = (ch.get("Name") or "").strip()
            if nm: names.append(nm)
        return names or None
    except Exception:
        return None

def get_channel_names(path):
    if path.lower().endswith(".czi"):
        n = chnames_from_czi_xml(path)
        if n: return n
    n = chnames_from_ome_tiff(path)
    return n

# ============================================================
# Shape helper
# ============================================================

def to_ZYXC(a):
    a = np.squeeze(a)
    if a.ndim == 5:
        a = a[0]
    if a.ndim == 4 and a.shape[0] <= 6:
        a = np.moveaxis(a, 0, -1)
    elif a.ndim == 3:
        a = a[..., np.newaxis]
    return a

# ============================================================
# Main
# ============================================================

def main(in_path, outdir):
    if not os.path.isfile(in_path):
        raise FileNotFoundError(in_path)

    name = os.path.splitext(os.path.basename(in_path))[0]

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # voxel sizes
    try:
        sz, sy, sx = get_voxel_sizes_um(in_path)
        print(f"Voxel sizes (µm): Z={sz:.6f}, Y={sy:.6f}, X={sx:.6f}")
    except Exception as e:
        print("[WARN] Could not parse voxel sizes:", e)

    names = get_channel_names(in_path)

    # load image
    arr = czifile.imread(in_path) if in_path.lower().endswith(".czi") else tiff.imread(in_path)
    arr = to_ZYXC(arr)
    print("Pixel array shaped to:", arr.shape, "(Z,Y,X,C)")

    nC = arr.shape[-1]
    print("\nAvailable channels:")
    for i in range(nC):
        nm = names[i] if names and i < len(names) else f"Ch{i}"
        print(f" [{i}] {nm}")

    # ============================================================
    # AUTO-PICK DAPI
    # ============================================================
    if names:
        nuc_idx = None
        for i, nm in enumerate(names):
            if "DAPI" in nm.upper():
                nuc_idx = i
                break
        if nuc_idx is None:
            print("[WARN] DAPI not found → using last channel")
            nuc_idx = nC - 1
    else:
        print("[WARN] No channel names → using last channel")
        nuc_idx = nC - 1

    print(f"Using channel index {nuc_idx} for nuclei.")

    blue = arr[..., nuc_idx]

    # smooth
    blue_smooth = np.zeros(blue.shape, dtype=np.float32)
    for z in range(blue.shape[0]):
        blue_smooth[z,:,:] = gaussian_filter(blue[z,:,:], 2)
    blue_smooth = (blue_smooth / np.max(blue_smooth) * 2**16).astype('uint16')

    # resize
    S = 1024
    blue_resized = resize(blue_smooth, (12, S, S), anti_aliasing=True)
    print("Resized volume:", blue_resized.shape)

    # Cellpose
    nuclei_model = models.CellposeModel(gpu=True)
    t1 = time.time()
    mask_blue,_,_ = nuclei_model.eval(
        blue_resized,
        diameter=30,
        flow_threshold=0.4,
        do_3D=True,
        z_axis=0
    )
    t2 = time.time()
    print("Cellpose done in %.2f minutes" % ((t2 - t1) / 60))

    # native grid
    Z_native, Y_native, X_native = blue.shape
    mask_native = zoom(mask_blue.astype(np.int32),
                       (Z_native/12, Y_native/1024, X_native/1024),
                       order=0)

    print("Mask rescaled:", mask_native.shape)

    # save to OUTDIR
    out_path = os.path.join(outdir, f"{name}_mask.tif")
    tiff.imwrite(out_path, mask_native.astype(np.int32))
    print("Saved mask:", out_path)

    print("Napari disabled. Continuing batch...")

# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input CZI or TIFF file")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    main(args.input, args.outdir)
