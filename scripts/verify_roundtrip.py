#!/usr/bin/env python3
"""Re-read everything nicloj wrote and check it against the original with nibabel.

Run after the Clojure suite, which leaves its output plus written.edn in
test-data/out/:

    python scripts/verify_roundtrip.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "test-data"
OUT = CORPUS / "out"

# Fields a faithful rewrite has to preserve. vox_offset, scl_slope, scl_inter
# and magic legitimately change with the layout, so they are left out: the
# voxel comparison covers the scaling, and magic is checked against the variant.
HEADER_FIELDS = [
    "dim", "datatype", "bitpix", "qform_code", "sform_code",
    "quatern_b", "quatern_c", "quatern_d", "qoffset_x", "qoffset_y", "qoffset_z",
    "srow_x", "srow_y", "srow_z", "xyzt_units", "dim_info", "intent_code",
    "intent_name", "descrip", "aux_file", "slice_code", "cal_max", "cal_min",
    "toffset", "slice_duration",
]

ENTRY = re.compile(r"\{([^}]*)\}")
PAIR = re.compile(r':(?P<key>[\w-]+)\s+"(?P<value>[^"]*)"')


def read_written():
    """Pull the {:source :written :variant} maps out of the Clojure-written EDN."""
    path = OUT / "written.edn"
    if not path.exists():
        sys.exit(f"{path} missing; run the Clojure test suite first")
    entries = [
        {m.group("key"): m.group("value") for m in PAIR.finditer(block.group(1))}
        for block in ENTRY.finditer(path.read_text())
    ]
    missing = [e for e in entries if not {"source", "written", "variant"} <= e.keys()]
    if missing:
        sys.exit(f"{path} has {len(missing)} unparsable entries")
    return entries


def compare(source: Path, written: Path, variant: str, problems: list[str]):
    def bad(msg):
        problems.append(f"{written.name} ({variant}): {msg}")

    original = nib.load(source)
    copy = nib.load(written)

    if original.shape != copy.shape:
        bad(f"shape {original.shape} != {copy.shape}")
        return
    if original.header.get_data_dtype() != copy.header.get_data_dtype():
        bad(f"dtype {original.header.get_data_dtype()} != {copy.header.get_data_dtype()}")
    if not np.allclose(original.affine, copy.affine, rtol=0, atol=1e-5):
        bad(f"affine differs:\n{original.affine}\n{copy.affine}")
    if not np.allclose(original.header.get_zooms(), copy.header.get_zooms(), atol=1e-5):
        bad(f"zooms {original.header.get_zooms()} != {copy.header.get_zooms()}")

    want, got = original.get_fdata(dtype=np.float64), copy.get_fdata(dtype=np.float64)
    if not np.allclose(want, got, rtol=0, atol=0, equal_nan=True):
        worst = np.nanmax(np.abs(want - got))
        bad(f"voxels differ, largest gap {worst}")

    for field in HEADER_FIELDS:
        a, b = original.header[field], copy.header[field]
        if not np.array_equal(np.asarray(a).astype(np.float64) if a.dtype.kind == "f" else a,
                              np.asarray(b).astype(np.float64) if b.dtype.kind == "f" else b):
            bad(f"header field {field}: {a!r} != {b!r}")

    if len(original.header.extensions) != len(copy.header.extensions):
        bad(f"{len(original.header.extensions)} extensions became "
            f"{len(copy.header.extensions)}")

    expected_version = 2 if variant == "as-nifti2" else (
        2 if original.header.sizeof_hdr == 540 else 1)
    got_version = 2 if copy.header.sizeof_hdr == 540 else 1
    if expected_version != got_version:
        bad(f"expected NIfTI-{expected_version}, got NIfTI-{got_version}")

    expected_single = variant != "as-pair"
    got_single = copy.header["magic"].item().startswith(b"n+")
    if expected_single != got_single:
        bad(f"magic {copy.header['magic']!r} does not match variant")


def main():
    entries = read_written()
    problems: list[str] = []
    for entry in entries:
        compare(CORPUS / entry["source"], OUT / entry["written"], entry["variant"], problems)

    print(f"checked {len(entries)} files nicloj wrote against their originals")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems:
            print("  -", p)
        sys.exit(1)
    print("all match")


if __name__ == "__main__":
    main()
