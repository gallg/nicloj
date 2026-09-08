#!/usr/bin/env python3
"""Build the NIfTI test corpus and the manifest of what nibabel reads from it.

Run with the project's uv venv active:

    python scripts/gen_testdata.py

Writes files into test-data/ plus test-data/manifest.edn, which the Clojure
test suite reads to check nicloj against nibabel field by field.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
import struct
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "test-data"
NILEARN_DATA = Path(nib.__file__).parent.parent / "nilearn" / "datasets" / "data"

# Sampled voxel positions are taken on this stride so the manifest stays small.
SAMPLE_COUNT = 24


# --------------------------------------------------------------------- EDN out

def edn(value, indent=0):
    pad = "  " * indent
    if isinstance(value, dict):
        items = "\n".join(
            f"{pad}  {edn(k)} {edn(v, indent + 1).lstrip()}" for k, v in value.items()
        )
        return f"{pad}{{\n{items}}}"
    if isinstance(value, (list, tuple)):
        inner = " ".join(edn(v).strip() for v in value)
        return f"{pad}[{inner}]"
    if isinstance(value, str):
        return pad + '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(value, Keyword):
        return pad + ":" + value.name
    if value is None:
        return pad + "nil"
    if isinstance(value, (bool, np.bool_)):
        return pad + ("true" if value else "false")
    if isinstance(value, (int, np.integer)):
        return pad + str(int(value))
    if isinstance(value, (float, np.floating)):
        v = float(value)
        if np.isnan(v):
            return pad + "##NaN"
        if np.isinf(v):
            return pad + ("##Inf" if v > 0 else "##-Inf")
        return pad + repr(v)
    raise TypeError(f"cannot write {type(value)} as EDN")


class Keyword:
    __slots__ = ("name",)

    def __init__(self, name):
        self.name = name


def kw(name):
    return Keyword(str(name).replace(" ", "-").replace("_", "-"))


# ------------------------------------------------------------------ inspection

def sample_indices(shape):
    """A deterministic spread of voxel indices across `shape`."""
    n = int(np.prod(shape))
    step = max(1, n // SAMPLE_COUNT)
    return [np.unravel_index(i, shape, order="F") for i in range(0, n, step)][:SAMPLE_COUNT]


def stats(data):
    finite = data[np.isfinite(data)]
    return {
        kw("min"): float(finite.min()) if finite.size else 0.0,
        kw("max"): float(finite.max()) if finite.size else 0.0,
        kw("sum"): float(finite.sum()),
        kw("mean"): float(finite.mean()) if finite.size else 0.0,
    }


def disk_fields(path):
    """Read the fields nibabel normalises away once a file is loaded.

    On load, nibabel moves scl_slope/scl_inter into its array proxy and blanks
    both them and vox_offset in the header it hands back, so only the raw bytes
    say what a file really holds.
    """
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    end = "<" if struct.unpack_from("<i", raw, 0)[0] in (348, 540) else ">"
    if struct.unpack_from(end + "i", raw, 0)[0] == 540:
        return {
            kw("vox-offset"): struct.unpack_from(end + "q", raw, 168)[0],
            kw("scl"): list(struct.unpack_from(end + "dd", raw, 176)),
        }
    return {
        kw("vox-offset"): int(struct.unpack_from(end + "f", raw, 108)[0]),
        kw("scl"): [float(v) for v in struct.unpack_from(end + "ff", raw, 112)],
    }


def geometry(img, path):
    hdr = img.header
    data = img.get_fdata(dtype=np.float64)
    ornt = nib.orientations.io_orientation(img.affine)
    canon = nib.as_closest_canonical(img)
    canon_data = canon.get_fdata(dtype=np.float64)
    xyz_unit, t_unit = hdr.get_xyzt_units()
    return {
        kw("version"): 2 if hdr.sizeof_hdr == 540 else 1,
        kw("byte-order"): kw("big" if hdr.endianness == ">" else "little"),
        kw("single-file?"): bool(hdr["magic"].item().startswith(b"n+")),
        kw("shape"): list(img.shape),
        kw("dtype"): kw(np.dtype(hdr.get_data_dtype()).name),
        kw("zooms"): [float(z) for z in hdr.get_zooms()],
        kw("affine"): [[float(v) for v in row] for row in img.affine],
        kw("qform-code"): kw(hdr.get_value_label("qform_code")),
        kw("sform-code"): kw(hdr.get_value_label("sform_code")),
        kw("qform"): [[float(v) for v in row] for row in hdr.get_qform()],
        kw("sform"): [[float(v) for v in row] for row in hdr.get_sform()],
        kw("axcodes"): list(nib.aff2axcodes(img.affine)),
        **disk_fields(path),
        kw("xyzt-units"): [kw(xyz_unit), kw(t_unit)],
        kw("dim-info"): list(hdr.get_dim_info()),
        kw("descrip"): hdr["descrip"].item().decode("latin-1"),
        kw("intent-code"): kw(hdr.get_value_label("intent_code")),
        kw("slice-code"): kw(hdr.get_value_label("slice_code")),
        kw("n-extensions"): len(hdr.extensions),
        kw("stats"): stats(data),
        kw("samples"): [
            [list(int(i) for i in idx), float(data[idx])] for idx in sample_indices(img.shape)
        ],
        kw("ornt"): [
            [None if np.isnan(v) else int(v) for v in row] for row in ornt
        ],
        kw("canonical"): {
            kw("shape"): list(canon.shape),
            kw("affine"): [[float(v) for v in row] for row in canon.affine],
            kw("axcodes"): list(nib.aff2axcodes(canon.affine)),
            kw("stats"): stats(canon_data),
            kw("samples"): [
                [list(int(i) for i in idx), float(canon_data[idx])]
                for idx in sample_indices(canon.shape)
            ],
        },
    }


# ------------------------------------------------------------------- fixtures

def ramp(shape, dtype, lo=0, hi=100):
    """A reproducible non-trivial volume spanning [lo, hi].

    The seed comes from sha256 rather than hash(), which Python salts per
    process -- that made the corpus, and therefore the whole suite, differ from
    run to run.
    """
    key = f"{tuple(shape)}|{np.dtype(dtype).name}|{lo}|{hi}".encode()
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(key).digest()[:8], "big"))
    base = np.linspace(lo, hi, int(np.prod(shape))).reshape(shape, order="F")
    out = base + rng.normal(0, (hi - lo) / 50, shape)
    if np.dtype(dtype).kind in "iu":
        # Clip before casting: noise around a ramp starting at 0 can go
        # negative, and an unsigned cast would wrap it to near the type maximum
        # instead of keeping the volume inside [lo, hi] as intended.
        info = np.iinfo(dtype)
        out = np.clip(np.rint(out), info.min, info.max)
    return np.ascontiguousarray(out.astype(dtype))


def oblique_affine():
    """A rotated, non-axis-aligned affine with a translation."""
    a, b = np.deg2rad(25), np.deg2rad(-12)
    rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rot = rz @ ry @ np.diag([1.5, 1.5, 3.0])
    return nib.affines.from_matvec(rot, [-31.5, 12.25, -8.0])


def build():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    plain = np.diag([2.0, 2.0, 2.5, 1.0])
    plain[:3, 3] = [-10.0, -20.0, -15.0]
    written = []

    def write(name, img):
        path = OUT / name
        nib.save(img, path)
        written.append(name)

    # NIfTI-1, one file, the common float32 case.
    img = nib.Nifti1Image(ramp((7, 8, 9), np.float32), plain)
    img.header.set_xyzt_units("mm", "sec")
    img.header.set_dim_info(freq=0, phase=1, slice=2)
    img.header["descrip"] = b"nicloj float32 fixture"
    write("nifti1-3d-float32.nii", img)

    # Gzipped int16 with explicit scaling, the usual scanner output.
    hdr = nib.Nifti1Header()
    hdr.set_data_dtype(np.int16)
    img = nib.Nifti1Image(ramp((6, 7, 5), np.int16), plain, hdr)
    img.header.set_slope_inter(0.25, -3.5)
    write("nifti1-3d-int16-scaled.nii.gz", img)

    # 4D time series.
    write("nifti1-4d-float64.nii.gz",
          nib.Nifti1Image(ramp((5, 6, 4, 3), np.float64), plain))

    # 5D, to exercise dim[5].
    write("nifti1-5d-float32.nii.gz",
          nib.Nifti1Image(ramp((4, 4, 3, 2, 2), np.float32), plain))

    # hdr/img pair, plain and gzipped.
    write("nifti1-pair.hdr", nib.Nifti1Pair(ramp((5, 5, 5), np.int16), plain))
    write("nifti1-pair-gz.hdr.gz", nib.Nifti1Pair(ramp((4, 5, 6), np.float32), plain))

    # Big-endian.
    be = nib.Nifti1Image(ramp((5, 4, 3), np.int16), plain,
                         nib.Nifti1Header(endianness=">"))
    write("nifti1-3d-int16-bigendian.nii", be)

    # NIfTI-2.
    write("nifti2-3d-float32.nii", nib.Nifti2Image(ramp((6, 6, 6), np.float32), plain))
    write("nifti2-4d-int32.nii.gz", nib.Nifti2Image(ramp((4, 5, 6, 2), np.int32), plain))

    # Every real datatype nicloj supports.
    for name in ["uint8", "int8", "int16", "uint16", "int32", "uint32",
                 "int64", "uint64", "float32", "float64"]:
        dt = np.dtype(name)
        lo, hi = (0, 200) if dt.kind == "u" else (-100, 100)
        hdr = nib.Nifti1Header()
        hdr.set_data_dtype(dt)
        write(f"dtype-{name}.nii",
              nib.Nifti1Image(ramp((4, 5, 3), dt, lo, hi), plain, hdr))

    # Oblique affine, both forms set.
    write("oblique.nii.gz", nib.Nifti1Image(ramp((6, 7, 8), np.float32), oblique_affine()))

    # qform only, sform only, and neither.
    base = ramp((5, 6, 7), np.float32)
    img = nib.Nifti1Image(base, plain)
    img.header.set_sform(None, code="unknown")
    img.header.set_qform(plain, code="scanner")
    write("qform-only.nii", nib.Nifti1Image(base, None, img.header))

    img = nib.Nifti1Image(base, plain)
    img.header.set_qform(None, code="unknown")
    write("sform-only.nii", nib.Nifti1Image(base, None, img.header))

    img = nib.Nifti1Image(base, plain)
    img.header.set_sform(None, code="unknown")
    img.header.set_qform(None, code="unknown")
    write("no-form.nii", nib.Nifti1Image(base, None, img.header))

    # Non-RAS orientations, for reorientation tests.
    for codes in [("L", "I", "A"), ("P", "S", "R"), ("A", "L", "S")]:
        ornt = nib.orientations.axcodes2ornt(codes)
        affine = plain @ nib.orientations.inv_ornt_aff(ornt, (6, 7, 8))
        data = nib.orientations.apply_orientation(ramp((6, 7, 8), np.float32), ornt)
        write(f"ornt-{''.join(codes).lower()}.nii.gz",
              nib.Nifti1Image(np.ascontiguousarray(data), affine))

    # A header extension, to check it survives a read.
    img = nib.Nifti1Image(ramp((4, 4, 4), np.float32), plain)
    img.header.extensions.append(
        nib.nifti1.Nifti1Extension(6, b"nicloj extension payload"))
    write("with-extension.nii", img)

    # Real-world images shipped with nilearn.
    for src, name in [("mni_icbm152_t1_tal_nlin_sym_09a_converted.nii.gz", "mni152-t1.nii.gz"),
                      ("image_10426.nii.gz", "nilearn-image-10426.nii.gz")]:
        path = NILEARN_DATA / src
        if path.exists():
            shutil.copy(path, OUT / name)
            written.append(name)

    manifest = [{kw("file"): name, **geometry(nib.load(OUT / name), OUT / name)}
                for name in written]
    (OUT / "manifest.edn").write_text(edn(manifest) + "\n")
    print(f"wrote {len(written)} images and manifest.edn into {OUT}")


if __name__ == "__main__":
    build()
