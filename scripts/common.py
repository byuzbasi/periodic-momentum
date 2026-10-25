"""File handling and access to the unchanged numerical kernels."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value, exclusive=True):
    path = Path(path)
    require(not exclusive or not path.exists(), f"Refusing to overwrite {path.name}")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def output_path(value, create=True):
    path = Path(value)
    path = (ROOT / path).resolve() if not path.is_absolute() else path.resolve()
    if create:
        require(not path.exists(), "Choose a new output directory")
        path.mkdir(parents=True)
    return path


def kernels():
    sys.path.insert(0, str(ROOT / "code/periodic_momentum_numerics_v06/code"))
    sys.path.insert(0, str(ROOT / "code/periodic_momentum_yacht_v09"))
    import numerics
    import dense_math
    require(importlib.metadata.version("mpmath") == "1.3.0", "Use mpmath 1.3.0")
    return numerics, dense_math


def configs():
    return read(ROOT / "configs/synthetic.json"), read(ROOT / "configs/yacht.json")


def yacht_units(cfg):
    combinations = ((p, m, d) for p in cfg["profiles"]
                    for m in cfg["methods"] for d in cfg["precisions"])
    return [dict(id=f"Y{i+1:02d}", profile=p, method=m, digits=d)
            for i, (p, m, d) in enumerate(combinations)]


def environment():
    return {"python": platform.python_version(), "mpmath": importlib.metadata.version("mpmath")}


def frozen_check():
    rows = read(ROOT / "frozen_inputs.json")["files"]
    for relative, expected in rows.items():
        path = ROOT / relative
        require(path.is_file() and path.stat().st_size == expected["bytes"]
                and digest(path) == expected["sha256"], f"Frozen file mismatch: {relative}")
    return len(rows)


def signature():
    paths = [p for directory in ["configs", "code", "scripts"]
             for p in (ROOT / directory).rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    paths += [ROOT / "data/raw/yacht_hydrodynamics.data", ROOT / "frozen_inputs.json"]
    return {"environment": environment(),
            "files": {str(p.relative_to(ROOT)): digest(p) for p in sorted(paths)}}
