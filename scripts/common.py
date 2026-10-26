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
    import numerics
    require(importlib.metadata.version("mpmath") == "1.3.0", "Use mpmath 1.3.0")
    return numerics, None


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
    paths += [ROOT / "frozen_inputs.json"]
    return {"environment": environment(),
            "files": {str(p.relative_to(ROOT)): digest(p) for p in sorted(paths)}}


def configs():
    return read(ROOT / "configs/synthetic.json")
