"""Validate accepted records without executing a trajectory, QR or eigensolver."""
import math
from pathlib import Path
from mpmath import mp
from common import ROOT, configs, read, require, frozen_check, digest, yacht_units


def load_records(directory=None):
    base = Path(directory) if directory else ROOT / "data/accepted"
    scfg, ycfg = configs()
    synthetic = [read(base / "synthetic" / (u["id"] + ".json")) for u in scfg["units"]]
    yacht = [read(base / "yacht" / (u["id"] + ".json")) for u in yacht_units(ycfg)]
    setups = {d: read(base / "yacht" / f"setup_{d}.json") for d in ycfg["precisions"]}
    return synthetic, yacht, setups


def validate(directory=None):
    count = frozen_check()
    scfg, ycfg = configs()
    srows, yrows, setups = load_records(directory)
    require(scfg["horizon"] == ycfg["horizon"] == 512, "Frozen horizon changed")
    require(scfg["precisions"] == ycfg["precisions"] == [80, 120], "Frozen precisions changed")
    require(scfg["seed_policy"] == "fully_deterministic_no_rng_no_seed"
            and ycfg["seed_policy"] == "deterministic_no_rng_no_seed", "Seed policy changed")
    require(len(srows) == 13 and len(yrows) == 16, "Missing work units")
    require(digest(ROOT / "data/raw/yacht_hydrodynamics.data") == ycfg["raw_sha256"], "Raw data changed")
    raw = [line.split() for line in (ROOT / "data/raw/yacht_hydrodynamics.data").read_text().splitlines() if line.strip()]
    require(len(raw) == 308 and all(len(row) == 7 for row in raw), "Unexpected raw data dimensions")
    precision, matrix, oracle = mp.mpf(0), mp.mpf(0), mp.mpf(0)
    with mp.workdps(120):
        for unit, row in zip(scfg["units"], srows):
            require(row["unit"] == unit and row["accepted"] is True, "Invalid synthetic record")
            if unit["kind"] == "scale":
                r = row["results"]
                require(mp.mpf(row["precision_relative_error"]) <= mp.mpf("1e-50")
                        and mp.mpf(r["coefficient_identity_error"]) <= mp.mpf("1e-50"), "Scalar precision gate")
                checks = r["rate_checks"]
                require(len(checks) == 3 and max(mp.mpf(v["agreement"]) for v in checks) <= mp.mpf("1e-50"), "Scalar rate gate")
                require(mp.mpf(checks[0]["q_difference"]) < 0 < mp.mpf(checks[2]["q_difference"])
                        and abs(mp.mpf(checks[1]["q_difference"])) <= mp.mpf("1e-50"), "Scalar threshold gate")
            else:
                require(set(row["methods"]) == set(scfg["methods"]), "Missing synthetic method")
                for r in row["methods"].values():
                    require(r["accepted"] is True and mp.mpf(r["precision_max_log_error"]) <= mp.mpf("1e-45"), "Synthetic precision gate")
                    require(mp.mpf(r["max_relative_error_violation"]) <= mp.mpf("1e-60")
                            and max(mp.mpf(v) for v in r["matrix_relative_errors"].values()) <= mp.mpf("1e-60"), "Synthetic acceptance gate")
                    require(set(r["checkpoints"]) == set(map(str, scfg["checkpoints"]))
                            and all(len(v) == 8 for v in r["checkpoints"].values()), "Synthetic checkpoints")
                    require(len(r["trace"]) == 513 and [v["k"] for v in r["trace"]] == list(range(513)), "Synthetic trace dimensions")
                    for v in r["trace"]:
                        require(all(math.isfinite(v[k]) for k in ["ln_relative_x", "ln_relative_state", "log10_x_over_qk"]), "Nonfinite synthetic trace")
        keyed = {}
        for unit, row in zip(yacht_units(ycfg), yrows):
            require(row["unit"] == unit and row["accepted"] is True
                    and row["gradient_updates"] == 512, "Invalid Yacht record")
            require(len(row["trace"]) == 513 and [v["k"] for v in row["trace"]] == list(range(513)), "Yacht trace dimensions")
            require(set(row["checkpoints"]) == set(map(str, ycfg["checkpoints"]))
                    and all(len(v) == 12 for v in row["checkpoints"].values()), "Yacht checkpoints")
            require(all(mp.isfinite(mp.mpf(v[k])) for v in row["trace"]
                        for k in ["ln_relative_x", "ln_relative_f", "ln_relative_state"]), "Nonfinite Yacht trace")
            matrix = max(matrix, max(mp.mpf(v) for v in row["matrix_relative_errors"].values()))
            oracle = max(oracle, mp.mpf(row["max_oracle_violation"]))
            keyed[(unit["profile"], unit["method"], unit["digits"])] = row
        for profile in ycfg["profiles"]:
            for method in ycfg["methods"]:
                lo, hi = [keyed[(profile, method, d)] for d in ycfg["precisions"]]
                precision = max(precision, max(abs(mp.mpf(a[k]) - mp.mpf(b[k]))
                    for a, b in zip(lo["trace"], hi["trace"])
                    for k in ["ln_relative_x", "ln_relative_f", "ln_relative_state"]))
        require(precision <= mp.mpf(ycfg["tolerances"]["precision_log"])
                and matrix <= mp.mpf(ycfg["tolerances"]["matrix_relative"])
                and oracle <= mp.mpf(ycfg["tolerances"]["oracle_relative"]), "Yacht acceptance gates")
        for setup in setups.values():
            require(len(setup["H"]) == 6 and all(len(r) == 6 for r in setup["H"]), "Hessian dimensions")
            require(mp.mpf(setup["normalized_reference_gradient"]) <= mp.mpf("1e-60"), "Stored reference residual")
            p = setup["p"]
            require(all(mp.mpf(p["mu"]) <= mp.mpf(v) <= mp.mpf(p["L"])
                        and (mp.mpf(v) <= mp.mpf(p["a"]) or mp.mpf(v) >= mp.mpf(p["b"]))
                        for v in setup["spectrum"]), "Stored spectral inclusion")
        result = dict(accepted=True, frozen_files=count, synthetic_units=13, synthetic_trajectory_arrays=16,
                      yacht_records=16, rows_per_trajectory=513, raw_dimensions=[308, 7],
                      yacht_max_precision_log_error=mp.nstr(precision, 120),
                      yacht_max_matrix_relative_error=mp.nstr(matrix, 120),
                      yacht_max_oracle_relative_violation=mp.nstr(oracle, 120),
                      new_gradient_updates=0, new_QR_solves=0, new_spectral_calculations=0)
    return result
