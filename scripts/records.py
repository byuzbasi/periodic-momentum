"""Check saved synthetic records without numerical trajectories or solvers."""
import math
from pathlib import Path
from mpmath import mp
from common import ROOT, configs, read, require, frozen_check

def load_records(directory=None):
    base=Path(directory) if directory else ROOT/'data/accepted'
    return [read(base/'synthetic'/(u['id']+'.json')) for u in configs()['units']]

def validate(directory=None):
    count=frozen_check();scfg=configs();srows=load_records(directory)
    require(len(srows)==13 and scfg['horizon']==512 and scfg['precisions']==[80,120], 'Synthetic design changed')
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
    return dict(accepted=True,synthetic_units=13,synthetic_trajectory_arrays=16,
        frozen_files=count,new_gradient_updates=0,new_QR_solves=0,new_spectral_calculations=0)
