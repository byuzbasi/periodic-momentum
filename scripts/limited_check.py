"""Two prescribed prefixes: four state updates plus four gradient updates."""
import copy
from mpmath import mp
from common import ROOT, configs, kernels, read, require, write
from records import validate, load_records


def smoke(out):
    validate()
    policy = read(ROOT / "configs/validation_policy.json")
    require(policy == {"precision": 120, "prefix_steps": 2, "max_updates": 8,
                       "synthetic_unit": "T01", "methods": ["CHB", "HB"],
                       "yacht_units": ["Y02", "Y04"]}, "Validation budget changed")
    numeric, dense = kernels()
    scfg, ycfg = configs()
    srows, yrows, setups = load_records()
    reserved, rows = 0, []
    with mp.workdps(120):
        sample = next(r for r in srows if r["unit"]["id"] == "T01")
        p = numeric.parameters(scfg, sample["unit"]["R"])
        for method in policy["methods"]:
            reserved += 2
            require(reserved <= 8, "Smoke budget exceeded")
            write(out / "budget.json", {"reserved_updates": reserved, "maximum": 8}, exclusive=False)
            # Exercise the existing independent transition kernel against saved direct-loop states.
            x0 = [numeric.number(v) for v in scfg["initial_x"]]
            z = mp.matrix(x0 + x0)
            differences = []
            for k in range(2):
                z = numeric.transition(p, method, "exact", mp.mpf(0), k) * z
                expected = mp.matrix([mp.mpf(v) for v in sample["methods"][method]["checkpoints"][str(k+1)]])
                differences.append(numeric.norm(list(z-expected))/numeric.norm(list(expected)))
            rows.append(dict(example="synthetic/T01", method=method, state_updates=2, gradient_updates=0,
                             max_relative_checkpoint_error=mp.nstr(max(differences), 120),
                             accepted=bool(max(differences) <= mp.mpf(scfg["tolerances"]["matrix_relative"]))))
        prefix = copy.deepcopy(ycfg)
        prefix["horizon"], prefix["checkpoints"] = 2, [1,2]
        for identity in policy["yacht_units"]:
            expected = next(r for r in yrows if r["unit"]["id"] == identity)
            reserved += 2
            require(reserved <= 8, "Smoke budget exceeded")
            write(out / "budget.json", {"reserved_updates": reserved, "maximum": 8}, exclusive=False)
            counter = {"updates": 0}
            setup = dense.unpack_setup(setups[120])
            actual = dense.trajectory(prefix, setup, expected["unit"], counter)
            require(counter["updates"] == 2, "Unexpected smoke gradient count")
            differences = []
            for k in [1,2]:
                a = mp.matrix([mp.mpf(v) for v in actual["checkpoints"][str(k)]])
                b = mp.matrix([mp.mpf(v) for v in expected["checkpoints"][str(k)]])
                differences.append(numeric.norm(list(a-b))/numeric.norm(list(b)))
            rows.append(dict(example="yacht/"+identity, method=expected["unit"]["method"], state_updates=0,
                             gradient_updates=counter["updates"], max_relative_checkpoint_error=mp.nstr(max(differences),120),
                             accepted=bool(actual["accepted"] and max(differences) <= mp.mpf(ycfg["tolerances"]["matrix_relative"]))))
    result = dict(accepted=all(r["accepted"] for r in rows), rows=rows, total_updates=8,
                  independent_synthetic_state_updates=4, yacht_gradient_updates=4,
                  full_trajectory_executions=0, QR_solves=0, spectral_calculations=0,
                  note="Cached Yacht preparation; only two-step prefixes of existing exact-gradient cases")
    write(out / "limited_check.json", result)
    require(result["accepted"], "Prefix comparison failed; do not repeat the numerical check automatically")
    return result
