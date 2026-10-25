"""Explicit, resumable replay of the existing design; never called by reporting."""
import json
import os
import time
from mpmath import mp
from common import (ROOT, configs, digest, kernels, output_path, read, require,
                    signature, write, yacht_units, frozen_check)


def recompute(destination, confirm=False, resume=False):
    require(confirm, "Numerical replay requires --confirm-frozen-replay")
    frozen_check()
    numeric, dense = kernels()
    scfg, ycfg = configs()
    out = output_path(destination, create=not resume)
    sign = signature()
    if resume:
        require(read(out / "signature.json") == sign, "Resume signature mismatch")
        receipts = read(out / "receipts.json")
    else:
        (out / "synthetic").mkdir()
        (out / "yacht").mkdir()
        write(out / "signature.json", sign)
        write(out / "receipts.json", {})
        receipts = {}
    for rel, expected in receipts.items():
        path = out / rel
        require(path.stat().st_size == expected["bytes"] and digest(path) == expected["sha256"], "Invalid replay checkpoint")
    attempted = {row["unit"] for row in [json.loads(line) for line in
        (out / "attempts.jsonl").read_text().splitlines()]} if (out / "attempts.jsonl").exists() else set()
    # A started but uncommitted unit is retained for inspection, never retried silently.
    require(attempted <= set(receipts), "An interrupted unit needs inspection; automatic retry is disabled")
    started = time.perf_counter()
    planned = [f"yacht/setup_{d}.json" for d in [80,120]]
    planned += ["synthetic/" + u["id"] + ".json" for u in scfg["units"]]
    planned += ["yacht/" + u["id"] + ".json" for u in yacht_units(ycfg)]

    def commit_unit(rel, compute):
        if rel in receipts:
            return read(out / rel)
        require(rel not in attempted, "Unit already attempted")
        with (out / "attempts.jsonl").open("a") as stream:
            json.dump({"unit": rel}, stream)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        attempted.add(rel)
        value = compute()
        require(value.get("accepted", True), f"Acceptance failed: {rel}")
        write(out / rel, value)
        receipts[rel] = {"bytes": (out / rel).stat().st_size, "sha256": digest(out / rel)}
        write(out / "receipts.json", receipts, exclusive=False)
        progress = dict(total_units=len(planned), durably_committed=len(receipts), last_unit=rel,
                        percentage=100*len(receipts)/len(planned), elapsed_seconds=time.perf_counter()-started)
        write(out / "progress.json", progress, exclusive=False)
        print(f"{len(receipts)}/{len(planned)} committed: {rel}", flush=True)
        return value

    for digits in ycfg["precisions"]:
        def prepare(digits=digits):
            with mp.workdps(digits):
                return dense.pack_setup(dense.preparation(ycfg, read(ROOT / "configs/yacht_bounds.json"),
                                                         ROOT / "data/raw/yacht_hydrodynamics.data"))
        commit_unit(f"yacht/setup_{digits}.json", prepare)
    for unit in scfg["units"]:
        commit_unit("synthetic/" + unit["id"] + ".json", lambda unit=unit: numeric.evaluate(scfg, unit))
    for unit in yacht_units(ycfg):
        def calculate(unit=unit):
            with mp.workdps(unit["digits"]):
                setup = dense.unpack_setup(read(out / "yacht" / f"setup_{unit['digits']}.json"))
                counter = {"updates": 0}
                result = dense.trajectory(ycfg, setup, unit, counter)
                require(counter["updates"] == 512, "Unexpected replay update count")
                return result
        commit_unit("yacht/" + unit["id"] + ".json", calculate)
    require(set(receipts) == set(planned), "Incomplete replay")
    from records import validate
    checked = validate(out)
    report = dict(accepted=True, validation=checked, committed_units=len(receipts),
                  trajectory_executions=48, gradient_updates=24576, QR_solves=2,
                  scope="Exact frozen-design numerical replay, not an expanded experiment")
    if not (out / "replay_check.json").exists():
        write(out / "replay_check.json", report)
    return report
