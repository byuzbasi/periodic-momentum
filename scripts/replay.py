"""Explicit, resumable replay of the existing design; never called by reporting."""
import json
import os
import time
from mpmath import mp
from common import (ROOT, configs, digest, kernels, output_path, read, require,
                    signature, write, frozen_check)


def recompute(destination, confirm=False, resume=False):
    require(confirm, "Numerical replay requires --confirm-frozen-replay")
    frozen_check()
    numeric, dense = kernels()
    scfg = configs()
    out = output_path(destination, create=not resume)
    sign = signature()
    if resume:
        require(read(out / "signature.json") == sign, "Resume signature mismatch")
        receipts = read(out / "receipts.json")
    else:
        (out / "synthetic").mkdir()
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
    planned = ["synthetic/" + u["id"] + ".json" for u in scfg["units"]]

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

    for unit in scfg["units"]:
        commit_unit("synthetic/" + unit["id"] + ".json", lambda unit=unit: numeric.evaluate(scfg, unit))
    require(set(receipts) == set(planned), "Incomplete replay")
    from records import validate
    checked = validate(out)
    report = dict(accepted=True, validation=checked, committed_units=len(receipts),
                  trajectory_executions=32, gradient_updates=16384, QR_solves=0,
                  scope="Exact frozen-design numerical replay, not an expanded experiment")
    if not (out / "replay_check.json").exists():
        write(out / "replay_check.json", report)
    return report
