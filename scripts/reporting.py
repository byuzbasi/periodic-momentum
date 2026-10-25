"""Reporting adapters from the original exporters; no trajectory or fit is run."""
import csv
import math
from mpmath import mp
from common import ROOT, configs, read, require, write
from records import load_records, validate


def csv_file(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export(out, directory=None):
    validation = validate(directory)
    scfg, ycfg = configs()
    srows, yrows, _ = load_records(directory)
    sout = out / "periodic_momentum_numerics_v06/analysis"
    yout = out / "periodic_momentum_yacht_v09/analysis"
    sout.mkdir(parents=True)
    yout.mkdir(parents=True)
    with mp.workdps(120):
        scales, cases, traces = [], [], []
        for row in srows[:5]:
            r, p = row["results"], row["results"]["parameters"]
            current = dict(id=row["unit"]["id"], R_label=row["unit"]["R"], R=float(p["R"]),
                epsilon_upper=float(p["epsilon_upper"]), epsilon_lower=float(p["epsilon_lower"]),
                coarse_upper=float(p["coarse_upper"]), ratio=float(r["ratio"]), limit=float(p["limit"]),
                relative_limit_error=float(r["relative_limit_error"]), local_exponent='',
                scalar_rate_half=float(r["rate_checks"][0]["rate"]),
                scalar_rate_at_upper=float(r["rate_checks"][1]["rate"]),
                scalar_rate_above=float(r["rate_checks"][2]["rate"]))
            if scales:
                current["local_exponent"] = math.log(scales[-1]["epsilon_upper"] / current["epsilon_upper"]) / math.log(2)
            scales.append(current)
        for row in srows[5:]:
            p = row["parameters"]
            cases.append(dict(id=row["unit"]["id"], R_label=row["unit"]["R"], R=float(p["R"]), profile=row["unit"]["profile"],
                epsilon=float(row["epsilon"]), epsilon_lower=float(p["epsilon_lower"]),
                within_proven_sufficient_region=float(row["epsilon"]) < float(p["epsilon_lower"]),
                epsilon_to_sufficient_bound=float(row["epsilon"]) / float(p["epsilon_lower"]),
                terminal_chb_to_hb_ratio=float(row["terminal_chb_to_hb_ratio"]),
                observed_suffix_win_start=row["observed_suffix_win_start"],
                chb_finite_window_factor=float(row["methods"]["CHB"]["late_window_factor"]),
                hb_finite_window_factor=float(row["methods"]["HB"]["late_window_factor"]),
                chb_mode_matrix_rate=row["methods"]["CHB"]["mode_matrix_rate"],
                hb_mode_matrix_rate=row["methods"]["HB"]["mode_matrix_rate"]))
            for method, result in row["methods"].items():
                for t in result["trace"]:
                    traces.append(dict(id=row["unit"]["id"], R_label=row["unit"]["R"],
                                       profile=row["unit"]["profile"], method=method, **t))
        csv_file(sout / "scale.csv", scales)
        csv_file(sout / "cases.csv", cases)
        csv_file(sout / "trajectories.csv", traces)
        primary = {(r["unit"]["profile"], r["unit"]["method"]): r
                   for r in yrows if r["unit"]["digits"] == 120}
        cases, traces, effects, phases = [], [], {}, {}
        for profile in ycfg["profiles"]:
            records = {m: primary[(profile, m)] for m in ycfg["methods"]}
            for method, r in records.items():
                for t in r["trace"]:
                    traces.append(dict(profile=profile, method=method, k=t["k"],
                        log10_relative_x=float(mp.mpf(t["ln_relative_x"]) / mp.log(10)),
                        log10_relative_f=float(mp.mpf(t["ln_relative_f"]) / mp.log(10))))
            a, b = [records[m]["trace"] for m in ["CHB", "HB"]]
            xdiff = [mp.mpf(x["ln_relative_x"]) - mp.mpf(y["ln_relative_x"]) for x, y in zip(a, b)]
            fdiff = [mp.mpf(x["ln_relative_f"]) - mp.mpf(y["ln_relative_f"]) for x, y in zip(a, b)]
            last = max(i for i, v in enumerate(xdiff) if v >= -mp.mpf(ycfg["tolerances"]["tie_log"]))
            row = dict(profile=profile, epsilon=float(records["CHB"]["epsilon"]),
                terminal_x_chb_to_hb_ratio=float(mp.exp(xdiff[-1])),
                terminal_f_chb_to_hb_ratio=float(mp.exp(fdiff[-1])),
                observed_suffix_x_win_start=last+1 if last < 512 else '')
            effects[profile] = {}
            for method, r in records.items():
                row[method.lower()+"_terminal_relative_x"] = float(mp.exp(mp.mpf(r["trace"][-1]["ln_relative_x"])))
                row[method.lower()+"_terminal_relative_f"] = float(mp.exp(mp.mpf(r["trace"][-1]["ln_relative_f"])))
                exact = primary[("exact", method)]
                effects[profile][method] = {k: mp.nstr(mp.exp(mp.mpf(r["trace"][-1][k]) - mp.mpf(exact["trace"][-1][k])) - 1, 100)
                                           for k in ["ln_relative_x", "ln_relative_f"]}
            cases.append(row)
            pairs = [("distance_CHB_HB", "ln_relative_x"), ("objective_CHB_HB", "ln_relative_f"), ("state_CHB_HB", "ln_relative_state")]
            checkpoints = [dict(k=k, **{key: float(mp.exp(mp.mpf(a[k][src]) - mp.mpf(b[k][src]))) for key, src in pairs})
                           for k in [1,2,16,64,128,256,509,510,511,512]]
            wins = {label: sum(mp.mpf(a[k][src]) < mp.mpf(b[k][src]) for k in range(1,513))
                    for label, src in [("distance_CHB_wins", "ln_relative_x"), ("objective_CHB_wins", "ln_relative_f")]}
            for src in ["ln_relative_x", "ln_relative_f"]:
                require(all(mp.mpf(a[k][src]) < mp.mpf(b[k][src]) for k in range(2,513,2)), "Even phase changed")
                require(all(mp.mpf(a[k][src]) > mp.mpf(b[k][src]) for k in range(1,512,2)), "Odd phase changed")
            phases[profile] = dict(saved_checkpoint_comparisons=checkpoints, **wins,
                description='Post-run phase diagnostics from existing traces; no new horizon or trajectory.')
        csv_file(yout / "cases.csv", cases)
        csv_file(yout / "trajectories.csv", traces)
        write(yout / "noise_effects.json", effects)
        write(yout / "phase_diagnostics.json", phases)
        adjacent = []
        a, b = [primary[("exact", m)]["trace"] for m in ["CHB", "HB"]]
        for k in [511,512]:
            adjacent.append(dict(k=k, D_CHB=float(mp.exp(mp.mpf(a[k]["ln_relative_x"]))),
                D_HB=float(mp.exp(mp.mpf(b[k]["ln_relative_x"]))),
                D_ratio=float(mp.exp(mp.mpf(a[k]["ln_relative_x"])-mp.mpf(b[k]["ln_relative_x"]))),
                F_CHB=float(mp.exp(mp.mpf(a[k]["ln_relative_f"]))), F_HB=float(mp.exp(mp.mpf(b[k]["ln_relative_f"])))))
        csv_file(yout / "adjacent_queries.csv", adjacent)
    matches = {}
    for group, folder in [("synthetic", sout), ("yacht", yout)]:
        for expected in sorted((ROOT / "reference" / group).iterdir()):
            actual = folder / expected.name
            equal = read(actual) == read(expected) if expected.suffix == ".json" else actual.read_bytes() == expected.read_bytes()
            matches[group + "/" + expected.name] = equal
            require(equal, f"Reporting mismatch: {group}/{expected.name}")
    report = dict(accepted=True, validation=validation, reference_matches=matches,
                  adjacent_queries=adjacent, new_gradient_updates=0)
    write(out / "reporting_check.json", report)
    return report
