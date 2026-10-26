"""Export only accepted saved synthetic results."""
import csv, math
from mpmath import mp
from common import ROOT, configs, read, require, write
from records import load_records, validate
def csv_file(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def export(out,directory=None):
    validation=validate(directory);srows=load_records(directory)
    sout=out/'periodic_momentum_numerics_v06/analysis';sout.mkdir(parents=True)
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
    for expected in (ROOT/'reference/synthetic').glob('*.csv'):
        require((sout/expected.name).read_bytes()==expected.read_bytes(),'Synthetic reporting mismatch')
    result=dict(accepted=True,validation=validation,reference_matches=3,new_gradient_updates=0)
    write(out/'reporting_check.json',result)
    return result
