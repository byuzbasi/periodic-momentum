"""Compare saved GD/NAG records with accepted CHB/HB data; no optimizer calls."""
import argparse
import csv
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from adapter import ROOT, PUBLIC, mp, numeric, read, write, require
from run_pilot import verify


def parsed_trace(trace):
    return {p['k']: {k:mp.mpf(str(v)) for k,v in p.items() if k!='k'
                    and k in ['ln_relative_x','ln_relative_f','ln_relative_state']} for p in trace}


def all_series(cfg):
    output = []
    for path in sorted((ROOT/'accepted/records').glob('*.json')):
        row = read(path)
        if row['unit']['digits'] == 120:
            u = row['unit']
            output.append(dict(case=u['case'], profile=u['profile'], method=u['method'],
                               trace=parsed_trace(row['trace']), source='new_120_digit_record'))
    for path in sorted((PUBLIC/'data/accepted/synthetic').glob('T*.json')):
        row = read(path); u = row['unit']
        if u['profile'] not in cfg['profiles']:
            continue
        case = 'synthetic_R_half' if u['R']=='1/2' else 'synthetic_R_eighth'
        H = mp.diag([mp.mpf(x) for x in row['parameters']['lambdas']])
        initial = mp.matrix([1]*4); initial_gap = (initial.T*H*initial)[0]/2
        for method, result in row['methods'].items():
            trace = parsed_trace(result['trace']); trace[0]['ln_relative_f'] = mp.mpf(0)
            for n, state in result['checkpoints'].items():
                state = [mp.mpf(x) for x in state]; n = int(n)
                for k, vector in [(n,state[:4]),(n-1,state[4:])]:
                    v = mp.matrix(vector)
                    trace[k]['ln_relative_f'] = mp.log((v.T*H*v)[0]/(2*initial_gap))
            output.append(dict(case=case, profile=u['profile'], method=method, trace=trace,
                               source='accepted_synthetic_float64_norm_logs_and_saved_states'))
    require(len(output)==16, 'Expected 2 cases x 2 profiles x 4 methods')
    require(len({(s['case'],s['profile'],s['method']) for s in output})==16, 'Duplicate series')
    return output


def csv_file(path, rows):
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def value(logvalue):
    return mp.nstr(mp.exp(logvalue),16) if logvalue is not None else ''


def table(rows, columns):
    return '\n'.join(['| '+' | '.join(columns)+' |',
        '| '+' | '.join(['---']*len(columns))+' |']+
        ['| '+' | '.join(str(row.get(c,'')) for c in columns)+' |' for row in rows])


def report(destination):
    cfg=read(ROOT/'config.json'); checked=verify(cfg)
    require(not destination.exists(),'Choose a new analysis directory; saved reports are preserved')
    with mp.workdps(120):
        series=all_series(cfg)
        traces=[]; endpoints=[]; pairs=[]; comparisons=[]
        for s in series:
            meta={k:s[k] for k in ['case','profile','method']}; trace=s['trace']
            require(sorted(trace)==list(range(513)), 'Incomplete trace')
            for k,p in trace.items():
                traces.append(meta | dict(k=k, log10_D=float(p['ln_relative_x']/mp.log(10)),
                    log10_F=float(p['ln_relative_f']/mp.log(10)) if 'ln_relative_f' in p else '',
                    log10_state=float(p['ln_relative_state']/mp.log(10)),source=s['source']))
            for k in range(2,513,2):
                odd,even=trace[k-1],trace[k]
                f=max(odd['ln_relative_f'],even['ln_relative_f']) if all('ln_relative_f' in p for p in [odd,even]) else None
                pairs.append(meta | dict(even_k=k,log10_max_D=float(max(odd['ln_relative_x'],even['ln_relative_x'])/mp.log(10)),
                    log10_max_F=float(f/mp.log(10)) if f is not None else ''))
            fcomplete=all('ln_relative_f' in p for p in trace.values())
            endpoints.append(meta | dict(D511=value(trace[511]['ln_relative_x']), D512=value(trace[512]['ln_relative_x']),
                F511=value(trace[511]['ln_relative_f']),F512=value(trace[512]['ln_relative_f']),
                pair_max_D=value(max(trace[k]['ln_relative_x'] for k in [511,512])),
                pair_max_F=value(max(trace[k]['ln_relative_f'] for k in [511,512])),
                peak_D=value(max(p['ln_relative_x'] for p in trace.values())),
                peak_F=value(max(p['ln_relative_f'] for p in trace.values())) if fcomplete else '',
                full_objective_trace=fcomplete))
        for case in cfg['cases']:
            for profile in cfg['profiles']:
                sub={s['method']:s for s in series if s['case']==case['id'] and s['profile']==profile}
                chb=sub['CHB']['trace']
                tol=mp.mpf(cfg['tolerances']['comparison_log_synthetic'])
                for method in ['HB','GD','NAG-Q']:
                    other=sub[method]['trace']
                    row=dict(case=case['id'],profile=profile,comparator=method,
                        CHB_to_other_D511=value(chb[511]['ln_relative_x']-other[511]['ln_relative_x']),
                        CHB_to_other_D512=value(chb[512]['ln_relative_x']-other[512]['ln_relative_x']),
                        CHB_to_other_pair_max_D=value(max(chb[k]['ln_relative_x'] for k in [511,512])-max(other[k]['ln_relative_x'] for k in [511,512])),
                        CHB_to_other_F512=value(chb[512]['ln_relative_f']-other[512]['ln_relative_f']),
                        CHB_to_other_pair_max_F=value(max(chb[k]['ln_relative_f'] for k in [511,512])-max(other[k]['ln_relative_f'] for k in [511,512])))
                    for label,indices in [('odd',range(1,512,2)),('even',range(2,513,2))]:
                        differences=[chb[k]['ln_relative_x']-other[k]['ln_relative_x'] for k in indices]
                        row[label+'_CHB_wins']=sum(x < -tol for x in differences)
                        row[label+'_CHB_losses']=sum(x > tol for x in differences)
                        row[label+'_ties']=sum(abs(x)<=tol for x in differences)
                    comparisons.append(row)
        destination.mkdir(parents=True)
        csv_file(destination/'trajectories.csv',traces)
        csv_file(destination/'paired_maxima.csv',pairs)
        csv_file(destination/'endpoint_metrics.csv',endpoints)
        csv_file(destination/'CHB_comparisons.csv',comparisons)
        write(destination/'numerical_quality.json',checked)
        write(destination/'comparison_scope.json',dict(series=16,trace_rows=len(traces),pair_rows=len(pairs),
            methods=['CHB','HB','GD','NAG-Q'],published_comparator_trajectories=16,reused_baseline_series=8,
            synthetic_baseline_objective_available_at=[0,1,2,15,16,63,64,511,512],
            baseline_recomputed=False,unreported_synthetic_objective_points='missing; no interpolation'))
        text=['# GD/NAG bounded comparison: saved results','',
              'NAG-Q is the constant-parameter Nesterov tuning for strongly convex quadratics. All trajectories use 512 gradient queries. Lower D and F are better.',
              '', '## Distance and objective at the final adjacent pair','',
              table(endpoints,['case','profile','method','D511','D512','F511','F512','pair_max_D','pair_max_F']),
              '', '## CHB / comparator ratios','',
              table(comparisons,['case','profile','comparator','CHB_to_other_D511','CHB_to_other_D512','CHB_to_other_pair_max_D','CHB_to_other_F512']),
              '', '## Odd/even distance counts (256 per phase)','',
              table(comparisons,['case','profile','comparator','odd_CHB_wins','odd_CHB_losses','odd_ties','even_CHB_wins','even_CHB_losses','even_ties']),
              '', '## Validation and scope','',
              'The published synthetic subset contains 16 accepted full trajectories and two prefixes: 8,192 + 4 = 8,196 historical updates. No new CHB/HB trajectory, QR solve or eigensolve.',
              'All new trajectories passed independent affine-state and relative-error checks, and 80/120-digit log comparisons. See numerical_quality.json.',
              'Synthetic CHB/HB objective values exist only at saved state points. Their full objective peaks are intentionally missing.',
              'Synthetic baseline norm logs were stored as float64: comparisons use a 1e-12 log tie tolerance. All reported distances refer to the same reference solution and initialization.',
              'No claim about arbitrary gradient-error paths, parameter-optimal noisy methods, out-of-sample prediction or general empirical superiority follows from this pilot.','']
        (destination/'results.md').write_text('\n'.join(text))
        print(json.dumps(dict(accepted=True,series=len(series),trace_rows=len(traces),endpoint_rows=len(endpoints),comparison_rows=len(comparisons)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',default='outputs/gd_nag_reports')
    args=parser.parse_args()
    target=Path(args.out)
    report(target if target.is_absolute() else PUBLIC/target)
