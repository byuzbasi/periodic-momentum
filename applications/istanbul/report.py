"""Analyze only accepted saved records; no numerical trajectory calls."""
import argparse
import csv
import json
from pathlib import Path
from adapter import ROOT,mp,read,write,require
from run_pilot import verify


def csv_file(path,rows):
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def display(logvalue): return mp.nstr(mp.exp(logvalue),16)


def markdown(rows,columns):
    return '\n'.join(['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']+
        ['| '+' | '.join(str(r.get(c,'')) for c in columns)+' |' for r in rows])


def report(destination):
    cfg=read(ROOT/'config.json'); checked=verify(cfg)
    require(not destination.exists(),'Choose a new analysis directory')
    with mp.workdps(160):
        series=[]
        for path in sorted((ROOT/'accepted/records').glob('*.json')):
            row=read(path);unit=row['unit']
            if unit['digits']!=160:continue
            trace={p['k']:{k:mp.mpf(p[k]) for k in ['ln_relative_x','ln_relative_f','ln_relative_state']}
                   for p in row['trace']}
            series.append(dict(profile=unit['profile'],method=unit['method'],trace=trace))
        require(len(series)==8 and len({(s['profile'],s['method']) for s in series})==8,'Primary-series count')
        traces,endpoints,pairs,comparisons=[],[],[],[]
        for s in series:
            meta={k:s[k] for k in ['profile','method']};t=s['trace']
            require(sorted(t)==list(range(513)),'Incomplete saved series')
            for k,p in t.items():
                traces.append(meta|dict(k=k,log10_D=float(p['ln_relative_x']/mp.log(10)),
                    log10_F=float(p['ln_relative_f']/mp.log(10)),log10_state=float(p['ln_relative_state']/mp.log(10))))
            for k in range(2,513,2):
                pairs.append(meta|dict(even_k=k,
                    log10_max_D=float(max(t[k-1]['ln_relative_x'],t[k]['ln_relative_x'])/mp.log(10)),
                    log10_max_F=float(max(t[k-1]['ln_relative_f'],t[k]['ln_relative_f'])/mp.log(10))))
            endpoints.append(meta|dict(D511=display(t[511]['ln_relative_x']),D512=display(t[512]['ln_relative_x']),
                F511=display(t[511]['ln_relative_f']),F512=display(t[512]['ln_relative_f']),
                pair_max_D=display(max(t[k]['ln_relative_x'] for k in [511,512])),
                pair_max_F=display(max(t[k]['ln_relative_f'] for k in [511,512])),
                peak_D=display(max(p['ln_relative_x'] for p in t.values())),
                peak_F=display(max(p['ln_relative_f'] for p in t.values()))))
        tol=mp.mpf(cfg['tolerances']['comparison_log'])
        for profile in cfg['profiles']:
            sub={s['method']:s['trace'] for s in series if s['profile']==profile};chb=sub['CHB']
            for method in ['HB','GD','NAG-Q']:
                other=sub[method];row=dict(profile=profile,comparator=method)
                for label,metric in [('D','ln_relative_x'),('F','ln_relative_f')]:
                    row[f'CHB_to_other_{label}511']=display(chb[511][metric]-other[511][metric])
                    row[f'CHB_to_other_{label}512']=display(chb[512][metric]-other[512][metric])
                    row[f'CHB_to_other_pair_max_{label}']=display(
                        max(chb[k][metric] for k in [511,512])-max(other[k][metric] for k in [511,512]))
                    for phase,indices in [('odd',range(1,512,2)),('even',range(2,513,2)),('all',range(1,513))]:
                        differences=[chb[k][metric]-other[k][metric] for k in indices]
                        row[f'{label}_{phase}_CHB_wins']=sum(x < -tol for x in differences)
                        row[f'{label}_{phase}_CHB_losses']=sum(x > tol for x in differences)
                        row[f'{label}_{phase}_ties']=sum(abs(x)<=tol for x in differences)
                    differences=[max(chb[j][metric] for j in [k-1,k])-max(other[j][metric] for j in [k-1,k])
                                 for k in range(2,513,2)]
                    row[f'{label}_pair_CHB_wins']=sum(x < -tol for x in differences)
                    row[f'{label}_pair_CHB_losses']=sum(x > tol for x in differences)
                    row[f'{label}_pair_ties']=sum(abs(x)<=tol for x in differences)
                    nonwins=[k for k in range(1,513) if chb[k][metric]-other[k][metric]>=-tol]
                    row[f'{label}_CHB_winning_suffix_start']=max(nonwins)+1 if nonwins and max(nonwins)<512 else (1 if not nonwins else None)
                comparisons.append(row)
        destination.mkdir(parents=True)
        for name,rows in [('trajectories',traces),('paired_maxima',pairs),('endpoint_metrics',endpoints),('CHB_comparisons',comparisons)]:
            csv_file(destination/(name+'.csv'),rows)
        write(destination/'numerical_quality.json',checked)
        scope=dict(primary_precision=160,arithmetic_check_precision=120,primary_series=8,trace_rows=len(traces),
            pair_rows=len(pairs),endpoint_rows=8,comparison_rows=6,all_methods_new_on_same_data=True,
            profile_policy=cfg['profiles'],random_repetitions=0,financial_forecast=False,
            source_response='ISE USD BASED, C column',new_work_after_accepted_run='saved-output analysis only')
        write(destination/'comparison_scope.json',scope)
        sections=['# Istanbul: accepted four-method pilot','',
            'All four methods use the same 536 observations, seven raw inputs, USD-based ISE response, rest initialization and 512 gradient queries. Lower D and F are better. The 160-digit records are primary; 120 digits are an arithmetic check.',
            'F is normalized quadratic energy about the computed QR reference. It approximates the exact objective gap; the residual correction and reference-error bounds are documented in ../REFERENCE_ACCURACY.md. The accepted metric values are unchanged.',
            '', '## Adjacent terminal queries and peaks','',
            markdown(endpoints,['profile','method','D511','D512','F511','F512','pair_max_D','pair_max_F','peak_D','peak_F']),
            '', '## CHB / comparator ratios','',
            markdown(comparisons,['profile','comparator','CHB_to_other_D512','CHB_to_other_F512','CHB_to_other_pair_max_D','CHB_to_other_pair_max_F']),
            '', '## Counts across all 512 positive queries','',
            markdown(comparisons,['profile','comparator','D_all_CHB_wins','D_all_CHB_losses','D_all_ties','F_all_CHB_wins','F_all_CHB_losses','F_all_ties']),
            '', '## Counts across 256 odd/even pairs','',
            markdown(comparisons,['profile','comparator','D_pair_CHB_wins','D_pair_CHB_losses','F_pair_CHB_wins','F_pair_CHB_losses','D_CHB_winning_suffix_start','F_CHB_winning_suffix_start']),
            '', '## Scope','',
            '16 full trajectories plus four two-step prefixes: 8,200 total updates. Two QR references and two high-precision eigensolves were performed once. No retuning, additional dataset, longer horizon or automatic retry occurred.',
            'Both metrics and both phases are reported, including peaks and early losses. This one fixed quadratic is an application illustration. It establishes neither broad empirical superiority, financial prediction accuracy, wall-clock speedups nor the fourth-power law by itself.','']
        (destination/'RESULTS.md').write_text('\n'.join(sections))
        print(json.dumps(scope,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='analysis_v01')
    args=parser.parse_args();dest=Path(args.out)
    report(dest if dest.is_absolute() else ROOT/dest)
