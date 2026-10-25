"""One bounded, resumable local pilot; persistent reservations prevent retries."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import signal
import sys
import time
sys.dont_write_bytecode = True
from adapter import ROOT, PUBLIC, mp, numeric, units, trajectory, read, write, digest, require, frozen_check


def signature():
    paths = [ROOT/n for n in ['config.json', 'protocol.md', 'adapter.py', 'run_pilot.py', 'report.py', 'plot.R']]
    paths += [p for folder in ['code', 'scripts', 'configs', 'data']
              for p in (PUBLIC/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    paths += [PUBLIC/'frozen_inputs.json', PUBLIC/'delivery_manifest.json']
    return dict(environment={'python': platform.python_version(), 'mpmath': importlib.metadata.version('mpmath')},
                files={str(p.relative_to(PUBLIC)): digest(p) for p in sorted(paths)})


def design(cfg):
    planned = units(cfg)
    require(len(planned) == cfg['limits']['full_trajectories'] == 24, 'Unexpected design')
    require(cfg['horizon'] == 512 and len(planned)*512 == 12288, 'Unexpected horizon')
    require(cfg['limits']['total_gradient_updates'] == 12296, 'Unexpected total cap')
    require(cfg['precisions'] == [80,120], 'Unexpected precision')
    return planned


class Counter:
    def __init__(self, reserved):
        self.reserved, self.used = reserved, 0

    def charge(self):
        require(self.used < self.reserved, 'Reserved gradient limit reached')
        self.used += 1


def quality(cfg, run):
    planned = units(cfg)
    rows = [read(run/'records'/f"{u['id']}.json") for u in planned]
    require(all(r['accepted'] and len(r['trace']) == 513 and r['gradient_updates'] == 512 for r in rows), 'Bad trajectory')
    require(all(r['unit'] == u and [p['k'] for p in r['trace']] == list(range(513))
                for r,u in zip(rows,planned)), 'Work-unit identity or trace order mismatch')
    maximum = mp.mpf(0)
    with mp.workdps(120):
        require(all(mp.isfinite(mp.mpf(p[key])) for row in rows for p in row['trace']
                    for key in ['ln_relative_x','ln_relative_f','ln_relative_state']), 'Nonfinite trace')
        for low, high in zip(rows[::2], rows[1::2]):
            require(low['unit']['digits'] == 80 and high['unit']['digits'] == 120, 'Bad precision pairing')
            for metric in ['ln_relative_x','ln_relative_f','ln_relative_state']:
                maximum = max(maximum, max(abs(mp.mpf(a[metric])-mp.mpf(b[metric]))
                    for a,b in zip(low['trace'],high['trace'])))
        require(maximum <= mp.mpf(cfg['tolerances']['precision_log']), 'Precision gate failed')
        return numeric.serialize(dict(accepted=True, full_trajectories=len(rows), main_updates=12288,
            max_80_120_log_difference=maximum,
            max_matrix_relative_error=max(mp.mpf(r['max_matrix_relative_error']) for r in rows),
            max_oracle_violation=max(mp.mpf(r['max_oracle_violation']) for r in rows),
            new_QR_solves=0, new_eigensolves=0, new_CHB_HB_trajectories=0))


def verify(cfg, directory=None):
    run = ROOT/'accepted' if directory is None else Path(directory)
    if directory is None:
        published = read(ROOT/'accepted/manifest.json')
        for relative, expected in published['files'].items():
            path = ROOT/relative
            require(path.stat().st_size == expected['bytes'] and digest(path) == expected['sha256'], 'Published comparison file changed: '+relative)
        for relative, expected in published['inputs'].items():
            require(digest(PUBLIC/relative) == expected, 'Frozen baseline input changed: '+relative)
        frozen_check()
    else:
        require(read(run/'signature.json') == signature(), 'Scientific/software signature changed')
    ledger = read(run/'ledger.json')
    require(ledger['reserved_updates'] <= 12296, 'Budget overrun')
    for unit, record in ledger['attempts'].items():
        require(record['status'] == 'complete', f'Incomplete attempted unit: {unit}')
        path = run/record['path']
        require(path.stat().st_size == record['bytes'] and digest(path) == record['sha256'], 'Receipt mismatch')
    expected = {u['id'] for u in units(cfg)} | {f'SM{i+1:02d}' for i in range(4)}
    require(set(ledger['attempts']) == expected, 'Missing work units')
    require(sum(r['used_updates'] for r in ledger['attempts'].values()) == 12296, 'Incorrect executed count')
    checked = quality(cfg, run)
    checked.update(prefix_checks=4, prefix_updates=8, total_updates=12296,
                   numerical_wall_seconds=ledger['elapsed_seconds'], originals_unchanged=True)
    return checked


def run(cfg, resume):
    planned = design(cfg)
    out = ROOT/'outputs/run_001'
    started = time.perf_counter()
    if resume:
        require(read(out/'signature.json') == signature(), 'Resume signature mismatch')
        ledger = read(out/'ledger.json')
        require(all(v['status'] == 'complete' for v in ledger['attempts'].values()), 'Interrupted unit; no automatic retry')
        for v in ledger['attempts'].values():
            require(digest(out/v['path']) == v['sha256'], 'Invalid checkpoint')
    else:
        require(not out.exists(), 'Output exists; use verified resume, never overwrite')
        out.mkdir(parents=True); (out/'records').mkdir(); (out/'prefix').mkdir()
        write(out/'signature.json', signature())
        ledger = dict(reserved_updates=0, elapsed_seconds=0, attempts={})
        write(out/'ledger.json', ledger)
    previous_elapsed = ledger['elapsed_seconds']
    remaining = cfg['limits']['cumulative_numerical_seconds']-previous_elapsed
    require(remaining > 0, 'Cumulative numerical time budget exhausted')
    def timeout(*_):
        raise TimeoutError('120-second cumulative pilot boundary reached')
    signal.signal(signal.SIGALRM, timeout); signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        frozen_check()
        from records import validate
        if not (out/'baseline_validation.json').exists():
            write(out/'baseline_validation.json', validate())
        prefixes = [dict(id=f'SM{i+1:02d}',case=c,profile='exact',method=m,digits=120)
                    for i,(c,m) in enumerate((c,m) for c in ['synthetic_R_half','yacht'] for m in cfg['methods'])]
        for unit in prefixes+planned:
            key = unit['id']; prefix = key.startswith('SM'); steps = 2 if prefix else 512
            if key in ledger['attempts']:
                continue
            require(ledger['reserved_updates']+steps <= 12296, 'No remaining query reservation')
            path = Path('prefix' if prefix else 'records')/(key+'.json')
            ledger['reserved_updates'] += steps
            ledger['attempts'][key] = dict(status='started', reserved_updates=steps, path=str(path))
            write(out/'ledger.json', ledger, exclusive=False)
            counter = Counter(steps); before = time.perf_counter()
            try:
                result = trajectory(cfg, unit, counter, steps=steps)
                require(counter.used == steps and result['accepted'], 'Numerical acceptance failed')
                result['wall_seconds'] = time.perf_counter()-before
                write(out/path, result)
            except BaseException as error:
                ledger['attempts'][key].update(status='failed', used_updates=counter.used, error=str(error))
                raise
            ledger['attempts'][key].update(status='complete', used_updates=counter.used,
                bytes=(out/path).stat().st_size, sha256=digest(out/path))
            ledger['elapsed_seconds'] = previous_elapsed+time.perf_counter()-started
            write(out/'ledger.json', ledger, exclusive=False)
            completed = sum(k.startswith('U') and v['status']=='complete' for k,v in ledger['attempts'].items())
            elapsed = time.perf_counter()-started
            progress = dict(total=24, completed=completed, pending=24-completed,
                failed=0, committed=completed, percent=100*completed/24,
                reserved_updates=ledger['reserved_updates'], executed_updates=sum(v.get('used_updates',0) for v in ledger['attempts'].values()),
                elapsed_seconds=ledger['elapsed_seconds'], last_unit=key,
                eta_seconds=(24-completed)*elapsed/completed if completed else None)
            write(out/'progress.json', progress, exclusive=False)
            with (out/'progress.tsv').open('a') as stream:
                if stream.tell()==0: stream.write('completed\ttotal\tactual_updates\telapsed_seconds\tlast_unit\n')
                stream.write(f"{completed}\t24\t{progress['executed_updates']}\t{ledger['elapsed_seconds']:.6f}\t{key}\n")
            print(f"{key}: accepted; complete={completed}/24, updates={progress['executed_updates']}, wall={result['wall_seconds']:.3f}s",flush=True)
        ledger['elapsed_seconds'] = previous_elapsed+time.perf_counter()-started
        write(out/'ledger.json', ledger, exclusive=False)
        checked = verify(cfg, out)
        if not (out/'quality.json').exists(): write(out/'quality.json', checked)
        print(json.dumps(checked, indent=2))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        ledger['elapsed_seconds'] = previous_elapsed+time.perf_counter()-started
        write(out/'ledger.json', ledger, exclusive=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['design','run','verify'])
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(); config = read(ROOT/'config.json')
    if args.command == 'design': print(json.dumps(dict(units=design(config), limits=config['limits']),indent=2))
    elif args.command == 'verify': print(json.dumps(verify(config),indent=2))
    else: run(config, args.resume)
