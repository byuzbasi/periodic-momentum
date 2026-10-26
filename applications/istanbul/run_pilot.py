"""One bounded pilot with persistent query/preparation reservations."""
import argparse
import fcntl
import importlib.metadata
import platform
import resource
import signal
import time
from adapter import ROOT, PUBLIC, mp, numeric, units, preparation, trajectory
from adapter import read, write, digest, require


def signature():
    paths=[ROOT/n for n in ['config.json','protocol.md','adapter.py','run_pilot.py','input_tokens.json']]
    paths += [PUBLIC/'scripts/common.py',PUBLIC/'code/periodic_momentum_numerics_v06/code/numerics.py',PUBLIC/'data/raw/data_akbilgic.xlsx']
    return dict(environment={'python':platform.python_version(),'mpmath':importlib.metadata.version('mpmath')},
        files={p.relative_to(PUBLIC).as_posix():digest(p) for p in paths})


def design(cfg):
    planned=units(cfg)
    require(len(planned)==cfg['limits']['full_trajectories']==16,'Unexpected design')
    require(cfg['horizon']==512 and len(planned)*512==8192,'Unexpected main budget')
    require(cfg['limits']['total_gradient_updates']==8200,'Unexpected total budget')
    require(cfg['precisions']==[120,160] and cfg['prefix_precision']==160,'Approved precision amendment')
    require(cfg['methods']==['CHB','HB','GD','NAG-Q'] and cfg['profiles']==['exact','certified_mixing'],'Unapproved methods/profiles')
    require(cfg['n']==536 and cfg['p']==7 and not cfg['retuning'],'Input model or tuning changed')
    return planned


class Counter:
    def __init__(self,reserved):
        self.reserved,self.used=reserved,0

    def charge(self):
        require(self.used<self.reserved,'Reserved gradient limit reached')
        self.used+=1


def quality(cfg,out):
    planned=units(cfg); rows=[read(out/'records'/f"{u['id']}.json") for u in planned]
    require(all(r['unit']==u and r['accepted'] and r['gradient_updates']==512 and
        len(r['trace'])==513 and [p['k'] for p in r['trace']]==list(range(513))
        for r,u in zip(rows,planned)),'Missing or invalid full trajectory')
    require(all(set(r['checkpoints'])==set(r['matrix_relative_errors'])==
        {str(k) for k in cfg['checkpoints']} for r in rows),'Incomplete independent state checks')
    with mp.workdps(160):
        maximum=mp.mpf(0)
        for low,high in zip(rows[::2],rows[1::2]):
            require(low['unit']['digits']==120 and high['unit']['digits']==160,'Bad precision pair')
            require(all(low['unit'][k]==high['unit'][k] for k in ['profile','method']),'Precision identity mismatch')
            for metric in ['ln_relative_x','ln_relative_f','ln_relative_state']:
                require(all(mp.isfinite(mp.mpf(p[metric])) for r in [low,high] for p in r['trace']),'Nonfinite log')
                maximum=max(maximum,max(abs(mp.mpf(a[metric])-mp.mpf(b[metric]))
                    for a,b in zip(low['trace'],high['trace'])))
        require(maximum<=mp.mpf(cfg['tolerances']['precision_log']),'120/160 precision gate failed')
        max_matrix=max(mp.mpf(r['max_matrix_relative_error']) for r in rows)
        max_oracle=max(mp.mpf(r['max_oracle_violation']) for r in rows)
        require(max_matrix<=mp.mpf(cfg['tolerances']['matrix_relative']),'Affine-state gate failed')
        require(max_oracle<=mp.mpf(cfg['tolerances']['oracle_relative']),'Relative-error gate failed')
        residuals=[]
        for digits in cfg['precisions']:
            s=read(out/f'setup_{digits}.json')
            require(s['working_digits']==s['stored_digits']==digits,'Truncated setup precision')
            require(len(s['H'])==len(s['star'])==7 and all(len(x)==7 for x in s['H']),'Setup dimensions')
            residuals.append(mp.mpf(s['normalized_reference_gradient']))
            p=s['p']; spectrum=list(map(mp.mpf,s['spectrum']))
            require(len(spectrum)==7 and all(mp.mpf(p['mu'])<=v<=mp.mpf(p['L']) and
                (v<=mp.mpf(p['a']) or v>=mp.mpf(p['b'])) for v in spectrum),'Saved spectral enclosure')
        require(max(residuals)<=mp.mpf(cfg['tolerances']['reference']),'Reference gate failed')
        return dict(accepted=True,full_trajectories=16,main_updates=8192,
            max_120_160_log_difference=mp.nstr(maximum,160),
            max_matrix_relative_error=mp.nstr(max_matrix,160),max_oracle_violation=mp.nstr(max_oracle,160),
            max_reference_gradient=mp.nstr(max(residuals),160))


def verify(cfg,directory=None):
    out=ROOT/'accepted' if directory is None else directory
    ledger=read(out/'ledger.json')
    if directory is not None:
        require(read(out/'signature.json')==signature(),'Scientific/software signature changed')
    else:
        for name,info in read(ROOT/'accepted/manifest.json')['files'].items():
            path=PUBLIC/name
            require(path.stat().st_size==info['bytes'] and digest(path)==info['sha256'],'Published input/record changed: '+name)
    prefixes={f'SM{i+1:02d}' for i in range(4)}
    expected={u['id'] for u in units(cfg)}|prefixes
    require(set(ledger['attempts'])==expected,'Missing work units')
    require(set(ledger['preparations'])=={'120','160'},'Missing preparation')
    for record in list(ledger['attempts'].values())+list(ledger['preparations'].values()):
        require(record['status']=='complete','Incomplete attempted unit; no automatic retry')
        path=out/record['path']
        require(path.stat().st_size==record['bytes'] and digest(path)==record['sha256'],'Receipt mismatch')
    require(ledger['reserved_updates']==8200,'Reservation count mismatch')
    require(sum(x['used_updates'] for x in ledger['attempts'].values())==8200,'Executed count mismatch')
    require(len(ledger['QR'])==len(ledger['eigen'])==2,'Preparation call budget mismatch')
    for i,method in enumerate(cfg['methods']):
        r=read(out/'prefix'/f'SM{i+1:02d}.json')
        require(r['accepted'] and r['unit']['method']==method and r['unit']['digits']==160 and
            r['steps']==r['gradient_updates']==2 and len(r['trace'])==3,'Prefix mismatch')
    checked=quality(cfg,out)
    checked.update(prefix_checks=4,prefix_updates=8,total_updates=8200,QR_reference_solves=2,
        high_precision_eigensolves=2,numerical_wall_seconds=ledger['elapsed_seconds'],
        peak_rss_bytes=ledger['peak_rss_bytes'])
    return checked


def run(cfg,resume):
    planned=design(cfg); out=ROOT/'outputs/run_001'
    import adapter
    adapter.RUN_DIRECTORY=out
    if not resume:
        require(not out.exists(),'Output exists; preserve it and use verified resume')
        out.mkdir(parents=True); (out/'records').mkdir(); (out/'prefix').mkdir()
        write(out/'signature.json',signature())
        ledger=dict(reserved_updates=0,elapsed_seconds=0,attempts={},preparations={},QR=[],eigen=[])
        write(out/'ledger.json',ledger)
    else:
        require(read(out/'signature.json')==signature(),'Resume signature mismatch')
        ledger=read(out/'ledger.json')
        for rec in list(ledger['attempts'].values())+list(ledger['preparations'].values()):
            require(rec['status']=='complete' and digest(out/rec['path'])==rec['sha256'],'Interrupted or invalid unit; inspect before retry')
    with (out/'run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        started=time.perf_counter(); previous_elapsed=ledger['elapsed_seconds']
        remaining=cfg['limits']['cumulative_numerical_seconds']-previous_elapsed
        require(remaining>0,'Cumulative numerical budget exhausted')
        def checkpoint():
            ledger['elapsed_seconds']=previous_elapsed+time.perf_counter()-started
            ledger['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            write(out/'ledger.json',ledger,exclusive=False)
        def timeout(*_): raise TimeoutError('120-second cumulative pilot boundary reached')
        signal.signal(signal.SIGALRM,timeout); signal.setitimer(signal.ITIMER_REAL,remaining)
        try:
            require(digest(PUBLIC/'data/raw/data_akbilgic.xlsx')==cfg['source_sha256'],'Raw source changed')
            def charge(kind,digits):
                require(kind in ['QR','eigen'] and len(ledger[kind])<2 and digits not in ledger[kind],'Preparation call reservation exhausted')
                ledger[kind].append(digits); checkpoint()
            for digits in cfg['precisions']:
                key=str(digits)
                if key in ledger['preparations']: continue
                path=f'setup_{digits}.json'
                ledger['preparations'][key]=dict(status='started',path=path); checkpoint()
                print(f'Preparing raw decimal data, QR and spectrum at {digits} digits',flush=True)
                try: result=preparation(cfg,digits,charge); write(out/path,result)
                except BaseException as error:
                    ledger['preparations'][key].update(status='failed',error=repr(error)); raise
                ledger['preparations'][key].update(status='complete',bytes=(out/path).stat().st_size,sha256=digest(out/path))
                checkpoint()
            prefixes=[dict(id=f'SM{i+1:02d}',profile='exact',method=m,digits=160) for i,m in enumerate(cfg['methods'])]
            for unit in prefixes+planned:
                key=unit['id']; prefix=key.startswith('SM'); steps=2 if prefix else 512
                if key in ledger['attempts']: continue
                require(ledger['reserved_updates']+steps<=8200,'Query reservation exhausted')
                path=f"{'prefix' if prefix else 'records'}/{key}.json"
                ledger['reserved_updates']+=steps
                ledger['attempts'][key]=dict(status='started',path=path,reserved_updates=steps); checkpoint()
                counter=Counter(steps); before=time.perf_counter()
                try:
                    result=trajectory(cfg,unit,counter,steps=steps)
                    require(counter.used==steps,'Query count mismatch')
                    result['wall_seconds']=time.perf_counter()-before
                    write(out/path,result)
                    require(result['accepted'],'Numerical acceptance failed; no automatic continuation')
                except BaseException as error:
                    ledger['attempts'][key].update(status='failed',used_updates=counter.used,error=repr(error)); raise
                ledger['attempts'][key].update(status='complete',used_updates=counter.used,
                    bytes=(out/path).stat().st_size,sha256=digest(out/path)); checkpoint()
                completed=sum(k.startswith('U') and x['status']=='complete' for k,x in ledger['attempts'].items())
                executed=sum(x.get('used_updates',0) for x in ledger['attempts'].values())
                progress=dict(phase='prefix' if prefix else 'main trajectories',total=16,completed=completed,
                    pending=16-completed,failed=0,durably_committed=completed,percent=100*completed/16,
                    executed_updates=executed,reserved_updates=ledger['reserved_updates'],last_unit=key,
                    elapsed_seconds=ledger['elapsed_seconds'],heartbeat_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                    remaining_trajectory_ETA_seconds=(16-completed)*ledger['elapsed_seconds']/completed if completed else None)
                write(out/'progress.json',progress,exclusive=False)
                with (out/'progress.tsv').open('a') as stream:
                    if stream.tell()==0: stream.write('completed\ttotal\texecuted_updates\telapsed_seconds\tlast_unit\n')
                    stream.write(f"{completed}\t16\t{executed}\t{ledger['elapsed_seconds']:.6f}\t{key}\n")
                    stream.flush()
                print(f"{key} {unit['method']} {unit['profile']} {unit['digits']}: accepted; {completed}/16, updates={executed}, wall={result['wall_seconds']:.3f}s",flush=True)
        finally:
            signal.setitimer(signal.ITIMER_REAL,0); checkpoint()
    checked=verify(cfg,out)
    if not (out/'quality.json').exists(): write(out/'quality.json',checked)
    print(__import__('json').dumps(checked,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['design','run','verify'])
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(); cfg=read(ROOT/'config.json')
    if args.command=='design': print(__import__('json').dumps(dict(units=design(cfg),limits=cfg['limits']),indent=2))
    elif args.command=='verify': print(__import__('json').dumps(verify(cfg),indent=2))
    else: run(cfg,args.resume)
