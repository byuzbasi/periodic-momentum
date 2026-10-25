"""Dense six-coordinate adapter; coefficients and oracle policies reuse v06."""
from pathlib import Path
import sys
import time
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'periodic_momentum_numerics_v06/code'))
from numerics import parameters, method_coefficients, error_vector, norm, serialize
from mpmath import mp


def preparation(cfg, bounds, raw):
    start=time.perf_counter()
    rows=[[mp.mpf(t) for t in line.split()] for line in raw.read_text().splitlines() if line.strip()]
    assert len(rows)==cfg['n'] and all(len(r)==7 for r in rows)
    X=mp.matrix([r[:6] for r in rows]); y=mp.matrix([r[6] for r in rows])
    H=X.T*X/cfg['n']; h=X.T*y/cfg['n']
    star,residual=mp.qr_solve(X,y)
    qr_gradient=H*star-h
    normalized=norm(qr_gradient)/(mp.norm(H)*norm(star)+norm(h))
    assert normalized <= mp.mpf(cfg['tolerances']['reference'])
    p=parameters({'mu':str(bounds['mu_safe']),'L':str(bounds['L_safe'])},cfg['R'])
    spectrum=list(mp.eigsy(H,eigvals_only=True))
    assert all(p['mu']<=v<=p['L'] and (v<=p['a'] or v>=p['b']) for v in spectrum)
    assert p['R']<=p['Rmax'] and 0<p['epsilon_lower']<=p['epsilon_upper']<1
    return dict(H=H,h=h,star=star,p=p,spectrum=spectrum,reference_gradient=qr_gradient,
        normalized_reference_gradient=normalized,residual_norm=residual,
        f_star=norm(X*star-y)**2/(2*cfg['n']),seconds=time.perf_counter()-start)


def pack_setup(obj):
    matrices=['H','h','star','reference_gradient']
    out={k:([[obj[k][i,j] for j in range(obj[k].cols)] for i in range(obj[k].rows)] if k in matrices else v)
         for k,v in obj.items()}
    return serialize(out)


def unpack_setup(obj):
    matrices=['H','h','star','reference_gradient']
    def cast(v):
        if isinstance(v,dict): return {k:cast(x) for k,x in v.items()}
        if isinstance(v,list): return [cast(x) for x in v]
        return mp.mpf(v) if isinstance(v,str) else v
    out=cast(obj)
    for k in matrices: out[k]=mp.matrix(out[k])
    return out


def budget(p,profile):
    if profile=='exact': return mp.mpf(0)
    return p['epsilon_lower']/2 if profile=='certified_mixing' else p['epsilon_upper']*(mp.mpf('0.5') if profile=='scalar_half' else mp.mpf('1.5'))


def state_matrix(s,method,profile,eps,k):
    d=6; Q=mp.zeros(d)
    if profile.startswith('scalar'):
        Q=mp.eye(d)
    elif profile=='certified_mixing':
        for i in range(d): Q[i,(i+1)%d if k%3==1 else (i-1)%d]=-1 if k%3==2 else 1
    alpha,beta=method_coefficients(s['p'],method)
    A=(1+beta)*mp.eye(d)-alpha[k%2]*(mp.eye(d)+eps*Q)*s['H']
    T=mp.zeros(2*d)
    for i in range(d):
        for j in range(d): T[i,j]=A[i,j]
        T[i,d+i]=-beta; T[d+i,i]=1
    return T


def independent_state(s,method,profile,eps,n):
    period=6 if profile=='certified_mixing' else 2
    block=mp.eye(12)
    for k in range(period): block=state_matrix(s,method,profile,eps,k)*block
    prefix=mp.eye(12)
    for k in range(n%period): prefix=state_matrix(s,method,profile,eps,k)*prefix
    initial=mp.matrix(list(-s['star'])+list(-s['star']))
    return prefix*(block**(n//period))*initial


def smoke(s):
    p=s['p']; first=p['alpha'][0]*s['h']
    second=(1+p['m'])*first-p['alpha'][1]*(s['H']*first-s['h'])
    actual=mp.matrix(list(second-s['star'])+list(first-s['star']))
    error=norm(actual-independent_state(s,'CHB','exact',mp.mpf(0),2))/norm(actual)
    assert error<mp.mpf('1e-60')
    return {'gradient_updates':2,'full_trajectories':0,'matrix_relative_error':mp.nstr(error,120)}


def trajectory(cfg,s,unit,counter):
    method,profile=unit['method'],unit['profile']; H,h,star=s['H'],s['h'],s['star']
    eps=budget(s['p'],profile); alpha,beta=method_coefficients(s['p'],method)
    w=mp.zeros(6,1); previous=w.copy(); initial=norm(star)
    initial_gap=(star.T*H*star)[0]/2; initial_state=mp.sqrt(2)*initial
    trace=[]; errors={}; points={}; largest=mp.mpf(0)
    for k in range(cfg['horizon']+1):
        x=w-star; xp=previous-star; gap=(x.T*H*x)[0]/2
        trace.append(dict(k=k,ln_relative_x=mp.log(norm(x)/initial),
            ln_relative_f=mp.log(gap/initial_gap),ln_relative_state=mp.log(norm(list(x)+list(xp))/initial_state)))
        if k in cfg['checkpoints']:
            z=mp.matrix(list(x)+list(xp)); reference=independent_state(s,method,profile,eps,k)
            errors[k]=norm(z-reference)/norm(z); points[k]=list(z)
        if k==cfg['horizon']: break
        counter['updates']+=1
        g=H*w-h
        mapped='mixing_half' if profile=='certified_mixing' else profile
        e=mp.matrix(error_vector(list(g),eps,mapped,k))
        if eps: largest=max(largest,max(mp.mpf(0),norm(e)/(eps*norm(g))-1))
        else: assert norm(e)==0
        following=(1+beta)*w-beta*previous-alpha[k%2]*(g+e)
        previous,w=w,following
    accepted=max(errors.values())<=mp.mpf(cfg['tolerances']['matrix_relative']) and largest<=mp.mpf(cfg['tolerances']['oracle_relative'])
    return serialize(dict(unit=unit,accepted=bool(accepted),epsilon=eps,
        trace=trace,checkpoints=points,matrix_relative_errors=errors,max_oracle_violation=largest,
        final_w=list(w),initial_distance=initial,initial_objective_gap=initial_gap,
        gradient_updates=cfg['horizon'],trajectory_executions=1))
