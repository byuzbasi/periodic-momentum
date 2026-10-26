"""Seven-coordinate adapter using the unchanged coefficients/error helpers."""
from pathlib import Path
import sys
import time
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT.parents[1]
sys.path.insert(0,str(PUBLIC/'scripts'))
from common import read,write,require,digest,kernels
numeric,_ = kernels()
from mpmath import mp
RUN_DIRECTORY=ROOT/'accepted'

def unpack_setup(obj):
    matrices=['H','h','star','reference_gradient']
    def cast(v):
        if isinstance(v,dict): return {k:cast(x) for k,x in v.items()}
        if isinstance(v,list): return [cast(x) for x in v]
        return mp.mpf(v) if isinstance(v,str) else v
    out=cast(obj)
    for k in matrices: out[k]=mp.matrix(out[k])
    return out


def pack(value,digits):
    """Retain the full approved working precision; old 120-digit writer is unchanged."""
    if isinstance(value,mp.matrix):
        return [[pack(value[i,j],digits) for j in range(value.cols)] for i in range(value.rows)]
    if isinstance(value,dict): return {str(k):pack(v,digits) for k,v in value.items()}
    if isinstance(value,list): return [pack(v,digits) for v in value]
    if isinstance(value,mp.mpf): return mp.nstr(value,digits)
    return value


def units(cfg):
    return [dict(id=f'U{i+1:02d}',profile=p,method=m,digits=d)
        for i,(p,m,d) in enumerate((p,m,d) for p in cfg['profiles']
                                 for m in cfg['methods'] for d in cfg['precisions'])]


def stationary_coefficients(p,method):
    if method == 'GD':
        return 2/(p['L']+p['mu']), mp.mpf(0)
    require(method == 'NAG-Q', 'Unapproved method')
    root = mp.sqrt(3*p['kappa']+1)
    return 4/(3*p['L']+p['mu']), (root-2)/(root+2)


def tuning(p,method,k):
    if method in ['CHB','HB']:
        alpha,beta=numeric.method_coefficients(p,method)
        return alpha[k%2],beta,mp.mpf(0)
    alpha,beta=stationary_coefficients(p,method)
    return alpha,beta,beta if method=='NAG-Q' else mp.mpf(0)


def preparation(cfg,digits,charge):
    with mp.workdps(digits):
        started=time.perf_counter(); raw=read(ROOT/'input_tokens.json')
        X=mp.matrix([[mp.mpf(t) for t in row] for row in raw['X']])
        y=mp.matrix([mp.mpf(t) for t in raw['y']])
        require(X.rows==cfg['n'] and X.cols==cfg['p'] and y.rows==cfg['n'],'Input dimensions')
        H=X.T*X/cfg['n']; h=X.T*y/cfg['n']
        charge('QR',digits)
        star,residual=mp.qr_solve(X,y)
        reference_gradient=H*star-h
        normalized=numeric.norm(reference_gradient)/(mp.norm(H)*numeric.norm(star)+numeric.norm(h))
        require(normalized<=mp.mpf(cfg['tolerances']['reference']),'QR reference gate failed')
        b=cfg['bounds']; p=numeric.parameters({'mu':b['mu_safe'],'L':b['L_safe']},b['R_use'])
        charge('eigen',digits)
        spectrum=list(mp.eigsy(H,eigvals_only=True))
        require(all(p['mu']<=v<=p['L'] and (v<=p['a'] or v>=p['b']) for v in spectrum),'Frozen spectral enclosure failed')
        require(p['R']<=p['Rmax'] and 0<p['epsilon_lower']<=p['epsilon_upper']<1,'Frozen theorem range')
        require(numeric.norm(star)>0,'Zero initial reference distance')
        return pack(dict(H=H,h=h,star=star,p=p,spectrum=spectrum,
            reference_gradient=reference_gradient,normalized_reference_gradient=normalized,
            residual_norm=residual,f_star=numeric.norm(X*star-y)**2/(2*cfg['n']),
            seconds=time.perf_counter()-started,working_digits=digits,stored_digits=digits),digits)


def affine_transition(s,method,profile,eps,k):
    """Independent centered affine dynamics, including the QR reference residual."""
    H,h,star=s['H'],s['h'],s['star']; d=H.rows
    Q=mp.zeros(d)
    if profile=='certified_mixing':
        for i in range(d): Q[i,(i+1)%d if k%3==1 else (i-1)%d]=-1 if k%3==2 else 1
    else: require(profile=='exact','Unapproved error profile')
    R=mp.eye(d)+eps*Q
    alpha,beta,eta=tuning(s['p'],method,k)
    RH=R*H; residual=R*(H*star-h)
    A=(1+beta)*mp.eye(d)-alpha*(1+eta)*RH
    B=-beta*mp.eye(d)+alpha*eta*RH
    T=mp.zeros(2*d+1)
    for i in range(d):
        for j in range(d): T[i,j]=A[i,j]; T[i,d+j]=B[i,j]
        T[i,2*d]=-alpha*residual[i]; T[d+i,i]=1
    T[2*d,2*d]=1
    return T


def independent_operators(s,method,profile,eps):
    period=6 if profile=='certified_mixing' else 2
    transitions=[affine_transition(s,method,profile,eps,k) for k in range(period)]
    prefixes=[mp.eye(2*s['H'].rows+1)]
    for T in transitions: prefixes.append(T*prefixes[-1])
    return period,prefixes[-1],prefixes


def trajectory(cfg,unit,counter,steps=None):
    with mp.workdps(unit['digits']):
        s=unpack_setup(read(RUN_DIRECTORY/f"setup_{unit['digits']}.json"))
        H,h,star,p=[s[k] for k in ['H','h','star','p']]; d=H.rows
        w=mp.zeros(d,1); previous=w.copy(); initial=numeric.norm(star)
        initial_gap=(star.T*H*star)[0]/2; initial_state=mp.sqrt(2)*initial
        eps=mp.mpf(0) if unit['profile']=='exact' else p['epsilon_lower']/2
        period,block,prefixes=independent_operators(s,unit['method'],unit['profile'],eps)
        z0=mp.matrix(list(-star)+list(-star)+[1])
        trace,checks,points=[],{},{}; largest=mp.mpf(0)
        horizon=cfg['horizon'] if steps is None else steps
        for k in range(horizon+1):
            x,xp=w-star,previous-star; distance=numeric.norm(x); gap=(x.T*H*x)[0]/2
            require(distance>0 and gap>0,'Unexpected exact zero; inspect before logging')
            trace.append(dict(k=k,ln_relative_x=mp.log(distance/initial),ln_relative_f=mp.log(gap/initial_gap),
                ln_relative_state=mp.log(numeric.norm(list(x)+list(xp))/initial_state)))
            if k in cfg['checkpoints']:
                actual=mp.matrix(list(x)+list(xp))
                reference=prefixes[k%period]*(block**(k//period))*z0
                checks[k]=numeric.norm(actual-reference[:2*d,:])/numeric.norm(actual)
                points[k]=list(actual)
            if k==horizon: break
            alpha,beta,eta=tuning(p,unit['method'],k)
            y=w+eta*(w-previous)
            counter.charge()
            gradient=H*y-h
            mapped='mixing_half' if unit['profile']=='certified_mixing' else unit['profile']
            perturbation=mp.matrix(numeric.error_vector(list(gradient),eps,mapped,k))
            if eps and numeric.norm(gradient):
                largest=max(largest,numeric.norm(perturbation)/(eps*numeric.norm(gradient))-1)
            else: require(numeric.norm(perturbation)==0,'Nonzero exact or zero-gradient error')
            following=w+beta*(w-previous)-alpha*(gradient+perturbation)
            previous,w=w,following
        maximum=max(checks.values())
        monotone=unit['method']!='GD' or unit['profile']!='exact' or all(
            b['ln_relative_f']<=a['ln_relative_f']+mp.mpf('1e-60') for a,b in zip(trace,trace[1:]))
        accepted=maximum<=mp.mpf(cfg['tolerances']['matrix_relative']) and largest<=mp.mpf(cfg['tolerances']['oracle_relative']) and monotone
        return pack(dict(unit=unit,accepted=bool(accepted),steps=horizon,gradient_updates=horizon,
            epsilon=eps,initial_distance=initial,initial_gap=initial_gap,trace=trace,checkpoints=points,
            matrix_relative_errors=checks,max_matrix_relative_error=maximum,max_oracle_violation=largest,
            exact_GD_objective_monotone=bool(monotone),final_w=list(w)),unit['digits'])
