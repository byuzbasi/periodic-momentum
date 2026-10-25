"""Frozen deterministic CHB/HB pilot mathematics; no optimizer or SDP solver."""
from mpmath import mp


def number(value):
    parts = str(value).split("/")
    return mp.mpf(parts[0]) / mp.mpf(parts[1]) if len(parts) == 2 else mp.mpf(parts[0])


def norm(v):
    return mp.sqrt(mp.fsum(t*t for t in v))


def parameters(cfg, gap):
    mu, L, R = number(cfg["mu"]), number(cfg["L"]), number(gap)
    k = L / mu
    q = (mp.sqrt(k)-1)/(mp.sqrt(k)+1)
    chi = (L+mu)/(L-mu)
    A, B = mp.sqrt(chi**2-R**2), mp.sqrt(chi**2-1)
    m = (1-R**2)/(A+B)**2
    a, b = (L+mu-R*(L-mu))/2, (L+mu+R*(L-mu))/2
    a0, a1 = (1+m)/a, (1+m)/b
    margin = q*q-m
    D, Br, V = a0*a1*L*(L-mu), a0*a1*L*L, margin**2/q**2
    upper = 2*V/(D+mp.sqrt(D*D+4*Br*V))
    dminus = mp.mpf(3)/4 * (2*q*q*mp.sqrt(k)/(k+1))
    T = mp.sqrt(4*k*k+2)
    K, G = mp.sqrt((4*k*k+3)**2+2), 4*k*(T+k)
    c, C = min(mp.mpf(1), dminus*dminus/(4*K*G)), 4*q*q/(k-1)
    return dict(mu=mu,L=L,R=R,kappa=k,q=q,m=m,a=a,b=b,alpha=[a0,a1],
                margin=margin,epsilon_upper=upper,epsilon_lower=c*R**4,
                coarse_upper=C*R**4,c=c,C=C,Rmax=min(mp.mpf('0.5'),(2*C)**(-mp.mpf('0.25'))),
                limit=(k-1)/(4*(k+1)**2),original_m=((A-B)/mp.sqrt(1-R**2))**2,
                lambdas=[mu,a,b,L],D=D,Br=Br)


def method_coefficients(p, method):
    if method == "CHB":
        return p["alpha"], p["m"]
    assert method == "HB"
    return [4/(mp.sqrt(p["L"])+mp.sqrt(p["mu"]))**2]*2, p["q"]**2


def scale_values(cfg, unit, digits):
    with mp.workdps(digits):
        p = parameters(cfg, unit["R"])
        assert 0 < p["R"] <= p["Rmax"]
        assert p["epsilon_lower"] <= p["epsilon_upper"] <= p["coarse_upper"] < 1
        tests = []
        for factor in [mp.mpf('0.5'),mp.mpf(1),mp.mpf('1.5')]:
            eps = factor*p["epsilon_upper"]
            c0,c1 = [1+p["m"]-a*p["L"]*(1+eps) for a in p["alpha"]]
            M = mp.matrix([[c1*c0-p["m"],-p["m"]*c1],[c0,-p["m"]]])
            root_rate = mp.sqrt(max(abs(v) for v in mp.eig(M,left=False,right=False)))
            u = p["D"]*eps+p["Br"]*eps*eps
            closed_rate = mp.sqrt(p["m"]+u/2+mp.sqrt(p["m"]*u+u*u/4))
            tests.append(dict(factor=factor,rate=root_rate,closed_rate=closed_rate,
                              agreement=abs(root_rate-closed_rate),q_difference=root_rate-p["q"]))
        return dict(parameters=p,ratio=p["epsilon_upper"]/p["R"]**4,
                    relative_limit_error=p["epsilon_upper"]/p["R"]**4/p["limit"]-1,
                    coefficient_identity_error=abs(p["m"]-p["original_m"]),rate_checks=tests)


def error_vector(gradient, eps, profile, k):
    if profile == "exact":
        return [mp.mpf(0)]*len(gradient)
    if profile.startswith("scalar"):
        return [eps*t for t in gradient]
    assert profile == "mixing_half"
    if k % 3 == 1:
        rotated = gradient[1:]+gradient[:1]
    else:
        rotated = gradient[-1:]+gradient[:-1]
    sign = -1 if k % 3 == 2 else 1
    return [sign*eps*t for t in rotated]


def transition(p, method, profile, eps, k):
    """Independent matrix realization, not a call to the gradient-error loop."""
    d = len(p["lambdas"])
    Q = mp.zeros(d)
    for i in range(d):
        if profile.startswith("scalar"):
            Q[i,i] = 1
        elif profile == "mixing_half":
            Q[i,(i+1)%d if k%3 == 1 else (i-1)%d] = -1 if k%3 == 2 else 1
    alpha, beta = method_coefficients(p, method)
    A = mp.eye(d)*(1+beta)-alpha[k%2]*(mp.eye(d)+eps*Q)*mp.diag(p["lambdas"])
    T = mp.zeros(2*d)
    for i in range(d):
        for j in range(d):
            T[i,j] = A[i,j]
        T[i,d+i] = -beta
        T[d+i,i] = 1
    return T


def matrix_state(p, method, profile, eps, n, x0):
    period = 6 if profile == "mixing_half" else 2
    block = mp.eye(2*len(x0))
    for k in range(period):
        block = transition(p,method,profile,eps,k)*block
    prefix = mp.eye(2*len(x0))
    for k in range(n%period):
        prefix = transition(p,method,profile,eps,k)*prefix
    return list(prefix*(block**(n//period))*mp.matrix(x0+x0))


def mode_rate(p, method, eps):
    """Two-step spectral rate on the four modes for identity errors only."""
    alpha,beta = method_coefficients(p,method)
    rates = []
    for lam in p["lambdas"]:
        c0,c1 = [1+beta-a*lam*(1+eps) for a in alpha]
        trace = c1*c0-2*beta
        discriminant = mp.sqrt(mp.mpc(trace*trace-4*beta*beta))
        rates.append(mp.sqrt(max(abs((trace+discriminant)/2),abs((trace-discriminant)/2))))
    return max(rates)


def trajectory(cfg, unit, method, digits):
    with mp.workdps(digits):
        p = parameters(cfg,unit["R"])
        profile = unit["profile"]
        eps = number(cfg["profiles"][profile]["factor"])*p["epsilon_upper"]
        alpha,beta = method_coefficients(p,method)
        x0 = [number(v) for v in cfg["initial_x"]]
        x,previous = x0[:],x0[:]
        initial,initial_state = norm(x0),norm(x0+x0)
        logs_x,logs_z = [mp.mpf(0)],[mp.mpf(0)]
        max_error,checkpoints = mp.mpf(0),{}
        for k in range(cfg["horizon"]):
            h = [lam*xi for lam,xi in zip(p["lambdas"],x)]
            e = error_vector(h,eps,profile,k)
            ratio = norm(e)/norm(h) if norm(h) else mp.mpf(0)
            if eps:
                max_error = max(max_error,max(mp.mpf(0),ratio/eps-1))
            else:
                assert norm(e) == 0
            following = [(1+beta)*xi-beta*xm-alpha[k%2]*(hi+ei)
                         for xi,xm,hi,ei in zip(x,previous,h,e)]
            previous,x = x,following
            logs_x.append(mp.log(norm(x)/initial))
            logs_z.append(mp.log(norm(x+previous)/initial_state))
            if k+1 in cfg["checkpoints"]:
                checkpoints[k+1] = x+previous
        return dict(parameters=p,epsilon=eps,log_x=logs_x,log_z=logs_z,
                    max_error_violation=max_error,checkpoints=checkpoints,
                    late_window_factor=mp.exp((logs_z[-1]-logs_z[cfg["horizon"]//2])/(cfg["horizon"]//2)),
                    mode_matrix_rate=None if profile == "mixing_half" else mode_rate(p,method,eps))


def serialize(obj):
    if isinstance(obj,dict): return {str(k):serialize(v) for k,v in obj.items()}
    if isinstance(obj,list): return [serialize(v) for v in obj]
    if isinstance(obj,mp.mpf): return mp.nstr(obj,120)
    return obj


def evaluate(cfg,unit):
    with mp.workdps(max(cfg["precisions"])):
        low_dps,high_dps = cfg["precisions"]
        if unit["kind"] == "scale":
            low,high = [scale_values(cfg,unit,d) for d in cfg["precisions"]]
            precision_error = abs(low["ratio"]/high["ratio"]-1)
            checks = high["rate_checks"]
            tol = number(cfg["tolerances"]["scale_relative"])
            accepted = (precision_error <= tol and high["coefficient_identity_error"] <= tol
                        and max(t["agreement"] for t in checks) <= tol
                        and abs(checks[1]["q_difference"]) <= tol
                        and checks[0]["q_difference"] < 0 < checks[2]["q_difference"])
            return serialize(dict(unit=unit,accepted=bool(accepted),precision_relative_error=precision_error,
                                  results=high,trajectory_executions=0,gradient_updates=0))
        methods,internal = {},{}
        for method in cfg["methods"]:
            low = trajectory(cfg,unit,method,low_dps)
            high = trajectory(cfg,unit,method,high_dps)
            diff = max(abs(a-b) for key in ["log_x","log_z"] for a,b in zip(low[key],high[key]))
            matrix_errors = {}
            x0 = [number(v) for v in cfg["initial_x"]]
            for n,state in high["checkpoints"].items():
                direct = matrix_state(high["parameters"],method,unit["profile"],high["epsilon"],n,x0)
                matrix_errors[n] = norm([a-b for a,b in zip(state,direct)])/norm(state)
            accepted = (diff <= number(cfg["tolerances"]["trajectory_log"])
                        and max(matrix_errors.values()) <= number(cfg["tolerances"]["matrix_relative"])
                        and max(low["max_error_violation"],high["max_error_violation"])
                            <= number(cfg["tolerances"]["error_relative"]))
            internal[method] = high
            methods[method] = dict(accepted=bool(accepted),precision_max_log_error=diff,
                                  max_relative_error_violation=max(low["max_error_violation"],high["max_error_violation"]),
                                  matrix_relative_errors=matrix_errors,
                                  late_window_factor=high["late_window_factor"],mode_matrix_rate=high["mode_matrix_rate"],
                                  checkpoints=high["checkpoints"],
                                  trace=[dict(k=k,ln_relative_x=float(lx),ln_relative_state=float(lz),
                                              log10_x_over_qk=float((lx-k*mp.log(high["parameters"]["q"]))/mp.log(10)))
                                         for k,(lx,lz) in enumerate(zip(high["log_x"],high["log_z"]))])
        differences = [a-b for a,b in zip(internal["CHB"]["log_x"],internal["HB"]["log_x"])]
        last_nonwin = max(k for k,d in enumerate(differences) if d >= -number(cfg["tolerances"]["crossover_log"]))
        crossover = last_nonwin+1 if last_nonwin < cfg["horizon"] else None
        p = internal["CHB"]["parameters"]
        return serialize(dict(unit=unit,accepted=all(v["accepted"] for v in methods.values()),parameters=p,
                              epsilon=internal["CHB"]["epsilon"],methods=methods,
                              terminal_chb_to_hb_ratio=mp.exp(differences[-1]),
                              observed_suffix_win_start=crossover,
                              trajectory_executions=4,gradient_updates=4*cfg["horizon"]))
