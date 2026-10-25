"""GD/NAG-Q adapter using the saved quadratic setup and scalar-error helpers."""
from pathlib import Path
import json
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT.parents[1]
sys.path.insert(0, str(PUBLIC / 'scripts'))
from common import kernels, read, require, write, digest, frozen_check
numeric, dense = kernels()
from mpmath import mp


def units(cfg):
    return [dict(id=f'U{i+1:02d}', case=c['id'], profile=p, method=m, digits=d)
            for i, (c, p, m, d) in enumerate((c, p, m, d) for c in cfg['cases']
                for p in cfg['profiles'] for m in cfg['methods'] for d in cfg['precisions'])]


def setup(cfg, unit):
    case = next(c for c in cfg['cases'] if c['id'] == unit['case'])
    if case['kind'] == 'synthetic':
        p = numeric.parameters(read(PUBLIC / 'configs/synthetic.json'), case['R'])
        H = mp.diag(p['lambdas']); h = mp.zeros(4, 1); star = mp.zeros(4, 1)
        initial = mp.matrix([1]*4)
    else:
        saved = read(PUBLIC / f"data/accepted/yacht/setup_{unit['digits']}.json")
        s = dense.unpack_setup(saved)
        H, h, star, p = [s[k] for k in ['H', 'h', 'star', 'p']]
        initial = mp.zeros(6, 1)
    epsilon = mp.mpf(0) if unit['profile'] == 'exact' else p['epsilon_upper']/2
    return H, h, star, initial, p, epsilon


def coefficients(p, method):
    if method == 'GD':
        return 2/(p['L']+p['mu']), mp.mpf(0)
    require(method == 'NAG-Q', 'Unapproved method')
    root = mp.sqrt(3*p['kappa']+1)
    return 4/(3*p['L']+p['mu']), (root-2)/(root+2)


def affine_transition(H, h, star, alpha, beta, epsilon):
    """Independent centered affine state matrix, including reference residual."""
    d = H.rows
    Q = mp.eye(d)-alpha*(1+epsilon)*H
    T = mp.zeros(2*d+1)
    residual = H*star-h
    for i in range(d):
        for j in range(d):
            T[i, j] = (1+beta)*Q[i, j]
            T[i, d+j] = -beta*Q[i, j]
        T[i, 2*d] = -alpha*(1+epsilon)*residual[i]
        T[d+i, i] = 1
    T[2*d, 2*d] = 1
    return T


def trajectory(cfg, unit, counter, steps=None):
    with mp.workdps(unit['digits']):
        H, h, star, x, p, epsilon = setup(cfg, unit)
        initial = x.copy(); previous = x.copy(); initial_error = initial-star
        alpha, beta = coefficients(p, unit['method'])
        initial_norm = numeric.norm(initial_error)
        initial_gap = (initial_error.T*H*initial_error)[0]/2
        initial_state = mp.sqrt(2)*initial_norm
        T = affine_transition(H, h, star, alpha, beta, epsilon)
        z0 = mp.matrix(list(initial_error)+list(initial_error)+[1])
        trace, points, checks = [], {}, {}
        oracle_violation = mp.mpf(0)
        horizon = cfg['horizon'] if steps is None else steps
        for k in range(horizon+1):
            error, prev_error = x-star, previous-star
            distance = numeric.norm(error)
            gap = (error.T*H*error)[0]/2
            require(distance > 0 and gap > 0, 'Unexpected exact zero; inspect before logging')
            trace.append(dict(k=k, ln_relative_x=mp.log(distance/initial_norm),
                ln_relative_f=mp.log(gap/initial_gap),
                ln_relative_state=mp.log(numeric.norm(list(error)+list(prev_error))/initial_state)))
            if k in cfg['checkpoints']:
                state = mp.matrix(list(error)+list(prev_error))
                reference = (T**k)*z0
                checks[k] = numeric.norm(state-reference[:2*H.rows, :])/numeric.norm(state)
                points[k] = list(state)
            if k == horizon:
                break
            y = x+beta*(x-previous)
            counter.charge()
            gradient = H*y-h
            perturbation = mp.matrix(numeric.error_vector(list(gradient), epsilon, unit['profile'], k))
            if epsilon and numeric.norm(gradient):
                oracle_violation = max(oracle_violation, numeric.norm(perturbation)/(epsilon*numeric.norm(gradient))-1)
            else:
                require(numeric.norm(perturbation) == 0, 'Nonzero exact/zero-gradient error')
            following = y-alpha*(gradient+perturbation)
            previous, x = x, following
        max_matrix = max(checks.values())
        monotone = (unit['method'] != 'GD' or unit['profile'] != 'exact' or all(
            b['ln_relative_f'] <= a['ln_relative_f']+mp.mpf('1e-60')
            for a, b in zip(trace, trace[1:])))
        accepted = (max_matrix <= mp.mpf(cfg['tolerances']['matrix_relative']) and
            oracle_violation <= mp.mpf(cfg['tolerances']['oracle_relative']) and monotone)
        return numeric.serialize(dict(unit=unit, accepted=bool(accepted), steps=horizon,
            gradient_updates=horizon, alpha=alpha, beta=beta, epsilon=epsilon,
            initial_distance=initial_norm, initial_gap=initial_gap, trace=trace,
            checkpoints=points, matrix_relative_errors=checks,
            max_matrix_relative_error=max_matrix, max_oracle_violation=oracle_violation,
            exact_GD_objective_monotone=bool(monotone)))
