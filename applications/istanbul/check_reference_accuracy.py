"""Audit QR-reference metrics using saved matrices/states only; no solver."""
from pathlib import Path
import argparse
import json
from mpmath import mp

ROOT = Path(__file__).resolve().parent


def audit():
    output = {}
    with mp.workdps(220):
        for digits in [120, 160]:
            setup = json.loads((ROOT / f'accepted/setup_{digits}.json').read_text())
            H, h, reference = [mp.matrix(setup[key]) for key in ['H', 'h', 'star']]
            mu = mp.mpf(setup['p']['mu'])
            residual = H * reference - h
            recorded = mp.matrix(setup['reference_gradient'])
            rn = mp.norm(residual)
            initial_energy = (reference.T * H * reference)[0] / 2
            initial_gap_bound = (rn * mp.norm(reference) + rn**2 / (2 * mu)) / initial_energy
            initial_distance_bound = rn / (mu * mp.norm(reference))
            maxima = {key: mp.mpf(0) for key in [
                'linear_correction_relative_to_energy',
                'recorded_residual_linear_correction',
                'gap_relative_to_exact_minimizer_bound',
                'distance_relative_to_exact_minimizer_bound',
                'normalized_gap_relative_difference_bound',
                'normalized_distance_relative_difference_bound']}
            count = 0
            for path in sorted((ROOT / 'accepted/records').glob('U*.json')):
                record = json.loads(path.read_text())
                if record['unit']['digits'] != digits:
                    continue
                states = {0: list(-reference)}
                for k, vector in record['checkpoints'].items():
                    states[int(k)] = vector[:7]
                    states[int(k) - 1] = vector[7:]
                for values in states.values():
                    difference = mp.matrix(values)
                    energy = (difference.T * H * difference)[0] / 2
                    norm = mp.norm(difference)
                    assert energy > 0 and norm > 0 and mu > 0
                    gap_bound = (rn * norm + rn**2 / (2 * mu)) / energy
                    distance_bound = rn / (mu * norm)
                    values = [abs((residual.T * difference)[0]) / energy,
                              abs((recorded.T * difference)[0]) / energy,
                              gap_bound, distance_bound,
                              (gap_bound + initial_gap_bound) / (1 - initial_gap_bound),
                              (distance_bound + initial_distance_bound) / (1 - initial_distance_bound)]
                    for key, value in zip(maxima, values):
                        maxima[key] = max(maxima[key], value)
                    count += 1
            assert count == 72
            output[str(digits)] = dict(saved_states_checked=count,
                saved_state_indices=[0, 1, 2, 15, 16, 63, 64, 511, 512],
                maxima={key: mp.nstr(value, 30) for key, value in maxima.items()})
    return dict(accepted=True,scope='Floating-point algebra on saved matrices and states; no interval certificate',
        diagnostic_digits=220,precisions=output,new_optimizer_updates=0,new_QR_solves=0,new_eigensolves=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), 'Refusing to overwrite an existing diagnostic'
    result = audit()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print('Reference metrics checked at 144 saved states; no optimizer or solver calls.')
