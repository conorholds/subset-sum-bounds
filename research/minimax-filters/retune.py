"""Bounded 17-block retuning with the predecessor's fixed positive reserves."""
import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
import hashlib
import json
import math
from pathlib import Path

import optimize
from filters import scores, moment_logs_from_scores

HERE = Path(__file__).resolve().parent
H = -.25*math.log2(.25)-.75*math.log2(.75)


def parameters_at(eta):
    rho = (eta-.5-.00001-.00001+.000016)/(H-.5)
    if not 0 < 17*rho < .999:
        raise ValueError('density lies outside the unclipped 17-block slice')
    tau = (.5+(H-1+.00003)*rho+(.5+.00002+.000026)/17+.000016)/(1+1/17)
    return {'eta': eta, 'rho': rho, 'tau': tau}


def retune(design):
    sc = scores(design)
    evaluations = 0

    def residual(eta):
        nonlocal evaluations
        row = parameters_at(eta)
        evaluations += 1
        return max(moment_logs_from_scores(sc, row))-(eta-.0001)

    lo, hi = .5179, .51805
    if not residual(lo) < 0 < residual(hi):
        raise ValueError('moment boundary is not bracketed in the retained slice')
    for _ in range(40):
        mid = (lo+hi)/2
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    row = parameters_at(lo)
    eta = Decimal(str(row['eta'])).quantize(Decimal('.00000001'), rounding=ROUND_FLOOR)
    rho = Decimal(str(parameters_at(float(eta))['rho'])).quantize(Decimal('.00000001'), rounding=ROUND_HALF_EVEN)
    tau = Decimal(str(parameters_at(float(eta))['tau'])).quantize(Decimal('.00000001'), rounding=ROUND_CEILING)+Decimal('.0000001')
    d = Decimal('.5')-tau+Decimal('.00002')
    delta = ((d+Decimal('.000026'))/(17*rho)).quantize(Decimal('.0000000001'), rounding=ROUND_CEILING)
    p = {'eta': str(eta), 'tau': str(tau), 'rho': str(rho), 'blocks': 17,
         'deficiency': str(d), 'local_deficiency': str(delta),
         'gamma': str(Decimal('.5')-delta-Decimal('.00003')), 'guard': '.00001',
         'preliminary': str(tau-Decimal('.00001')), 'alpha': str(tau+Decimal('.00001')),
         'epsilon': '.00000001', 'moment_ceiling': str(eta-Decimal('.0001'))}
    rounded_logs = moment_logs_from_scores(sc, p)
    return {'parameters': p, 'bracket': [lo, hi], 'evaluations': evaluations+1,
            'rounded_log2_moments': rounded_logs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--witnesses-output', type=Path, required=True)
    args = parser.parse_args()
    baseline = optimize.BASELINE['witness']
    witnesses = [dict(baseline, id='predecessor_block17')]
    results = []
    for degree in [10, 14]:
        path = HERE/f'discovery-degree{degree}.json'
        design = json.loads(path.read_text())['design']
        tuned = retune(design)
        results.append({'degree': degree, 'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), **tuned})
        witnesses.append({'id': f'degree{degree}_filter_only', 'design': design, 'parameters': baseline['parameters']})
        witnesses.append({'id': f'degree{degree}_retuned', 'design': design, 'parameters': tuned['parameters']})
    report = {'status': 'EXPLORATORY', 'scope': 'A fixed-reserve 17-block slice, not global optimality.',
              'eta_bracket': ['.5179', '.51805'], 'bisections': 40, 'results': results,
              'retuner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    args.witnesses_output.write_text(json.dumps({'scope': 'Conditional minimax filter candidates; independent review pending.',
                                                'witnesses': witnesses}, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
