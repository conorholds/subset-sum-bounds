"""Bounded fixed-filter frontier search after eliminating auxiliary coverage budgets."""
import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_HALF_EVEN
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'
sys.path.insert(0, str(CORE))
from filters import scores, moment_logs_from_scores
from verify import load_witnesses

PREDECESSOR = 'f25a303525df636c079e6a7c574527cccf56e611'
ENTROPY = -.25*math.log2(.25)-.75*math.log2(.75)
RESERVES = {'deficiency': '.00002', 'coverage': '.000026', 'representations': '.00003',
            'filtered_lists': '.000016', 'visited_pairs': '.00001', 'block_space': '.001',
            'moments': '.0001', 'guard': '.00001', 'epsilon': '.00000001', 'assembly': '.00001'}


def balanced_rho(eta, blocks):
    if type(blocks) is not int or not 1 <= blocks <= 17 or not math.isfinite(eta):
        raise ValueError('finite eta and a block count from 1 through 17 required')
    return min((eta-.5-float(RESERVES['guard'])-float(RESERVES['visited_pairs'])
                +float(RESERVES['filtered_lists']))/(ENTROPY-.5),
               (1-float(RESERVES['block_space']))/blocks)


def screen(eta, blocks, moment_floor):
    rho = balanced_rho(eta, blocks)
    if not 0 < rho < .5:
        return None
    hrho = -rho*math.log2(rho)-(1-rho)*math.log2(1-rho)
    if eta <= hrho+rho+float(RESERVES['guard']):
        return None
    coverage, representations = float(RESERVES['coverage']), float(RESERVES['representations'])
    constant = .5+float(RESERVES['deficiency'])+coverage
    lists = (.5+(ENTROPY-1+representations)*rho+constant/blocks
             +float(RESERVES['filtered_lists']))/(1+1/blocks)
    pairs = (1+(2*ENTROPY-1.5+representations)*rho+constant/blocks-eta
             +float(RESERVES['guard'])+float(RESERVES['visited_pairs']))/(1+1/blocks)
    moment = moment_floor(eta)
    tau = max(lists, pairs, moment, rho+.251)
    if tau >= .49999:
        return None
    d = .5-tau+float(RESERVES['deficiency'])
    delta = (d+coverage)/(blocks*rho)
    if .5-delta-representations <= 0:
        return None
    return {'blocks': blocks, 'eta': eta, 'rho': rho, 'tau': tau,
            'caps': {'filtered_lists': lists, 'visited_pairs': pairs, 'moments': moment}}


def search(blocks=tuple(range(1, 18)), coarse_points=401, refinement_points=201):
    if (not blocks or len(set(blocks)) != len(blocks)
            or any(type(j) is not int or not 1 <= j <= 17 for j in blocks)):
        raise ValueError('distinct block counts from 1 through 17 required')
    if (type(coarse_points) is not int or not 3 <= coarse_points <= 401
            or type(refinement_points) is not int or not 3 <= refinement_points <= 201):
        raise ValueError('grid sizes exceed the bounded campaign')
    started = time.perf_counter()
    design = next(w for w in load_witnesses() if w['id'] == 'combined')['design']
    sc = scores(design)
    if min(sc[1]) <= 0:
        raise ValueError('bisection requires positive sampled v scores')
    moment_evaluations = 0

    @lru_cache(None)
    def moment_floor(eta):
        def evaluate(tau):
            nonlocal moment_evaluations
            moment_evaluations += 1
            return max(moment_logs_from_scores(sc, {'eta': eta, 'tau': tau}))
        lo, hi = .45, .499999
        if evaluate(hi) >= eta-float(RESERVES['moments']):
            return math.inf
        for _ in range(40):
            mid = (lo+hi)/2
            if evaluate(mid) >= eta-float(RESERVES['moments']):
                lo = mid
            else:
                hi = mid
        return hi+1e-9

    rows = []
    for count in blocks:
        candidates = [screen(float(eta), count, moment_floor) for eta in np.linspace(.510, .550, coarse_points)]
        best = min((row for row in candidates if row), key=lambda row: row['tau'], default=None)
        if best is None:
            raise ValueError(f'no feasible screened point for block count {count}')
        radius = .040/(coarse_points-1)
        for _ in range(2):
            candidates = [screen(float(eta), count, moment_floor)
                          for eta in np.linspace(max(.510, best['eta']-radius),
                                                 min(.550, best['eta']+radius), refinement_points)]
            best = min((row for row in candidates if row), key=lambda row: row['tau'])
            radius *= 2/(refinement_points-1)
        rows.append(best)
    return {
        'status': 'EXPLORATORY', 'scope': 'Bounded floating-point screens; exact witnesses require interval replay.',
        'predecessor_commit': PREDECESSOR, 'design_id': 'combined', 'reserves': RESERVES,
        'bounds': {'blocks': list(blocks), 'eta': ['.510', '.550'], 'coarse_points': coarse_points,
                   'refinement_points': refinement_points, 'refinement_rounds': 2,
                   'moment_bisections': 40, 'grid': 128},
        'parameter_screens': len(blocks)*(coarse_points+2*refinement_points),
        'moment_evaluations': moment_evaluations,
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
        'inputs_sha256': {name: hashlib.sha256((CORE/name).read_bytes()).hexdigest()
                          for name in ['filters.py', 'verify.py', 'certify.py', 'SOURCE.json', 'witnesses.json']},
        'search_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'frontier': rows, 'seconds': time.perf_counter()-started}


def witness(row, design, rounded=False):
    tick = Decimal('.00001' if rounded else '.0000001')
    eta = Decimal(str(row['eta'])).quantize(tick, rounding=ROUND_HALF_EVEN)
    rho = Decimal(str(row['rho'])).quantize(tick, rounding=ROUND_HALF_EVEN)
    tau = Decimal(str(row['tau'])).quantize(tick, rounding=ROUND_CEILING)
    if not rounded:
        tau += tick
    d = Decimal('.5')-tau+Decimal(RESERVES['deficiency'])
    delta = ((d+Decimal(RESERVES['coverage']))/(row['blocks']*rho)).quantize(
        Decimal('.000001' if rounded else '.000000001'), rounding=ROUND_CEILING)
    parameters = {'eta': str(eta), 'tau': str(tau), 'rho': str(rho), 'blocks': row['blocks'],
                  'deficiency': str(d), 'local_deficiency': str(delta),
                  'gamma': str(Decimal('.5')-delta-Decimal(RESERVES['representations'])),
                  'guard': RESERVES['guard'], 'preliminary': str(tau-Decimal(RESERVES['assembly'])),
                  'alpha': str(tau+Decimal(RESERVES['assembly'])), 'epsilon': RESERVES['epsilon'],
                  'moment_ceiling': str(eta-Decimal(RESERVES['moments']))}
    return {'id': f'block{row["blocks"]}_{"rounded" if rounded else "retuned"}',
            'design': design, 'parameters': parameters}


def retained_witnesses(result):
    design = next(w for w in load_witnesses() if w['id'] == 'combined')['design']
    rows = {row['blocks']: row for row in result['frontier']}
    return {'scope': 'Conditional fixed-filter block-count frontier; independent mathematical review pending.',
            'witnesses': [witness(rows[16], design), witness(rows[17], design),
                          witness(rows[17], design, rounded=True)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='smaller grids for smoke testing, not the retained campaign')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--witnesses-output', type=Path)
    args = parser.parse_args()
    result = search(coarse_points=41 if args.quick else 401,
                    refinement_points=21 if args.quick else 201)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    if args.witnesses_output:
        args.witnesses_output.write_text(json.dumps(retained_witnesses(result), indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'parameter_screens', 'moment_evaluations', 'seconds']}))


if __name__ == '__main__':
    main()
