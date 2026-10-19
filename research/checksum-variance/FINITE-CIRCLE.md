# An effective finite-circle moment threshold

For the exact degree-twenty asymmetric candidate in
[FINITE-INPUT.json](FINITE-INPUT.json), the existing continuous log-moment
ceiling `chi=0.51836094` transfers to the finite-circle ceiling

```text
xi=0.51841070 for every integer P >= 2^351.
```

Consequently `n>=681` suffices for **this transfer alone** whenever
`|log2 P-eta n|<=2` with `eta=0.51846094`. This is not a full algorithm
cutoff. Coverage, cap probabilities, bounded random draws and numerical
precision have other requirements. No practical crossover is claimed.

[finite_circle.py](finite_circle.py) computes the conservative constants
with exact rational arithmetic and retains them in [FINITE.json](FINITE.json).
It assumes the pinned continuous ceiling; it does not re-prove that ceiling.
The following argument supplies the analytic implication.

## 1. Finite inner scores

Use raw real trigonometric polynomials `r,s` and `F=2^r,G=2^s`. Write

```text
B=||r||_infinity+||s||_infinity,
A=||r'||_infinity+||s'||_infinity.
```

The calculator replaces these norms by the upper bounds obtained by summing
absolute cosine and sine coefficients, with harmonic index factors for
derivatives. Put `h=2pi/P`. Fix a grid difference `x` and define

```text
l(w)=r(x+w)+s(w), C(x)=integral 2^l(w),
T(x)=integral 2^l(w) l(w) / C(x),
Z=log2 integral 2^(kappa r) + log2 integral 2^(kappa s),
L(x)=log2 C(x).
```

All integrals use normalized circle measure. The source scores reduce to

```text
u=(kappa-1)T+L-Z,
v=(2-kappa)T-2L+Z.                                          (10)
```

Indeed, the original normalizers of `f,g` cancel from the conditional kernel,
and their `kappa` powers cancel from the normalized marginal tests. The same
identities hold for finite-circle averages, with every normalizer replaced.

Compare each circle cell `[w_j,w_j+h)` with its left endpoint. Since
`|l(w)-l(w_j)|<=Ah`, the positive inner integrand changes by a factor between
`2^(-Ah)` and `2^(Ah)`. Thus

```text
|L_P-L|<=Ah, |Z_P-Z|<=kappa Ah.
```

The step density obtained by normalizing the endpoint values differs from
the continuous conditional density by a factor between `2^(-2Ah)` and
`2^(2Ah)`: one factor comes from the point value and one from normalization.
Using `|l|<=B`, the difference of the two expectations is at most

```text
|T_P-T| <= Ah+B(2^(2Ah)-1) <= Ah(1+2B),                      (11)
```

provided `Ah<=1/2`. The last inequality follows from convexity:
`2^t<=1+t` for `0<=t<=1`. This compares two normalized positive measures;
no unproved absolute lower bound for `C` is used in (11).

Since `2pi<8`, it is enough that `P>=16A`. Equations (10)-(11) give

```text
|u_P-u| <= (8A/P)[|kappa-1|(1+2B)+1+kappa],
|v_P-v| <= (8A/P)[|2-kappa|(1+2B)+2+kappa].                  (12)
```

For `R1=lambda1(Delta-u)` and `R2=lambda2(av-u)`, let `C_inner` be eight
times `A` times the larger of the two corresponding brackets in (12),
including their slope and `a` factors. At every outer grid point,

```text
|Rr_P-Rr| <= C_inner/P.
```

## 2. Continuous outer quadrature

Let `R_r` now denote the continuous exponents. The same raw formulas imply

```text
|u|<=B(|kappa-1|+1+kappa),
|v|<=B(|2-kappa|+2+kappa).
```

Choose an integer `K` bounding both `|R1|` and `|R2|` using these inequalities.
For their derivatives, only `r(x+w)` depends on `x`. Differentiating a
normalized expectation gives

```text
T'=E_K[r']+ln(2) Cov_K(l,r'),
|T'|<=||r'||_infinity(1+2B), |L'|<=||r'||_infinity.
```

Here `|l-E_K l|<=2B` and `ln 2<1`. Substitute these bounds in (10) and
then in the two `R_r` to obtain a common derivative upper bound `D`.

For `H_r(x)=(2+2cos x)2^Rr(x)`,

```text
||H_r'||_infinity <= 2^K(2+4D),
I_r=integral H_r >= 2^(1-K).
```

The normalized grid error of a periodic function with derivative bound `J`
is at most `hJ/2<4J/P`, by integrating the endpoint-to-point Lipschitz
error within every cell. Therefore its relative error is at most

```text
C_outer/P, C_outer=2^(2K+1)(2+4D).
```

Using `log2(1+t)<=2t` for `t>=0` (because `ln 2>1/2`), and then the finite
inner error from (12), proves

```text
log2 I_r,P <= log2 I_r + (C_inner+2C_outer)/P.               (13)
```

The zero of `2+2cos x` at `pi` causes no difficulty: the proof uses an
absolute derivative bound and the positive integral lower bound, not a
logarithmic derivative of that vanishing factor.

## 3. Exact retained cutoff

Put `s=eta-chi-48epsilon`. Choose the smallest dyadic modulus threshold not
below

```text
max(2, ceil(16A), ceil[2(C_inner+2C_outer)/s]).
```

The retained exact calculation gives `K=162` and threshold `2^351`.
Equation (13) then bounds each finite log moment by `chi+s/2=0.51841070`.
Finally `eta*681-2>=351`, converting the modulus threshold to the stated
input-size component bound. The calculator's tests also reproduce this
cutoff with independent, slightly looser rational envelopes.

The bounded calculator accepts at most 32 harmonics per polynomial, total
amplitude at most 32, derivative sum at most 512, `0<kappa<=4`,
`0<lambda1,lambda2<=16`, `eta<=1`, and `0<a<=1`. These are computation limits for this
receipt, not impossibility bounds on filters outside that domain.

Replay from the repository root:

```sh
python3 -B research/checksum-variance/finite_circle.py --author conorholds --output /tmp/checksum-finite.json
```

The source is OpenAI's pinned filtering construction and its finite-circle
definitions. This explicit transfer bound is part of the checksum audit;
independent review of the analytic argument remains pending.
