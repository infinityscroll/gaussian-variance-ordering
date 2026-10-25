#!/usr/bin/env python3
"""Exact finite algebra checks for the three-dimensional VI counterexample.

Requires only Python 3's standard library. No floating-point arithmetic or
numerical optimization is used. This is NOT a formal proof checker. It does
not verify the Gaussian integration, strict-convexity/global-optimality
argument, intermediate value theorem, arbitrary-alpha construction, or
general entropy theorem. Those are analytic arguments in the accompanying
paper. It checks the explicit rational inputs and algebra supporting the
finite example at alpha=9 and the root bracket at alpha=10000.
"""

from fractions import Fraction as Q
from itertools import combinations
import json


checks = []


def check(condition, name):
    if not condition:
        raise RuntimeError("FAILED: " + name)
    checks.append(name)


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def inverse(a):
    n = len(a)
    z = [list(row) + ident for row, ident in zip(a, eye(n))]
    for c in range(n):
        p = next((r for r in range(c, n) if z[r][c]), None)
        if p is None:
            raise RuntimeError("singular exact matrix")
        z[c], z[p] = z[p], z[c]
        pivot = z[c][c]
        z[c] = [x / pivot for x in z[c]]
        for r in range(n):
            if r != c:
                factor = z[r][c]
                z[r] = [x - factor * y for x, y in zip(z[r], z[c])]
    return [row[n:] for row in z]


def det(a):
    if len(a) == 1:
        return a[0][0]
    return sum((-1) ** j * a[0][j] *
               det([[row[k] for k in range(len(a)) if k != j]
                    for row in a[1:]]) for j in range(len(a)))


def spd_minors(a, label):
    n = len(a)
    check(all(a[i][j] == a[j][i] for i in range(n) for j in range(n)),
          label + ": symmetry")
    result = {}
    for count in range(1, n + 1):
        for ix in combinations(range(n), count):
            value = det([[a[i][j] for j in ix] for i in ix])
            key = ",".join(str(i + 1) for i in ix)
            check(value > 0, label + ": principal minor " + key)
            result[key] = str(value)
    return result


I = eye(3)
Sigma = [[Q(x, 12648) for x in row] for row in
         [[5673, 4650, 4650], [4650, 9854, -1825],
          [4650, -1825, 9854]]]
K = [[Q(1), Q(3, 4), Q(3, 4)],
     [Q(3, 4), Q(1), Q(19, 50)],
     [Q(3, 4), Q(19, 50), Q(1)]]
K_inverse_printed = [[Q(x, 1581) for x in row] for row in
                     [[8556, -4650, -4650], [-4650, 4375, 1825],
                      [-4650, 1825, 4375]]]

sigma_minors = spd_minors(Sigma, "Sigma")
k_shift_minors = spd_minors(add(K, scale(Q(-1, 9), I)), "K-I/9")
check(det(K) == Q(1581, 10000), "printed determinant of K")
check(inverse(K) == K_inverse_printed, "printed inverse of K")
check(mul(K, K_inverse_printed) == I, "K inverse multiplication")
check(Sigma == scale(Q(1, 8), add(scale(Q(9), I),
                                scale(Q(-1), K_inverse_printed))),
      "Sigma=(9I-K^{-1})/8")

H9 = add(scale(Q(9), I), scale(Q(-8), Sigma))
check(H9 == K_inverse_printed, "H_9(I)=K^{-1}")
check(inverse(H9) == K, "H_9(I) inverse")
check(all(K[i][i] == 1 for i in range(3)), "D(9)=I stationarity equations")
phi9 = add(scale(Q(9), inverse(Sigma)), scale(Q(-8), I))
phi9_minors = spd_minors(phi9, "Phi_9(I)")

hadamard = [[x * x for x in row] for row in K]
M = add(scale(Q(9), hadamard), scale(Q(-1), I))
rhs = [[(sum(row) - 1) / 8] for row in hadamard]
v = [[Q(x, 14808952)] for x in [264375, -3209, -3209]]
derivative_minors = spd_minors(M, "9(K o K)-I")
check(mul(M, v) == rhs, "printed derivative satisfies linear system")
check(mul(inverse(M), rhs) == v, "unique exact derivative solution")
check(v[0][0] > 0 and v[1][0] == v[2][0] < 0,
      "one positive and two negative derivative entries")
check(sum(row[0] for row in v) > 0, "logdet derivative positive at alpha=9")

A, B, S, T = (Q(x, 12648) for x in [5673, 4650, 8029, 11679])
y0 = Q(19, 20)
q0, d0 = y0 - S, y0 - T
finite_bounds = [
    (T > S > 0, "T>S>0"),
    (Q(3, 10) < q0 < Q(1, 3), "3/10<q0<1/3"),
    (d0 > Q(1, 40), "d0>1/40"),
    (0 < A < Q(1, 2), "0<A<1/2"),
    (0 < B * B < Q(1, 7), "0<B^2<1/7"),
    (A > Q(56, 125), "A>56/125"),
    (S < Q(127, 200), "S<127/200"),
    (B * B < Q(17, 125), "B^2<17/125"),
    (y0 - 2 * S < 0, "y0-2S<0"),
    (q0 > Q(63, 200), "q0>63/200"),
    (y0 - 2 * S > Q(-8, 25), "y0-2S>-8/25"),
]
for condition, name in finite_bounds:
    check(condition, name)


def f(y):
    return A * (y - S) ** 2 + B * B * (y - 2 * S)


f_lower = Q(56, 125) * Q(63, 200) ** 2 - Q(17, 125) * Q(8, 25)
den0 = q0 * (A * q0 + 2 * B * B)
check(f_lower == Q(583, 625000), "printed f lower bound")
check(f(y0) > f_lower > Q(9, 10000), "strict f(y0) lower bounds")
check(0 < den0 < Q(19, 126) < Q(1, 6), "strict denominator upper bounds")
check(f(y0) / den0 > Q(27, 5000), "positive uniform sign margin")
check(f(T) < 0 < f(y0), "limiting determinant root bracket signs")


def scalar(a, y):
    r, k = a - 1, a / (a - 1)
    q, d = k * y - S, k * y - T
    x = A + 2 * B * B / q
    e = y / (2 * r) * (1 / q + 1 / d) + y * B * B / (q * (A * q + 2 * B * B)) - 1
    return r, k, q, d, x, e


check(scalar(Q(9), Q(1))[-1] == 0, "scalar equation vanishes at (9,1)")
a = Q(10000)
k = a / (a - 1)
ylo = T / k + Q(1, 1000000)
trial_values = {}
for label, y in [("lower", ylo), ("upper", y0)]:
    r, k, q, d, x, e = scalar(a, y)
    check(q > 0 and d > 0 and x > 0, label + ": positive scalar domain")
    check((k * x - A) * q - 2 * B * B == x * q / r,
          label + ": exact block determinant identity")
    D = [[x, Q(0), Q(0)], [Q(0), y, Q(0)], [Q(0), Q(0), y]]
    H = add(scale(a, D), scale(-r, Sigma))
    spd_minors(H, label + ": H_a(D)")
    Hinv = inverse(H)
    check(Hinv[0][0] == 1 / x, label + ": exact center stationarity")
    check(y * Hinv[1][1] - 1 == e == y * Hinv[2][2] - 1,
          label + ": exact leaf residual identity")
    trial_values[label] = {"y": str(y), "E": str(e)}

E_lower = scalar(a, ylo)[-1]
E_upper = scalar(a, y0)[-1]
check(T / k < ylo < y0, "nonempty rational root bracket")
check(E_lower > 0 > E_upper, "opposite exact scalar signs at alpha=10000")
check(E_upper < Q(65, 3 * 9999) - Q(27, 5000), "printed scalar upper bound")
check(Q(65, 29997) < Q(1, 400), "65/29997<1/400")
check(Q(1, 400) - Q(27, 5000) < 0, "uniform final margin negative")

print(json.dumps({
    "status": "PASS_EXACT_FINITE_ALGEBRA",
    "check_count": len(checks),
    "arithmetic": "Python standard-library fractions.Fraction; no floating point",
    "scope": "Finite rational matrix identities, positive-definiteness minors, derivative system and scalar bracket/sign checks only; not formal or complete analytic proof verification.",
    "sigma_principal_minors": sigma_minors,
    "k_minus_identity_over_9_principal_minors": k_shift_minors,
    "phi9_principal_minors": phi9_minors,
    "derivative_system_principal_minors": derivative_minors,
    "derivative": [str(row[0]) for row in v],
    "alpha10000_bracket": trial_values,
    "checks": checks,
}, indent=2))
