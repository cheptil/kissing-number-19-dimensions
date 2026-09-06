#!/usr/bin/env python3
"""Exact certificate for a restricted 19-dimensional kissing-number problem.

Proves, for the explicit binary code D of Ho (2026),
    max{|A| : A is a subset of D, minimum Hamming distance >= 5} = 1280.
Also proves alpha(Cay(K,S)) = 320, where S comprises the weight-3/4
words of D and K = span(S).

This does NOT prove that the unrestricted kissing number in dimension 19
is 11948. Nor does this script establish priority of the argument.

Source of the input generators and the known lower-bound construction:
    https://arxiv.org/html/2603.10425v2
The five-word upper-bound certificate was found in this conversation
using the dual of a translation-invariant Lovasz-theta linear program.

Usage (Python 3.10+):
    python verify_golay_1280_optimality.py
    python verify_golay_1280_optimality.py --discover

The default verification uses only the standard library and exact integers.
The optional discovery step requires numpy and scipy, uses floating point,
and is NOT relied upon for the proof.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from collections import Counter
from typing import Iterable, Sequence


def require(condition: bool, message: str) -> None:
    """Do not use assert: verification must remain active under python -O."""
    if not condition:
        raise ValueError(message)


def mask(support: Iterable[int]) -> int:
    support = tuple(support)
    require(len(support) == len(set(support)), "Repeated coordinate")
    require(all(1 <= i <= 19 for i in support), "Coordinate outside 1..19")
    return sum(1 << (i - 1) for i in support)


def support(word: int) -> tuple[int, ...]:
    return tuple(i + 1 for i in range(19) if word & (1 << i))


def span(generators: Iterable[int]) -> set[int]:
    result = {0}
    for g in generators:
        result |= {x ^ g for x in result}
    return result


def ordered_span(basis: Sequence[int]) -> list[int]:
    result = [0]
    for g in basis:
        result += [x ^ g for x in result]
    require(len(set(result)) == len(result), "Dependent basis")
    return result


# Section 2 of Ho's version 2: six m_i, four s_i, and two r_i.
GENERATOR_SUPPORTS = (
    (1, 8, 9, 12, 16, 17, 18, 19),
    (2, 10, 11, 14, 15, 17, 18),
    (3, 7, 9, 13, 15, 16, 17, 19),
    (4, 7, 8, 10, 12, 15, 16, 19),
    (5, 10, 12, 13, 15, 16, 17, 18),
    (6, 7, 8, 9, 10, 13, 16, 18),
    (1, 4, 7, 9),
    (1, 5, 6, 18),
    (1, 3, 12, 15),
    (1, 10, 13, 19),
    (1, 3, 5, 6, 7, 13, 14, 15, 18),
    (2, 4, 6, 7, 8, 13, 14, 16, 17, 18),
)

# These five words are all in Table 1 of that paper.
CERTIFICATE_SUPPORTS = (
    (1, 4, 7, 9),
    (4, 8, 10, 12),
    (7, 8, 13, 18),
    (1, 10, 13, 19),
    (9, 12, 18, 19),
)


def independence_number(vertices: Sequence[int], forbidden: set[int]) -> int:
    """Exact subset-DP: alpha(U)=max(alpha(U-v),1+alpha(U-N[v]))."""
    n = len(vertices)
    require(n <= 20, "This verifier is intended for the small certificate")
    adjacency = [0] * n
    for i, j in combinations(range(n), 2):
        if vertices[i] ^ vertices[j] in forbidden:
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
    dp = bytearray(1 << n)
    for subset in range(1, 1 << n):
        first = subset & -subset
        i = first.bit_length() - 1
        rest = subset ^ first
        dp[subset] = max(dp[rest], 1 + dp[rest & ~adjacency[i]])
    return dp[-1]


def discover(k_order: Sequence[int], forbidden: set[int]) -> None:
    """Solve the translation-invariant theta relaxation on K ~= F_2^10.

    A character distribution p_a >= 0 satisfies sum p_a = 1 and
    sum_a p_a (-1)^(a dot s) = 0 for each forbidden difference s.
    Maximizing |K| p_0 is the theta relaxation.
    """
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError as exc:
        raise RuntimeError("--discover requires numpy and scipy") from exc
    n = len(k_order)
    s_labels = [i for i, w in enumerate(k_order) if w in forbidden]
    aeq = np.array([[1] * n] + [
        [1 - 2 * ((a & s).bit_count() % 2) for a in range(n)]
        for s in s_labels
    ], dtype=float)
    c = np.zeros(n)
    c[0] = -n
    result = linprog(c, A_eq=aeq, b_eq=np.r_[1, np.zeros(len(s_labels))],
                     bounds=(0, None), method="highs")
    require(bool(result.success), f"Discovery LP failed: {result.message}")
    selected = [k_order[s] for s, multiplier in
                zip(s_labels, result.eqlin.marginals[1:])
                if abs(float(multiplier)) > 1e-7]
    print(f"Discovery LP theta value (numerical): {-result.fun:g}")
    print("Dual-supported forbidden differences:", [support(x) for x in selected])
    h_found = sorted(span(selected))
    if len(h_found) <= 20:
        a = independence_number(h_found, forbidden)
        print(f"Extracted subgroup: {len(h_found)} vertices, exact alpha={a}")
        print(f"Its exact coset bound on K: {(n // len(h_found)) * a}")
    else:
        print("An alternative LP optimum was returned; the fixed certificate "
              "below is verified independently.")


def verify(use_discovery: bool = False) -> None:
    generators = [mask(s) for s in GENERATOR_SUPPORTS]
    d_order = ordered_span(generators)
    d = set(d_order)
    m = span(generators[:6])
    k_order = ordered_span(generators[:10])
    k = set(k_order)
    require(len(d) == 4096 and len(k) == 1024 and len(m) == 64,
            "Unexpected code dimensions")
    require(len({x & ((1 << 12) - 1) for x in d}) == 4096,
            "The first twelve coordinate columns are not nonsingular")
    require(min(x.bit_count() for x in d if x) == 3, "Unexpected minimum weight")
    forbidden = {x for x in d if 0 < x.bit_count() < 5}
    require(len(forbidden) == 21 and span(forbidden) == k,
            "Unexpected forbidden-difference structure")

    if use_discovery:
        discover(k_order, forbidden)

    h_generators = [mask(s) for s in CERTIFICATE_SUPPORTS]
    m1, m2, m3, m4, m5, m6, s1, s2, s3, s4, r1, r2 = generators
    identities = [s1, m1 ^ m4 ^ m5 ^ m6 ^ s2, m1 ^ m3 ^ s3, s4,
                  m3 ^ m4 ^ m5 ^ m6 ^ s1 ^ s2 ^ s3 ^ s4]
    require(h_generators == identities, "A displayed membership identity is wrong")
    require(all(x in forbidden for x in h_generators), "Invalid certificate word")
    require(all((x & y).bit_count() == 1 for x, y in combinations(h_generators, 2)),
            "Certificate supports do not intersect pairwise in one position")
    total = 0
    for x in h_generators:
        total ^= x
    require(total == 0, "Certificate words do not sum to zero")
    h = span(h_generators)
    require(len(h) == 16 and h <= k, "Invalid certificate subgroup")
    require(Counter(x.bit_count() for x in h) == {0: 1, 4: 5, 6: 10},
            "Incorrect weight distribution in the certificate subgroup")
    require(h & m == {0}, "The claimed direct sum K = M + H is not direct")

    # Exact validation of the Fourier dual certificate in the manuscript.
    k_labels = {word: i for i, word in enumerate(k_order)}
    h_labels = [k_labels[word] for word in h_generators]
    for character in range(len(k_order)):
        q = 384 + 128 * sum(1 - 2 * ((character & label).bit_count() % 2)
                            for label in h_labels)
        require(q in (0, 512, 1024), "Unexpected dual polynomial value")
        require(q >= (1024 if character == 0 else 0),
                "The exact dual certificate is infeasible")
    require(h & forbidden == set(h_generators), "Unexpected induced graph")
    alpha_h = independence_number(sorted(h), forbidden)
    require(alpha_h == 5, "The small graph does not have independence number 5")

    # K and D partition into cosets of H. Every coset has exactly the same
    # induced forbidden-difference graph. Therefore each contains <= alpha_h
    # words of any distance-5 code, irrespective of cross-coset constraints.
    upper_k = (len(k) // len(h)) * alpha_h
    upper_d = (len(d) // len(h)) * alpha_h
    require(upper_k == 320 and upper_d == 1280, "Unexpected upper bounds")

    # Independently rebuild the known matching lower-bound construction.
    s_reps = generators[6:10] + [mask((3, 5, 7, 10))]
    b = {s ^ x for s in s_reps for x in m}
    r1, r2 = generators[10:12]
    a = {x ^ r for x in b for r in (0, r1, r2, r1 ^ r2)}
    require(b <= k and len(b) == 320, "Invalid 320-word construction")
    require(min((x ^ y).bit_count() for x, y in combinations(b, 2)) == 6,
            "Unexpected minimum distance in B (it is 6, not 5)")
    require(a <= d and len(a) == 1280, "Invalid 1280-word construction")
    require(all((x ^ s) not in a for x in a for s in forbidden),
            "Lower-bound construction violates minimum distance 5")
    require(min((x ^ y).bit_count() for x, y in combinations(a, 2)) == 5,
            "Unexpected exact minimum distance")

    print("PASS: all five generator identities and the first-12-column rank certificate.")
    print("PASS: weight distribution of H is 1 + 5 z^4 + 10 z^6; K = M direct-sum H.")
    print("PASS: the Fourier dual polynomial is feasible in exact arithmetic.")
    print("PASS: |D|=4096, |K|=1024, |H|=16; 21 forbidden differences.")
    print("PASS: exact independence number of the 16-vertex certificate graph is 5.")
    print("PASS: 64 H-cosets in K give alpha(K) <= 320.")
    print("PASS: 256 H-cosets in D give alpha(D) <= 1280.")
    print("PASS: explicit matching constructions of sizes 320 and 1280.")
    print("PASS: their exact minimum distances are 6 and 5, respectively.")
    print("PROVED for these explicit codes: alpha(K)=320 and alpha(D)=1280.")
    print("No assertion of an upper bound on the unrestricted kissing number k(19).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discover", action="store_true",
                        help="also run the optional floating-point discovery LP")
    args = parser.parse_args()
    verify(args.discover)
