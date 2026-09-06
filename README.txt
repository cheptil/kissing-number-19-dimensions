Verification material for
"Optimality of the added-vector code in the 19-dimensional kissing construction"
Alexey Mikhailovich Kolosov

CONTENTS
  verify_golay_1280_optimality.py -- exact verification and optional LP discovery
  certificate.json              -- explicit mathematical input and output data
  verification.txt              -- output of the default exact verification
  discovery.txt                 -- output of the optional numerical run
  requirements-discovery.txt    -- package versions used for the numerical run

EXACT VERIFICATION
  python verify_golay_1280_optimality.py

Python 3.10 or later is required. Only the standard library is needed.
All validation conditions remain active under python -O.
The verifier reconstructs the code from the twelve generator supports;
checks dimensions and the first-twelve-column rank certificate;
checks the five explicit generator identities, the induced sixteen-vertex
forbidden-distance graph, and its independence number by exact subset DP;
checks the exact Fourier dual certificate; and reconstructs both known
matching lower-bound codes.

The 320-word code B has exact minimum distance 6. The 1280-word code A
has exact minimum distance 5. Both are admissible for the requirement
of minimum distance at least 5.

OPTIONAL NUMERICAL DISCOVERY
  python -m pip install -r requirements-discovery.txt
  python verify_golay_1280_optimality.py --discover

The recorded run used Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0.
The numerical LP is a discovery aid, not a correctness certificate.
An LP with multiple optima can return a different dual support on another
solver/version. The fixed certificate is always checked separately.

INPUT SOURCES
  B. S. Ho, A new lower bound for the kissing number in 19 dimensions,
  arXiv:2603.10425v2, Section 2 and Proposition 9:
  https://arxiv.org/html/2603.10425v2

  V. Gonzalez, Maximality of the added-vector codes in the Cohn-Li
  kissing constructions, version 3, Section 7, Conjecture 16:
  https://f-keys.com/papers/maximality-added-vector-codes-cohn-li/
  Version consulted on 6 September 2026.

The input generators and the matching lower-bound construction are prior
work of Ho. The five-word induced-subgraph certificate supplies the
matching upper bound. The proof itself is given in the manuscript.
No upper bound on the unrestricted kissing number k(19) is asserted.
