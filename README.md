# Maximal Robust Forward Invariant Set

Recursive feasibility and stability of robust model predictive control (MPC) usually depends on the notion of robust forward invariant sets, also known as robust postively invariant sets. This matlab script allows the computation of such a step under the following conditions:

- The plant dynamics are stabilizable
- The constraint sets and bounded process noise are defined by compact convex polytopes described by a system of linear inequalities.

It combines algorithms presented in Gilbert, E. G., & Tan, K. (1991). "Linear systems with state and control constraints: The theory and 
application of maximal output admissible sets", IEEE Transactions on Automatic Control, 36, 1008–1020, and Borrelli, F.,Bemporad, A. & Morari M., "Predictive Control for Linear and Hybrid Systems", Cambridge University Press, 2017, ISBN 1107016886, 9781107016880.
