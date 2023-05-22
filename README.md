# dodo
## Simultaneous diagonalization of nearly commuting Hermitian matrices: do-one-then-do-the-other

Commuting Hermitian matrices $A$ and $B$ may be simultaneously diagonalized by a common unitary matrix $U$ of eigenvectors:

$U^H A U = \Lambda$, $U^H B U = M$.

`dodo.m` computes $U$, $\Lambda$, and $M$. The routines under `test/` run the experiments described in the article

"Simultaneous diagonalization of nearly commuting Hermitian matrices: do-one-then-do-the-other"

published in _IMA Journal of Numerical Analysis_.
