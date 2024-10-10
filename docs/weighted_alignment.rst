Weighted Alignment
========================


.. autofunction:: GWProt.weighted_alignment.weighted_RMSD


Explicitly it finds a special orthogonal matrix *S* which minimizes

    .. math:: \sum_{i,j} |(x_i - x')- S(y_j - y')|^2 * T_{i,j} 

where *x'* is the weighted mean of the *x_i* and *y'* is the weighted mean of the *y_j*.

Note - in general there may not be a unique solution matrix *S* which minimizes this. 




Thus ``(Y - y_mean) @ rot + x_mean`` and ``X`` should be superimposed.
Note that if ``n == m`` and ``T`` is the identity matrix, this is just the usual Kabsch algorithm for minimizing RMSD between ``X`` and ``Y``.
