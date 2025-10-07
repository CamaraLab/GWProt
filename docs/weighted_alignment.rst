Weighted Alignment
========================

Weighted alignment methods in GWProt allow for optimal superposition of two point clouds (e.g., protein structures) using a transport plan as weights. This is particularly useful for comparing structures after GW or FGW alignment, where the transport plan encodes the correspondence between residues.

**Weighted RMSD Alignment**
--------------------------

.. autofunction:: GWProt.weighted_alignment.weighted_RMSD

This function finds a special orthogonal matrix :math:`S` that minimizes the weighted sum of squared distances:

    .. math::
        \sum_{i,j} |(x_i - x') - S(y_j - y')|^2 \cdot T_{i,j}

where :math:`x'` is the weighted mean of the :math:`x_i` and :math:`y'` is the weighted mean of the :math:`y_j`, and :math:`T_{i,j}` is the transport plan (e.g., from GW alignment).

Note: In general, there may not be a unique solution matrix :math:`S` that minimizes this expression.

After alignment, ``(Y - y_mean) @ rot + x_mean`` and ``X`` should be superimposed. If ``n == m`` and ``T`` is the identity matrix, this reduces to the classic Kabsch algorithm for minimizing RMSD between ``X`` and ``Y``.

This approach is especially useful for visualizing and quantifying structural similarity after optimal transport-based alignments, and can be used in conjunction with LGD analysis.
