Local Geometric Distortion
==========================

Overview
--------

When comparing two proteins with GW, the **local geometric distortion (LGD)** of a residue 
quantifies its contribution to the overall GW distance. This provides a residue-level measure 
of structural conservation.

Definition
----------

Formally, when comparing proteins :math:`X` and :math:`Y`, the local geometric distortion of 
residue :math:`x_i` is defined as:

.. math::
   LGD(x_i) = \sum_{j, k, l} |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2 T_{i, k} T_{j, l}

where :math:`T` is the optimal transport plan between residues of :math:`X` and :math:`Y`.

Relationship to GW Distance
---------------------------

The sum of local geometric distortions relates directly to the GW distance:

.. math::
   \frac{1}{2} \sqrt{\sum_i LGD(x_i)} = GW(X, Y)

and similarly for residues in :math:`Y`.

Interpretation and Applications
------------------------------

- **Low LGD:** Indicates that a residue is structurally well-conserved relative to the other protein.
- **High LGD:** Indicates lower conservation, which may correspond to flexible or functionally 
important regions (e.g., catalytic sites or switch regions).

LGD can be used to:

- Identify structurally conserved regions in evolutionarily related proteins.
- Detect flexible or variable regions by comparing multiple conformations of the same or similar proteins.

.. note::

   For robust results, it is recommended to average LGD values across multiple protein comparisons rather than relying on a single pairwise comparison. Avoid using LGD values from downsampled proteins.

Fused Local Geometric Distortion
-------------------------------

When using fused GW, the **fused local geometric distortion** extends the concept to 
include biochemical data. It is defined as:

.. math::
   LGD(x_i) = \sum_{j, k, l} \left[ \alpha |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2 + (1 - \alpha) \delta(x_i, y_j) \right] T_{i, k} T_{j, l}

where :math:`\delta(x_i, y_j)` measures the biochemical difference between residues, 
and :math:`\alpha` controls the balance between geometric and biochemical contributions.