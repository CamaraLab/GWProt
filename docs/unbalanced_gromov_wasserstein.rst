Unbalanced Gromov-Wasserstein Correspondences
=========================================

Overview
--------
Unbalanced Gromov-Wasserstein (UGW) [1]_ is a variant of the Gromov-Wasserstein distance that permits some of the mass to be discarded instead of transported between two proteins. This can be more robust to proteins of different lengths and indel mutations.


Mathematical Formulation
------------------------
 The *unbalanced Gromov-Wasserstein distance* between proteins :math:`X` and :math:`Y` is defined as:

.. math::
   UGW(X, Y) &= \min_{T \geq 0} \frac{1}{2} \left( \sum_{i, j, k, l} \left[   |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2   \cdot T_{i, k} T_{j, l} \right] 
    &+ \rho \cdot KL^\otimes(\pi_1(T) | \mu_X) + \rho \cdot KL^\otimes(\pi_2(T) | \mu_Y)  \right)^{1/2}

   
Where :math:`KL` denotes Kullback-Leibler divergence, :math:`\mu_1, \mu_2` are the probability distributions on :math:`X` and :math:`Y`, and :math:`\pi_1(T), \pi_2(T)` are the projection of :math:`T` onton :math:`X` and :math:`Y`. Note that :math:`T` is non-negative, but does not have the strict marginal constraints as in GW and FGW.
Optimal values of :math:`\rho` may vary depending on the proteins, and higher :math:`\rho` means more mass is preserved. 

We similarly have *fused unbalanced Gromov-Wasserstein distance* defined as: 

.. math::
   UGW(X, Y) &= \min_{T \geq 0} } \frac{1}{2} \left( \sum_{i, j, k, l} \left[   |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2    \cdot T_{i, k} T_{j, l} \right]  
    &+ \rho \cdot KL^\otimes(\pi_1(T) | \mu_X) + \rho \cdot KL^\otimes(\pi_2(T) | \mu_Y)   
    &+ (1 - \alpha) \cdot \sum_{i,j} \delta(x_i, y_k) \cdot T_{i, k}\right)^{1/2}

Where :math:`\delta(x, y)` and :math:`\alpha` are as in FGW.

References
----------
.. [1] Séjourné, T.,  Vialard, F., and Peyré, G. (2021) The Unbalanced Gromov Wasserstein Distance: Conic Formulation and Relaxation. Neural Information Processing Systems, 35.