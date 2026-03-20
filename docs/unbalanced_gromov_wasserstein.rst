Unbalanced Gromov-Wasserstein Correspondences
=================================================

Overview
------------
Unbalanced Gromov-Wasserstein (UGW) [1]_ is a variant of the Gromov-Wasserstein distance that permits some of the mass to be discarded instead of transported between two proteins. This can be more robust to proteins of different lengths and indel mutations. However it does not define a distance metric.


Mathematical Formulation
----------------------------
 The *unbalanced Gromov-Wasserstein distance* between proteins :math:`X` and :math:`Y` is defined as:

.. math::
    UGW(X, Y) 
    = \min_{T \geq 0} \frac{1}{2} \Big( \sum_{i, j, k, l} [   |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2   \cdot T_{i, k} T_{j, l} ] 

    + \rho \cdot KL^\otimes(\pi_1(T) | \mu_X) + \rho \cdot KL^\otimes(\pi_2(T) | \mu_Y)  \Big)^{1/2}

   
Where :math:`KL` denotes Kullback-Leibler divergence,  :math:`KL^\otimes(\mu|\nu)` denotes :math:`KL(\mu \otimes \mu|\nu \otimes \nu)` , :math:`\mu_1, \mu_2` are the probability distributions on :math:`X` and :math:`Y`, and :math:`\pi_1(T), \pi_2(T)` are the projections of :math:`T` onto :math:`X` and :math:`Y`, in this case uniform distributions.
Note that :math:`T` is non-negative, but does not have the strict marginal constraints as in GW and FGW. Instead the marginal constraints are weakly enforced via Kullback-Leibler divergence.
Optimal values of :math:`\rho` may vary depending on the proteins, and higher :math:`\rho` means more mass is preserved. 

We similarly have *fused unbalanced Gromov-Wasserstein distance* defined as: 

.. math::
    UGW(X, Y) 
    = \min_{T \geq 0}  \frac{1}{2} \Bigg( \sum_{i, j, k, l} [   |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2    \cdot T_{i, k} T_{j, l}]  \\
    + \rho \cdot KL^\otimes(\pi_1(T) | \mu_X) + \rho \cdot KL^\otimes(\pi_2(T) | \mu_Y)   \\
      (1 - \alpha) \cdot \sum_{i,j} \delta(x_i, y_k) \cdot T_{i, k} \Bigg)^{1/2}

Where :math:`\delta(x, y)` and :math:`\alpha` are as in FGW.





References
-------------
.. [1] Séjourné, T.,  Vialard, F., and Peyré, G. (2021) The Unbalanced Gromov Wasserstein Distance: Conic Formulation and Relaxation. Neural Information Processing Systems, 35.