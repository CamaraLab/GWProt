Stress
=========

When comparing two proteins with GW, the **stress** of a residue is defined as its contribution to the GW distance.

Formally, when comparing proteins :math:`X` and :math:`Y`, the stress of residue :math:`x_i` is

.. math:: stress(x_i) = \sum_{j,k,l} |d_X(x_i,x_j) - d_Y(y_k,y_l)|^2 T_{i,k}T_{j,l} .



Note that :math:`\frac{1}{2}\sqrt{\sum_i stress(x_i)} = GW(X,Y)` and similarly for :math:`Y`. 

Intuitively this quantifies how well the individual residue :math:`x_i` is structurally conserved compared to :math:`Y.`
A higher stress level indicates lower conservation.
For instance this can be used to identify structurally conserved regions in evolutionarily related proteins, 
where lower stress regions could indicate active sites.
In another example, higher stress regions can indicate flexible switch regions when comparing multiple files of the same or similar proteins. 

In practice, it is better to average the stresses associated to GW computations across multiple proteins 
than using the stresses associated to a single comparison. It is not recommended to use stresses of downsampled proteins.

.. warning:: 
	This notion of stress is unrelated to chemical forces acting within a protein or the frustration of a protein conformation.

We can also define the **fused stress** when using fused GW. Similarly it is a single residue's contribution to the FGW distance.
Formally:

.. math:: stress(x_i) = \sum_{j,k,l} (\alpha \cdot |d_X(x_i,x_j) - d_Y(y_k,y_l)|^2  + (1 - \alpha) \delta(x_i,y_j))T_{i,k}T_{j,l} .