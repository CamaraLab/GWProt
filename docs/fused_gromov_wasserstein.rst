Fused Gromov Wasserstein
========================

Fused Gromov Wasserstein is a variant of Gromov Wasserstein that allow biochemical data to be incorporated.

In addition to its coordinates each residue now has an added datapoint which can be used to compare it to residues in other proteins. 
For instance this could be the difference in isoelectric points, or a substitution score based on the BLOSUM62 matrix. 
Fused GW now adds a penalty term for aligning dissimilar residues. 

Formally, we write :math: `\delta(x,y)` to denote the difference in the data associated to residues :math:`x` and :math:`y`. Then the fused Gromov Wasserstein distance between proteins  :math:`X` and :math:`Y` is defined as


.. math::  FGW(X,Y) = \frac{1}{2} \sqrt{ \sum_{i,j,k,l} (\alpha |d_X(x_i,x_j) - d_Y(y_k,y_l)|^2  + (1 - \alpha) \delta(x_i,y_j))T_{i,k}T_{j,l}}.


where :math:`\alpha` in the interval [0,1] determines the weight of the penalty relative to the usual geometric cost. We found 0.05 to work well for isoelectic points though good values will vary depending on the type of biochemical data used.

[https://arxiv.org/pdf/1811.02834]
