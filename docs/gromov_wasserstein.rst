.. -*- coding: utf-8 -*-

Gromov-Wasserstein
==================

Gromov-Wasserstein (GW) metric is a method  of quantifying the difference between shapes or point clouds. In GWProt we treat proteins as shapes in 3-dimensional space defined by their alpha-Carbons. 

Intuitively GW compares proteins by structurally aligning their residues in a way that minimizes distortion. 

As a starting point we can form a one to one pairing :math:`f` between the residues in protein :math:`X` with those in protein :math:`Y` that minimizes the largest distortion - the difference between a distance within :math:`X` and the corresponding distance in :math:`Y`. This defines the *Gromov-Hausdorff distance* between :math:`X` and :math:`Y`:

.. math::  d_{GH}(X,Y) = \min_{f :X\cong Y} \max_{i,j \in X} \lvert d_X(i,j) - d_Y(f(i),f(j)) \rvert .

However this in not computable in practice as the number of possibilities for :math:`f` is :math:`|X|!`. 


Thus we adjust our approach by turning it into a continuous problem which can be efficiently approximated[^1]. 
[^1] In almost all cases the calculated approximation is equivalent to the precise GW metric for practical purposes. The key exception is that the calculated approximations do not satisfy the mathematical formulae that the actual metric does. For instance the triangle inequality may not be satisfied.
To do this we give each protein a mass of 1 and distribute it evenly among its residues, we call this a *distribution*. 
Aligning two proteins of lengths *n* and *m* then amounts to transferring the mass of one protein the the other. We call this assignment a *transport plan* and represent it as a *n x m* matrix where each column sums to *1/m* and each row sums to *1/n*. The *ij*th entry is thus the amount of mass tranported from the *i*th residue of one protein to the *j*th residue of the other. 
Finding the best alignment is now finding the best transport plan.


We define the Gromov-Wasserstein distance based on the sum of all distortions weighted by the optimal tranport plan:

.. math::  GW(X,Y) = \min_T \frac{1}{2} \sqrt{ \sum_{i,j,k,l} |d_X(x_i,x_j) - d_Y(y_k,y_l)|^2  T_{i,k}T_{j,l}}.


Along with calculating the GW distance, we also have the optimal transport plan.

(The squaring and squareroot is not mathematically necessary, but is needed for efficient computation.)






