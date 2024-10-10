.. -*- coding: utf-8 -*-

Gromov-Wasserstein
==================

Gromov-Wasserstein (GW) metric is a method  of quantifying the difference between shapes or point clouds. In GWProt we treat proteins as shapes in 3-dimensional space defined by their alpha-Carbons. 

Intuitively GW compares proteins by structurally aligning their residues in a way that minimizes distortion. 

As a starting point we can form a one to one pairing :math:`f` between the residues in protein :math:`X` with those in protein :math:`Y` that minimizes the largest distortion - the difference between a distance within :math:`X` and the corresponding distance in :math:`Y`. This defines the *Gromov-Hausdorff distance* between :math:`X` and :math:`Y`:

.. math::  d_{GH}(X,Y) = \min_{f :X\cong Y} \max_{i,j \in X} \lvert d_X(i,j) - d_Y(f(i),f(j)) \rvert .

However this in not computable in practice as the number of possibilities for :math:`f` is :math:`|X|!`. 


Thus we adjust our approach by turning it into a continuous problem which can be efficiently approximated. 




GWProt relies on intra-protein distances as opposed the the raw coordinates, thus is not affected by rotations or translations.



.. math::  FGW(X,Y) = \frac{1}{2} \sqrt{ \sum_{i,j,k,l} |d_X(x_i,x_j) - d_Y(y_k,y_l)|^2  T_{i,k}T_{j,l}}.




Gromov-Wasserstein defines a metric, implying that the GW-distance between two shapes is zero if and only if the shapes are identical and that for three shapes :math:`X` , :math:`Y`, and :math:`Z`, 

.. math::  GW(X,Y) + GW(Y,Z) \geq GW(X,Z)

However as the computer does not output the exact GW distance but an approximation of it, this inequality is not guaranteed to be satisfied. Instead we have :math:`GW_{computed}(X,Y) + GW_{computed}(Y,Z) + \epsilon \geq GW_{computed}(X,Z)` for a small value of :math:`epsilon`.  




