.. -*- coding: utf-8 -*-

Gromov-Wasserstein
==================

Gromov-Wasserstein is a method  of quantifying the difference between shapes or point clouds. In GWProt we treat proteins as shapes in 3-dimensional space defined by their alpha-Carbons. 


Gromov-Wasserstein defines a metric, implying that the GW-distance between two shapes is zero if and only if the shapes are identical and that for three shapes :math:`X` , :math:`Y`, and :math:`Z`, 

.. math::  GW(X,Y) + GW(Y,Z) \geq GW(X,Z)

However as the computer does not output the exact GW distance but an approximation of it, this inequality is not guaranteed to be satisfied. Instead we have :math:`GW_{computed}(X,Y) + GW_{computed}(Y,Z) + \epsilon \geq GW_{computed}(X,Z)` for a small value of :math:`epsilon`.  




