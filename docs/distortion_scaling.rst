Distortion Scaling
=====================


When comparing the shapes of molecules, we are generally more interested in atoms and substructures which are (geometrically) nearby one another as those are more likely to interact chemically. This is the key idea of the template-modelling (TM) score and the alignment software TM-align. To implement this idea in GWProt, we replace the intra-protein distance matrice of our proteins with related ones, causing the GW computation to give greater weight to nearby residues. We found that this generally improves accuracy with no impact on runtime.






We choose a *scaling function* :math:`f` such that :math:`f(0) = 0`, :math:`f` is strictly monotonic increasing, and :math:`f` is concave down. We found the square root function works well. Then for all proteins we apply :math:`f` to all entries in their intra-protein distance matrices. Finally we run GW using the new intra-protein distance matrices.

We can use distortion scaling with FGW the same way as with GW. The only difference is we must use a larger value of ``alpha``; we found 0.5 to work well with isoelectric points.


Warning - Do not use GW or FGW to compare a scaled protein and and unscaled one. Similarly do not compare two proteins scaled with different scaling functions. The resulting numbers will be meaningless.

Warning - The distances given by GW or FGW after scaling will be smaller than those without scaling. We cannot directly compare distances between scaled proteins with distances between unscaled ones; similarly we cannot compare distances from two different scaling functions.



