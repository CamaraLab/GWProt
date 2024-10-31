Stress Comparison
========================


.. module:: GWProt.stress_comparison

.. autoclass:: GWProt.stress_comparison.Stress_Comparison

This module stores the transport plans between all pairs of proteins, thus is very memory intensive.
Setting ``RAM = False`` saves these to files, but significantly slows down computations due to file I/O. 


**Calculating Stresses and Distances**
------------------------------------------

These three methods run all pairwise computations. For each pair of proteins they run GW/FGW, store the distance, transport plan, and associated stresses. One of these must be run before any of the latter methods are called. 

As these involve a large number of computations, they have the potential to be time consuming on large datasets.

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.GW_compute_stresses

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.FGW_compute_stresses_data_lists

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.FGW_compute_stresses_dict


These two method must be run after computing the stresses:

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.get_GW_dmat



.. autofunction:: GWProt.stress_comparison.normalize_stress_dict

This method combines the raw stresses which have multiple stress array for each protein to a single stress array for each protein.

For each protein it return the ``normalized_stress`` calculated below, where the rows of ``mat`` are its different stresses in the raw stress dict:

.. code-block:: python

	a, b, c, d, e = code
	normalized_stress = np.sum(mat**a * (np.sum(mat**b, axis=1) ** c)[:, np.newaxis] * mat.shape[1] ** d, axis=0) * mat.shape[0] ** e


**Transferring Stresses**
----------------------------------



.. autofunction:: GWProt.stress_comparison.Stress_Comparison.raw_transferred_stresses

This is in the same format as ``self.raw_stress_dict`` so must be normalized before further use. 



**Analysis Helper Methods**
--------------------------------

We also havehelper methods not part of the ``Stress_Comparison`` class but can be useful for analyzing stresses.



.. autofunction:: GWProt.stress_comparison.get_percentile_of_dict


.. autofunction:: GWProt.stress_comparison.get_AP_scores







