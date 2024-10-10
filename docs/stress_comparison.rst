Stress Comparison
========================


.. module:: GWProt.stress_comparison

.. autoclass:: GWProt.stress_comparison.Stress_Comparison

This module stores the transport plans between all pairs of proteins, thus is very memory intensive.
Setting ``RAM = False`` saves these to files, but significantly slows down computations due to file I/O. 


These two methods run all pairwise computations. For each pair of proteins they run GW/FGW, store the distance, transport plan, and associated stresses. One of these must be run before any of the latter methods are called. 

As these involve a large number of computations, they have the potential to be time consuming on large datasets.

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.GW_compute_stresses

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.FGW_compute_stresses_data_lists

.. autofunction:: GWProt.stress_comparison.Stress_Comparison.FGW_compute_stresses_dict




.. autofunction:: GWProt.stress_comparison.Stress_Comparison.get_GW_dmat

This method returns the distance matrix computed with ``GW_compute_stresses`` or ``FGW_compute_stresses``, so must be run after one of them. The ijth entry is the (F)GW distance between the ith and jth entries in ``self.prot_list``.


.. autofunction:: GWProt.stress_comparison.Stress_Comparison.raw_transferred_stresses





We also have several helper methods not part of the ``Stress_Comparison`` class but can be useful for analyzing stresses:

.. autofunction:: GWProt.stress_comparison.normalize_stress_dict

.. autofunction:: GWProt.stress_comparison.get_AP_scores

.. autofunction:: GWProt.stress_comparison.mean_AP_scores


.. autofunction:: GWProt.stress_comparison.single_threshold_AP_score

.. autofunction:: GWProt.stress_comparison.avgd_single_threshold_AP_score

.. autofunction:: GWProt.stress_comparison.z_single_threshold_AP_score


.. autofunction:: GWProt.stress_comparison.get_percentile_of_dict

