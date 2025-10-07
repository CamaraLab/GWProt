LGD Comparison
========================


.. module:: GWProt.stress_comparison

.. autoclass:: GWProt.stress_comparison.LGD_Comparison

This module stores the transport plans and local geometric distortion (LGD) values between all pairs of proteins, which can be memory intensive. Setting ``RAM = False`` saves these to files, but significantly slows down computations due to file I/O.


**Calculating LGDs and Distances**
------------------------------------------

These methods run all pairwise computations. For each pair of proteins, they run GW/FGW, store the distance, transport plan, and associated LGD values. One of these must be run before any of the latter methods are called. 

As these involve a large number of computations, they can be time consuming on large datasets.

.. autofunction:: GWProt.stress_comparison.LGD_Comparison.GW_compute_lgd

.. autofunction:: GWProt.stress_comparison.LGD_Comparison.FGW_compute_lgd_data_lists

.. autofunction:: GWProt.stress_comparison.LGD_Comparison.FGW_compute_lgd_dict


These methods must be run after computing the LGDs:

.. autofunction:: GWProt.stress_comparison.LGD_Comparison.get_GW_dmat



.. autofunction:: GWProt.stress_comparison.normalize_lgd_dict

This method combines the raw LGD arrays (one for each protein pair) into a single LGD array for each protein.

For each protein, it returns the ``normalized_lgd`` calculated below, where the rows of ``mat`` are its different LGD arrays in the raw LGD dict:

.. code-block:: python

	a, b, c, d, e = code
	normalized_lgd = np.sum(mat**a * (np.sum(mat**b, axis=1) ** c)[:, np.newaxis] * mat.shape[1] ** d, axis=0) * mat.shape[0] ** e


**Transferring LGDs**
----------------------------------



.. autofunction:: GWProt.stress_comparison.LGD_Comparison.raw_transferred_lgd

This is in the same format as ``self.raw_lgd_dict`` so must be normalized before further use. 



**Analysis Helper Methods**
--------------------------------

We also provide helper methods not part of the ``LGD_Comparison`` class but useful for analyzing LGD values.



.. autofunction:: GWProt.stress_comparison.get_percentile_of_dict


.. autofunction:: GWProt.stress_comparison.get_AP_scores







