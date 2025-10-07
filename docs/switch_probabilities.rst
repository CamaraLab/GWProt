Switch Probabilities
========================

These methods calculate and help visualize the probability of pairs of residues having their order switched when one protein is aligned to another using the Gromov-Wasserstein (GW) framework. This analysis is most meaningful when the proteins are morphologically similar (typically with a GW distance under 2).

**Switch Probability Calculation**
----------------------------------

.. autofunction:: GWProt.switch_probabilities.get_switch_probabilities

This function uses sparse matrices for efficiency, but may encounter issues if the transport plan ``T`` has many non-zero entries. In practice, this is rarely a problem.

**Visualization and Utilities**
------------------------------

.. autofunction:: GWProt.switch_probabilities.visualize_switch_probabibilities

.. autofunction:: GWProt.switch_probabilities.preprocess

.. autofunction:: GWProt.switch_probabilities.max_rectangle_diagonal

The maximal rectangle algorithm is based on the solution to the classic maximal rectangle problem, `as written up by David Vandevoorde <https://drdobbs.com/database/the-maximal-rectangle-problem/184410529>`_.

Note: While these methods are related to GW alignment and LGD analysis, they focus on the order and mapping of residues rather than direct geometric distortion.
