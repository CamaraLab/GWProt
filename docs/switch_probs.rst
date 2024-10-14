Switch Probs
==============

These methods calculate and help visualize the probability of pairs of residues having their order switched when one protein is aligned to another. This works best when the proteins in question are morphologically similar, with a GW distance under 2 between them.

.. autofunction:: GWProt.switch_probs.get_switch_prob

This uses sparse matrices so may cause problems if ``T`` has many non-zero entries. In practice this is rarely the case.

.. autofunction:: GWProt.switch_probs.visualize_switch_probabibilities


.. autofunction:: GWProt.switch_probs.preprocess

.. autofunction:: GWProt.switch_probs.max_rectangle_diagonal

 This algorithm is based on the solution to the classic maximal rectangle problem `as written up by David Vandevoorde <https://drdobbs.com/database/the-maximal-rectangle-problem/184410529>`_.
