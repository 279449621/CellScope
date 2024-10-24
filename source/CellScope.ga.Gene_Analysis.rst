CellScope.ga.Gene_Analysis
===========================

CellScope.ga. ``Gene_Analysis`` (**fea** : :class:`np.ndarray`, **layer** : :class:`list` | :class:`np.ndarray`)

The ``Gene_Analysis`` function calculates Wasserstein distances [PZ19]_ between gene expressions in different clusters defined by `layer`, 
classifying genes as Housekeeper Gene, Type-Related Gene, or Type-Determining Gene based on their roles across layers.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray`): 

   The feature matrix of shape (n_cell, n_gene), where rows correspond to cells and columns to genes.

- **layer** (:class:`list` | :class:`np.ndarray`): 

   List or array indicating clusters. Must contain an even number of elements representing paired clusters for comparison.

.. rubric:: Return

- **Res** (:class:`np.ndarray`): 

   Wasserstein distances of shape (num_Gene, len_flow), where `len_flow` is half the length of `layer`.

- **label** (:class:`np.ndarray`): 

   Integer matrix indicating gene categories for each flow, with values: 0 ('Housekeeper Gene'), 1 ('Type-Related Gene'), or 2 ('Type-Determining Gene').

- **label_str** (:class:`np.ndarray`): 

   String matrix with the same shape as `label`, containing gene type labels.

- **flow_labels** (:class:`np.ndarray`): 

   Array indicating flow labels for genes based on their category across layers.
