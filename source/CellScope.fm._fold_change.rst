CellScope.fm._fold_change
==========================

CellScope.fm. ``_fold_change`` (**fea_1** : :class:`np.ndarray`, **fea_2** : :class:`np.ndarray`, **selected_number** : :class:`int`)

The ``_fold_change`` function computes the fold change of gene expression between two clusters and selects the top genes with the largest fold changes based on the specified number of genes.

.. rubric:: Parameters

- **fea_1** (:class:`np.ndarray`): 

   Feature matrix for the first cluster.

- **fea_2** (:class:`np.ndarray`): 

   Feature matrix for the second cluster.

- **selected_number** (:class:`int`): 

   Number of top markers to select.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the top markers based on fold change.

- **FC_values** (:class:`ndarray`): 

   Fold change values of the selected markers.
