CellScope.fm._diff_pct
=======================

CellScope.fm. ``_diff_pct`` (**fea_1** : :class:`np.ndarray`, **fea_2** : :class:`np.ndarray`, **selected_number** : :class:`int`)

The ``_diff_pct`` function calculates the differential percentage of gene expression between two clusters, 
selecting the top genes with the highest percentage differences according to the specified number of genes.

.. rubric:: Parameters

- **fea_1** (:class:`np.ndarray`): 

   Feature matrix for the first cluster.

- **fea_2** (:class:`np.ndarray`): 

   Feature matrix for the second cluster.

- **selected_number** (:class:`int`): 

   Number of top markers to select.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the top markers based on differential percentage.

- **marker_values** (:class:`ndarray`): 

   Differential percentage values of the selected markers.
