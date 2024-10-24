CellScope.fm._diff_mean
========================

CellScope.fm. ``_diff_mean`` (**fea_1** : :class:`np.ndarray`, **fea_2** : :class:`np.ndarray`, **selected_number** : :class:`int`)

The ``_diff_mean`` function calculates the differential mean expression of genes between two clusters 
and selects the top genes with the highest mean differences, based on the specified number of genes to select.

.. rubric:: Parameters

- **fea_1** (:class:`np.ndarray`): 

   Feature matrix for the first cluster.

- **fea_2** (:class:`np.ndarray`): 

   Feature matrix for the second cluster.

- **selected_number** (:class:`int`): 

   Number of top markers to select.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the top markers based on differential mean.

- **marker_values** (:class:`ndarray`): 

   Differential mean values of the selected markers.
