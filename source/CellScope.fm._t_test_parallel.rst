CellScope.fm._t_test_parallel
==============================

CellScope.fm. ``_t_test_parallel`` (**fea_1** : :class:`np.ndarray`, **fea_2** : :class:`np.ndarray`, **selected_number** : :class:`int`)

The ``_t_test_parallel`` function performs a t-test [KT15]_ in parallel to compare gene expression between two clusters and selects the top genes based on the specified number of genes with the most significant p-values.

.. rubric:: Parameters

- **fea_1** (:class:`np.ndarray`): 

   Feature matrix for the first cluster.

- **fea_2** (:class:`np.ndarray`): 

   Feature matrix for the second cluster.

- **selected_number** (:class:`int`): 

   Number of top markers to select.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the top markers based on t-test.

- **p_values** (:class:`ndarray`): 

   P-values of the selected markers.
