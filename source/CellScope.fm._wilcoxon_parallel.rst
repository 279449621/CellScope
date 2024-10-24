CellScope.fm._wilcoxon_parallel
================================

CellScope.fm. ``_wilcoxon_parallel`` (**fea_1** : :class:`np.ndarray`, **fea_2** : :class:`np.ndarray`, **selected_number** : :class:`int`)

The ``_wilcoxon_parallel`` function conducts a Wilcoxon rank-sum test [GE65]_ in parallel to evaluate differential gene expression between two clusters, selecting the specified number of top genes based on the smallest p-values.

.. rubric:: Parameters

- **fea_1** (:class:`np.ndarray`): 

   Feature matrix for the first cluster.

- **fea_2** (:class:`np.ndarray`): 

   Feature matrix for the second cluster.

- **selected_number** (:class:`int`): 

   Number of top markers to select.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the top markers based on Wilcoxon rank-sum test.

- **p_values** (:class:`ndarray`): 

   P-values of the selected markers.
