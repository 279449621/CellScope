CellScope.fm.FindMarker
========================

CellScope.fm. ``FindMarker`` (**fea** : :class:`np.ndarray` | :class:`csr_matrix`, **cluster1** : :class:`np.ndarray`, **cluster2** : :class:`np.ndarray`, **selected_number** : :class:`int` = 5, **selected_method** : :class:`str` = 'diff pct')

The ``FindMarker`` function identifies marker genes that differentiate two clusters. 
It supports various selection methods, such as 'diff pct', 'diff mean', 'FC' (Fold Change), 'Wilcoxon' [GE65]_, and 't-test' [KT15]_, and returns the top specified number of marker genes.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray` | :class:`csr_matrix`): 

   The feature matrix of shape (n_cell, n_gene), where rows correspond to cells and columns to genes.

- **cluster1** (:class:`np.ndarray`): 

   Indices of cells belonging to the first cluster.

- **cluster2** (:class:`np.ndarray`): 

   Indices of cells belonging to the second cluster.

- **selected_number** (:class:`int`, optional, default=5): 

   Number of top markers to select.

- **selected_method** (:class:`str`, optional, default='diff pct'): 

   The method for selecting markers. Options include 'diff pct', 'diff mean', 'FC', 'Wilcoxon', and 't-test'.

.. rubric:: Return

- **marker_indices** (:class:`ndarray`): 

   Indices of the selected marker genes.

- **marker_values** (:class:`ndarray`): 

   Corresponding values of the selected markers based on the chosen method.
