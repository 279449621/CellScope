CellScope.cm.best_map
==========================

CellScope.cm. ``best_map`` (**L1** : :class:`np.ndarray`, **L2** : :class:`np.ndarray`)

``best_map`` is a utility function that leverages the Hungarian algorithm [Mu57]_ to find the optimal mapping between two sets of numerical labels.
It adjusts the predicted labels (`L2`) to best align with the true labels (`L1`), thus minimizing the mismatch between the two sets. 

.. rubric:: Parameters

- **L1** (:class:`np.ndarray`): 

   A numerical label vector with a length of `num_cell`, representing the true labels for each cell.

- **L2** (:class:`np.ndarray`): 
 
   A numerical label vector with a length of `num_cell`, representing the predicted labels for each cell.

.. rubric:: Return

- **new_L2** (:class:`ndarray`): 
 
   A numerical label vector of size `(num_cell)`, obtained by permuting `L2` using the Hungarian algorithm, making it align with `L1` as closely as possible.

.. rubric:: Notes

The `best_map` function utilizes the Hungarian algorithm  to perform label alignment, which is commonly used in clustering tasks to measure the similarity between two label assignments. The method ensures that the predicted labels (L2) are reordered to match the true labels (L1), making it a valuable step in cluster evaluation.
