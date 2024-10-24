CellScope.cs.Manifold_Fitting_2
===============================

**CellScope.cs.** ``Manifold_Fitting_2`` (**fea_selected** : :class:`np.ndarray`, **num_neighbor** : :class:`int` = 5, **fitting_prop** : :class:`float` = 0.05, **coeff** : :class:`float` = 0.1, **op_Outlier** : :class:`bool` = False)

``Manifold_Fitting_2`` reduces the impact of ambiguous manifold boundaries on downstream analysis by pulling low-density cells closer to higher-density ones. This process involves identifying low-density cells, adjusting potential outliers, and further fitting a subset of cells to minimize noise.

.. rubric:: Parameters

- **fea_selected** (:class:`np.ndarray`): 

   The selected feature matrix from the first step of Manifold Fitting, with dimensions (n_cell, num_Selected_Gene).

- **num_neighbor** (:class:`int`, optional, default=5): 

   The number of nearest neighbors to consider when calculating local density in the second manifold fitting step.

- **fitting_prop** (:class:`float`, optional, default=0.05): 

   The proportion of cells selected for the second manifold fitting step, representing a fraction of the total cells. It must be a float between 0 and 1.

- **coeff** (:class:`float`, optional, default=0.1): 

   The coefficient used to update the fitting cells by averaging with their nearest neighbors. It must be a float between 0 and 1.

- **op_Outlier** (:class:`bool`, optional, default=False): 

   Specifies whether to identify and remove outliers among the fitting cells during the second manifold fitting step.

.. rubric:: Return

- **fea_selected** (:class:`np.ndarray`): 

   The updated feature matrix after the second step of manifold fitting, with dimensions (n_remaining_cells, num_Selected_Gene), where n_remaining_cells is the number of cells left after fitting and optional outlier removal.

- **fitting_index** (:class:`ndarray`): 

   The indices of the cells selected for the second manifold fitting step, representing the most significant fitting cells.

- **index** (:class:`ndarray`): 

   The indices of the remaining cells after removing outliers, representing the final set of cells included in the adjusted feature matrix.
