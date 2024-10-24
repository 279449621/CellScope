CellScope.cs.Manifold_Fitting_1
===============================

CellScope.cs. ``Manifold_Fitting_1`` (**fea** : :class:`np.ndarray` | :class:`csr_matrix`, **num_pca** : :class:`int` = 100, **num_Selected_Gene** : :class:`int` = 500, **knn** : :class:`int` = 20, **num_center** : :class:`int` = 0, **random_seed** : :class:`int` = 83)

``Manifold_Fitting_1`` fits the data onto a low-dimensional manifold by identifying significant genes. 
This process involves dimensionality reduction through PCA [AW10]_, selecting manifold seeds [RL14]_, identifying high-confidence cliques, and selecting key genes that define the signal space based on these cliques.


.. rubric:: Parameters

- **fea** (:class:`np.ndarray` | :class:`csr_matrix`): 

   The feature matrix of shape (n_cell, n_gene), where rows correspond to cells and columns to genes. Supports both dense (NumPy array) and sparse (CSR matrix) formats.

- **num_pca** (:class:`int`, optional, default=100): 

   The number of principal components to retain after performing PCA on `fea`. This parameter determines the dimensionality of the PCA-transformed data.

- **num_Selected_Gene** (:class:`int`, optional, default=500): 

   The number of selected genes representing the signal space, based on the lowest p-values.

- **knn** (:class:`int`, optional, default=20): 

   The number of nearest neighbors used to compute the local density (`rho`) and distance to higher density neighbors (`delta`).

- **num_center** (:class:`int`, optional, default=0): 

   The number of cluster centers to select. If set to 0, the cluster centers are automatically determined using the `findCenters` function based on the product of `rho` and `delta`.

- **random_seed** (:class:`int`, optional, default=83): 

   The random seed for ensuring reproducibility of PCA results and Truncated SVD.

.. rubric:: Return

- **fea_selected** (:class:`ndarray` | :class:`csr_matrix`): 

   A matrix of shape (n_cell, num_Selected_Gene), representing the gene-selected data matrix.

- **significant_features_index** (:class:`ndarray`): 

   A vector containing the indices of the selected genes, corresponding to the most significant features based on p-values.

- **id_max_deltarho** (:class:`ndarray`): 

   An array of cluster center indices, representing the cells identified as cluster centers.
