CellScope.cs.GraphCluster
==========================

CellScope.cs. ``GraphCluster`` (**fea_selected** : :class:`np.ndarray`, **metric** : :class:`str` = None, 
**num_cell_thre** : :class:`int` = 100000, **index** : :class:`list` = [])

`GraphCluster` constructs a cell similarity graph using UMAP’s fuzzy simplicial set [MMH18]_ and performs hierarchical clustering to group cells, 
returning clustering results at multiple levels of granularity. 
The function is adaptable to large-scale datasets by sampling a subset of cells to enhance computational efficiency.

.. rubric:: Parameters

- **fea_selected** (:class:`np.ndarray`): 

   A 2D array representing the selected features for clustering, where rows correspond to cells and columns to genes.

- **metric** (:class:`str`, optional, default=None): 

   The distance metric used for clustering. If not specified, it is automatically determined based on the dataset size: Euclidean is selected when the dataset has fewer than 10,000 cells, and Jaccard is used otherwise.

- **num_cell_thre** (:class:`int`, optional, default=100000): 

   The threshold for the number of cells. If the dataset exceeds this threshold, a subset of cells is selected for clustering.

- **index** (:class:`list`, optional, default=[]): 

   An optional list of indices for selecting specific cells. If left empty and the number of cells exceeds `num_cell_thre`, a random subset of cells will be selected.

- **random_seed** (:class:`int`, optional, default=83): 

   The random seed for ensuring reproducibility when constructing the graph and selecting the subset of indices.

.. rubric:: Return

- **T_all** (:class:`np.ndarray`): 

   A matrix of shape (num_cell, 49), where each column represents the clustering results from different steps of the agglomerative clustering process.




