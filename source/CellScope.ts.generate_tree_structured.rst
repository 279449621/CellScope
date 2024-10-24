CellScope.ts.generate_tree_structured
=====================================

CellScope.ts. ``generate_tree_structured`` (**fea**, **T**, **step0** = None, **step1** = None, **cell_type** = None)

The ``generate_tree_structured`` function constructs a tree-like representation of cell clustering, 
identifying hierarchical relationships between clusters across different levels, 
and generates UMAP-based visualizations for each cluster.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray`): 

   The feature matrix of shape (n_cells, n_genes).

- **T** (:class:`np.ndarray`): 

   A matrix of clustering labels for different levels of clustering.

- **step0** (:class:`int`, optional, default=None): 

   The starting level for the tree-structured visualization. If None, it is automatically determined.

- **step1** (:class:`int`, optional, default=None): 

   The ending level for the tree-structured visualization. If None, it is determined based on the highest ARI score.

- **cell_type** (:class:`np.ndarray`, optional, default=None): 

   A 1D array of cell type labels. Required if `step1` is not provided.

.. rubric:: Return

- **Y_initial** (:class:`np.ndarray`): 

   The initial 2D UMAP embedding of all cells.

- **label_step0** (:class:`np.ndarray`): 

   Clustering labels at the starting level (`step0`).

- **Y_1** (:class:`list` of :class:`np.ndarray`): 

   A list of UMAP embeddings for clusters at the first level.

- **Title_1** (:class:`list`): 

   A list of titles for the first-level clusters.

- **Y_all** (:class:`list` of :class:`np.ndarray`): 

   A list of UMAP embeddings for clusters at all levels.

- **Title_all** (:class:`list`): 

   A list of titles for clusters at all levels.

- **index_1** (:class:`list`): 

   A list of indices corresponding to cells in each cluster at the first level.

- **index_all** (:class:`list`): 

   A list of indices corresponding to cells in each cluster at all levels.

- **step0** (:class:`int`): 

   The final starting level for the tree-structured visualization.

- **step1** (:class:`int`): 

   The final ending level for the tree-structured visualization.
