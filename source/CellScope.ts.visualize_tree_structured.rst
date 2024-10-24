CellScope.ts.visualize_tree_structured
======================================

CellScope.ts. ``visualize_tree_structured`` (**Y_initial**, **label_step0**, **Y_1**, **Title_1**, **Y_all**, **Title_all**,  
**index_1**, **index_all**, **step0**, **step1**, **T**, **save_fig**=False, **save_path**='')

The ``visualize_tree_structured`` function creates a series of UMAP-based scatter plots representing hierarchical 
clustering across different levels. 
It visually captures initial clusters, subclusters, and refined clusters, 
supporting tree-structured analysis of cell data. Optionally, it saves the figures to a specified directory.

.. rubric:: Parameters

- **Y_initial** (:class:`np.ndarray`): 

   Initial 2D UMAP embedding of all cells.

- **label_step0** (:class:`np.ndarray`): 

   Clustering labels at the starting level (`step0`).

- **Y_1** (:class:`list` of :class:`np.ndarray`): 

   List of UMAP embeddings for clusters at the first level.

- **Title_1** (:class:`list`): 

   List of titles for the first-level clusters.

- **Y_all** (:class:`list` of :class:`np.ndarray`): 

   List of UMAP embeddings for clusters at all levels.

- **Title_all** (:class:`list`): 

   List of titles for clusters at all levels.

- **index_1** (:class:`list`): 

   List of indices corresponding to cells in each cluster at the first level.

- **index_all** (:class:`list`): 

   List of indices corresponding to cells in each cluster at all levels.

- **step0** (:class:`int`): 

   Final starting level for the tree-structured visualization.

- **step1** (:class:`int`): 

   Final ending level for the tree-structured visualization.

- **T** (:class:`np.ndarray`): 

   Matrix of clustering labels for different levels of clustering.

- **save_fig** (:class:`bool`, optional, default=False): 

   If True, saves the figures to the specified path.

- **save_path** (:class:`str`, optional, default=''): 

   Directory path to save the figures. Required if `save_fig=True`.

.. rubric:: Return

- This function does not return any values. It displays a grid of scatter plots showing UMAP embeddings for different clusters at different levels.
