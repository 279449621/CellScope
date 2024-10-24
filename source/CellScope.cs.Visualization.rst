CellScope.cs.Visualization
==========================

CellScope.cs. ``Visualization`` (**fea** : :class:`np.ndarray`, **Visualization_Method** : :class:`str` = "UMAP")

The `Visualization` function generates a two-dimensional representation of high-dimensional data, offering UMAP [MMH18]_, PCA [AW10]_, and t-SNE [VH08]_ as available visualization methods.


.. rubric:: Parameters

- **fea** (:class:`np.ndarray`): 

   The result from either the first or second Manifold Fitting step, with dimensions (num_cell, num_Selected_Gene).

- **Visualization_Method** (:class:`str`, optional, default="UMAP"): 

   The method used for visualizing the data. The default is UMAP, but available options include PCA, tSNE, and UMAP.


.. rubric:: Return

- **Y** (:class:`ndarray`): 

   A matrix of size (n_cell, 2), representing the visualization result with two-dimensional coordinates for each cell.
