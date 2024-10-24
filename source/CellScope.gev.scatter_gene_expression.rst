CellScope.gev.scatter_gene_expression
=====================================

CellScope.gev. ``plot_mean_expression_heatmap`` (**fea** : :class:`np.ndarray`, **Y** : :class:`np.ndarray`, **marker_gene_indices** : :class:`list`, **marker_gene_name** : :class:`list`, **figsize** : :class:`tuple` = (15, 9), **subplot_size** : :class:`tuple` = None, **save_fig** : :class:`bool` = False, **save_path** : :class:`str` = None)

The ``scatter_gene_expression`` function generates scatter plots for 
visualizing the expression of specified marker genes in a two-dimensional space defined by `Y`. 
It supports customization of figure size, subplot layout, and saving of the resulting figure.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray` | :class:`csr_matrix`): 

   The feature matrix for gene expression.

- **Y** (:class:`np.ndarray`): 

   2D array of coordinates for scatter plot.

- **marker_gene_indices** (:class:`list`): 

   List of indices for marker genes to be visualized.

- **marker_gene_name** (:class:`list`): 

   List of marker gene names corresponding to the indices.

- **figsize** (:class:`tuple`, optional, default=(15, 9)): 

   Size of the figure.

- **subplot_size** (:class:`tuple`, optional): 

   Number of rows and columns for subplots.

- **save_fig** (:class:`bool`, optional, default=False): 

   If set to True, saves the figure to the specified path.

- **save_path** (:class:`str`, optional): 

   File path to save the figure. Required if `save_fig=True`.

.. rubric:: Return

- This function does not return any values. It directly displays a scatter plot of marker gene expression.
