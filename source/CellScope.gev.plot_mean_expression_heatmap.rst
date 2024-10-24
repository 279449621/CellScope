CellScope.gev.plot_mean_expression_heatmap
==========================================

CellScope.gev. ``plot_mean_expression_heatmap`` (**fea** : :class:`np.ndarray`, **marker_gene_indices** : :class:`list`, **marker_gene_name** : :class:`list`, **cell_type** : :class:`np.ndarray`, **save_fig** : :class:`bool` = False, **save_path** : :class:`str` = None)

The ``plot_mean_expression_heatmap`` function creates a heatmap displaying the mean expression levels of marker genes 
across different cell types. 
It allows for the visualization of gene expression patterns and 
provides an option to save the heatmap as a file for further analysis or reporting.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray` | :class:`csr_matrix`): 

   Feature matrix for gene expression.

- **marker_gene_indices** (:class:`list`): 

   List of indices for marker genes.

- **marker_gene_name** (:class:`list`): 

   List of marker gene names.

- **cell_type** (:class:`np.ndarray`): 

   Array of cell type labels.

- **save_fig** (:class:`bool`, optional, default=False): 

   If set to True, saves the figure to the specified path.

- **save_path** (:class:`str`, optional): 

   File path to save the figure. Required if `save_fig=True`.

.. rubric:: Return

- This function does not return any values. It displays a heatmap of mean expression for marker genes across cell types.
