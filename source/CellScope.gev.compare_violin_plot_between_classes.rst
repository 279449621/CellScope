CellScope.gev.compare_violin_plot_between_classes
=================================================

CellScope.gev. ``compare_violin_plot_between_classes`` (**fea** : :class:`np.ndarray`, **marker_gene_indices** : :class:`list`, **marker_gene_name** : :class:`list`, **cluster1** : :class:`np.ndarray`, **cluster2** : :class:`np.ndarray`, **class_name** : :class:`list` = ['1', '2'], **figsize** : :class:`tuple` = (10, 3), **save_fig** : :class:`bool` = False, **save_path** : :class:`str` = None)

The ``compare_violin_plot_between_classes`` function visualizes and compares the expression distributions of marker genes between two specified cell clusters using violin plots. 
It allows for customization of figure size, class names, and the option to save the generated plots to a file.

.. rubric:: Parameters

- **fea** (:class:`np.ndarray` | :class:`csr_matrix`): 

   Feature matrix for gene expression.

- **marker_gene_indices** (:class:`list`): 

   List of indices for marker genes.

- **marker_gene_name** (:class:`list`): 

   List of marker gene names.

- **cluster1** (:class:`np.ndarray`): 

   Indices of the first cluster.

- **cluster2** (:class:`np.ndarray`): 

   Indices of the second cluster.

- **class_name** (:class:`list`, optional, default=['1', '2']): 

   Names of the classes for labeling in the violin plot.

- **figsize** (:class:`tuple`, optional, default=(10, 3)): 

   Size of the figure.

- **save_fig** (:class:`bool`, optional, default=False): 

   If set to True, saves the figure to the specified path.

- **save_path** (:class:`str`, optional): 

   File path to save the figure. Required if `save_fig=True`.

.. rubric:: Return

- This function does not return any values. It displays a violin plot comparing gene expression between two classes.
