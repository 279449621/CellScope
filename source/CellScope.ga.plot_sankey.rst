CellScope.ga.plot_sankey
=========================

CellScope.ga. ``plot_sankey`` (**label_str** : :class:`np.ndarray`, **width** : :class:`int` = 1100, **height** : :class:`int` = 600, **save_fig** : :class:`bool` = False, **save_path** : :class:`str` = './sankey_diagram.png')

The ``plot_sankey`` function generates a Sankey diagram to visualize gene type transitions across different layers. 
The diagram illustrates flows among Housekeeper Gene, Type-Related Gene, and Type-Determining Gene, 
based on the provided gene labels (`label_str`). 

.. rubric:: Parameters

- **label_str** (:class:`np.ndarray`): 

   A 2D NumPy array where each column contains gene labels, such as 'Housekeeper Gene', 'Type-Related Gene', or 'Type-Determining Gene'.

- **width** (:class:`int`, optional, default=1100): 

   Width of the Sankey diagram.

- **height** (:class:`int`, optional, default=600): 

   Height of the Sankey diagram.

- **save_fig** (:class:`bool`, optional, default=False): 

   If set to True, saves the figure as an image.

- **save_path** (:class:`str`, optional, default='./sankey_diagram.png'): 

   File path where the Sankey diagram is saved. Must have a valid image extension: '.png', '.pdf', or '.jpeg'.

.. rubric:: Return

- **gene_counts** (:class:`pd.DataFrame`): 

   DataFrame showing the number of 'Housekeeper Gene', 'Type-Related Gene', and 'Type-Determining Gene' at each layer.
