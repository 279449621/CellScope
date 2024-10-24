CellScope.cm.process_labels
============================

CellScope.cm. ``process_labels`` (**labels** : :class:`np.ndarray`)

``CellScope.cm.process_labels`` is a utility function that converts string labels into numerical labels. 
Additionally, the function provides a mapping between original string labels and numerical labels.

.. rubric:: Parameters

- **labels** (:class:`np.ndarray`): 
   A string-based label vector with a length of num_cell, representing the original labels for each cell.

.. rubric:: Return

- **labels** (:class:`ndarray`): 
   A numerical label vector of size (num_cell), representing the converted numerical labels.

- **label_mapping** (:class:`dict`): 
   A mapping between the original string labels and the converted numerical labels.

