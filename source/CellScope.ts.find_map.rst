CellScope.ts.find_map
=====================

CellScope.ts. ``find_map`` (**label1** : :class:`np.ndarray`, **label2** : :class:`np.ndarray`)

The ``find_map`` function creates a mapping matrix to identify overlaps between two sets of labels, 
helping to determine the correspondence between clusters across different labeling schemes.

.. rubric:: Parameters

- **label1** (:class:`np.ndarray`): 

   First set of clustering labels.

- **label2** (:class:`np.ndarray`): 

   Second set of clustering labels.

.. rubric:: Return

- **mapping** (:class:`np.ndarray`): 

   Mapping matrix indicating the overlap between two sets of labels.
