CellScope: ``cs``
==========================

The main process of CellScope includes: Data normalization, First Manifold Fitting, Second step of Manifold Fitting and Hierarchical clustering based on similarity graph.

.. list-table::

   * - ``CellScope.cs.Normalization``
     - Normalize the raw scRNA data.
   * - ``CellScope.cs.findCenters``
     - Adaptive selection of cell types' centers
   * - ``CellScope.cs.Manifold_Fitting_1``
     - Perform manifold fitting to identify Signal Space.
   * - ``CellScope.cs.Manifold_Fitting_2``
     - Perform second manifold fitting.
   * - ``CellScope.cs.GraphCluster``
     - Perform graph-based clustering on the selected features.
   * - ``CellScope.cs.Visualization``
     -  Visualize the feature data using specified visualization method


.. toctree::
   :maxdepth: 2

   CellScope.cs.Normalization
   CellScope.cs.findCenters
   CellScope.cs.Manifold_Fitting_1
   CellScope.cs.Manifold_Fitting_2
   CellScope.cs.GraphCluster
   CellScope.cs.Visualization