ClusteringMetrics: ``cm``
==========================

Computing the clustering accuracy between predict cluster and real type.

.. list-table::

   * - ``CellScope.cm.process_labels``
     - Convert string-based labels into numerical labels for further processing.
   * - ``CellScope.cm.best_map``
     - Map predicted labels to the corresponding ground truth labels.
   * - ``CellScope.cm.calculate_metrics``
     - Compute clustering accuracy metrics between predicted and true labels.


.. toctree::
   :maxdepth: 2

   CellScope.cm.process_labels
   CellScope.cm.best_map
   CellScope.cm.calculate_metrics