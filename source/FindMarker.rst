FindMarker: ``fm``
==========================

Identifies and visualizes marker genes between two specific cell clusters.

.. list-table::

   * - ``CellScope.fm.FindMarker``
     - Find the Marker Gene between two clusters.

   * - ``CellScope.fm._diff_pct``
     - Calculate the difference in the percentage of gene expression between two clusters.

   * - ``CellScope.fm._diff_mean``
     - Calculate the difference in the mean expression of genes between two clusters.

   * - ``CellScope.fm._fold_change``
     - Compute the fold change in gene expression between two clusters, adjusting for cases with zero values.

   * - ``CellScope.fm._wilcoxon_parallel``
     - Perform the Wilcoxon rank-sum test in parallel to identify statistically significant marker genes between two clusters.

   * - ``CellScope.fm._t_test_parallel``
     - Perform a t-test in parallel to identify statistically significant marker genes between two clusters.

     
.. toctree::
   :maxdepth: 2
  
   CellScope.fm.FindMarker
   CellScope.fm._diff_pct
   CellScope.fm._diff_mean
   CellScope.fm._fold_change
   CellScope.fm._wilcoxon_parallel
   CellScope.fm._t_test_parallel
