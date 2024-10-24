CellScope.cs.findCenters
==========================

CellScope.cs. ``findCenters`` (**rho** : :class:`np.ndarray`, **delta** : :class:`np.ndarray`)

`findCenters` adaptively identifies multiple cluster centers based on the product of local density (`rho`) 
and the distance to the nearest higher-density cell (`delta`).

.. rubric:: Parameters

- **rho** (:class:`ndarray`): 
   A vector of length `n_cell`, representing the local density for each cell.

- **delta** (:class:`ndarray`): 
   A vector of Euclidean distances between each cell and the nearest higher-density cell.



.. rubric:: Return

- **centers** (:class:`ndarray`): 
   A vector containing the indices of the selected cluster centers, with the length equal to the number of centers.

- **deltarho** (:class:`ndarray` | :class:`csr_matrix`): 
   A vector of length `n_cell`, representing the product of the normalized `rho` (local density) and `delta` (distance to higher-density cells).
