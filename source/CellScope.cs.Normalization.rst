CellScope.cs.Normalization
==========================

CellScope.cs. ``Normalization`` (**self**, **fea_raw** : :class:`np.ndarray` | :class:`csr_matrix`, **data_type** : :class:`str`)

``Normalization`` performs logarithmic transformation on the raw data and applies row normalization, ensuring that the norm of each cell equals 1.

.. rubric:: Parameters

- **fea_raw** (:class:`ndarray` | :class:`csr_matrix`): 
    The data matrix of shape (n_cell, n_gene), where rows correspond to cells and columns to genes.


.. rubric:: Return

- **fea_raw** (:class:`ndarray` | :class:`csr_matrix`): 
   A data matrix with dimensions (n_cell, n_gene), where rows represent cells and columns represent genes.

- **fea_log** (:class:`ndarray` | :class:`csr_matrix`): 
   A transformed version of `fea_raw`. If the maximum value in `fea_raw` is less than 1000, then `fea_log = log2(fea_raw + 1)`; otherwise, `fea_log = fea_raw`.

- **fea** (:class:`ndarray` | :class:`csr_matrix`): 
   A row-normalized version of `fea_log`, where the norm of each row is scaled to 1.
