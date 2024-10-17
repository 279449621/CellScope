Manifold Fitting and Clustering
============================

This dataset consists of 4,714 cells from the Red Nucleus region of the Midbrain, as part of the Human Brain Cell Atlas. It is freely available on the CELLxGENE website in h5ad format and can be downloaded via the following link: https://datasets.cellxgene.cziscience.com/5488ff72-58ed-4f0d-913c-1b6d4d8412b1.h5ad.

Data Download and Reading
-------------------------

.. code-block:: python

    import requests
    import anndata
    url = "https://datasets.cellxgene.cziscience.com/5488ff72-58ed-4f0d-913c-1b6d4d8412b1.h5ad"
    file_path = "Siletti-1.h5ad"
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    adata = anndata.read_h5ad("Siletti-1.h5ad")
    fea_raw = adata.X
    cell_types = adata.obs['cell_type']
    label = np.array(cell_types)
    data_type = 'csr'

Normalization
-------------------------

.. code-block:: python

    import CellScope 
    CS = CellScope()

.. code-block:: python

    fea_raw, fea_log, fea = CS.Normalization(fea_raw, data_type)

Manifold Fitting Step 1
-------------------------
.. code-block:: python

    fea_Fitting_1, p_values, Signal_Space = CS.Manifold_Fitting_1(fea, num_pca=100, num_Selected_Gene=500, knn=20, num_center=0)

Manifold Fitting Step 2  
-------------------------
.. code-block:: python

    fea_Fitting_2,fitting_index,index = CS.Manifold_Fitting_2(fea_Fitting_1, num_neighbor=5, fitting_prop=0.05, coeff=0.1, op_Outlier=False)

Manifold Embedding in Graph and Clustering
-------------------------
.. code-block:: python

    T_all_1 = CS.GraphCluster(fea_Fitting_1, metric='ST', num_cell_thre=100000, index=[])
    T_all_2 = CS.GraphCluster(fea_Fitting_2, metric='ST', num_cell_thre=100000, index=[])
