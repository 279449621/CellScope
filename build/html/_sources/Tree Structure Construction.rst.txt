Tree Structure Construction
============================

This dataset consists of 4714 cells from the Red Nucleus region of the Midbrain, as part of the Human Brain Cell Atlas. It is freely available on the CELLxGENE website in h5ad format and can be downloaded via the following link: https://datasets.cellxgene.cziscience.com/5488ff72-58ed-4f0d-913c-1b6d4d8412b1.h5ad.

Hierarchical Clustering of the Data
-------------------------

.. code-block:: python

    import anndata
    import CellScope
    import CellScope.CellScope as CS
    adata = anndata.read_h5ad("Siletti-1.h5ad")
    fea_raw = adata.X
    cell_types = adata.obs['cell_type']
    label = np.array(cell_types)
    data_type = 'csr'
    CS = CellScope()
    fea_raw, fea_log, fea = CS.Normalization(fea_raw, data_type)
    fea_Fitting_1, p_values, Signal_Space = CS.Manifold_Fitting_1(fea, num_pca=100, num_Selected_Gene=500, knn=20, num_center=0)
    fea_Fitting_2,fitting_index,index = CS.Manifold_Fitting_2(fea_Fitting_1, num_neighbor=5, fitting_prop=0.05, coeff=0.1, op_Outlier=False)
    T_all_1 = CS.GraphCluster(fea_Fitting_1, metric='ST', num_cell_thre=100000, index=[])
    T_all_2 = CS.GraphCluster(fea_Fitting_2, metric='ST', num_cell_thre=100000, index=[])


Tree Structure Visualization
-------------------------

.. code-block:: python

    import CellScope.TreeStructured as TS
    label = T_all_1[:,8]
    TS = TS(fea_Fitting_1, T_all_1, cell_type = label,step0 = 'ST', step1=8)
    Y_initial, label_step0, Y_1, Title_1, Y_all, Title_all, index_1, index_all, step0, step1 = TS.visualize_tree_structured()

.. image:: _static/TreeStructure.png

Marker Gene Between Each clusters
-------------------------

.. code-block:: python

    import CellScope.FindMarker as FM
    import CellScope.GeneExpressionVisualizer as GEV
    FM = FM(fea_log)

    n = Y_initial.shape[0]
    selected_number = 5
    Gene_Name = adata.var['feature_name']
    Gene_Name = np.array(Gene_Name)
    for ii in range(len(Title_1)):
        cluster1 = index_1[ii]
        cluster2 = np.setdiff1d(range(n),cluster1)
        Marker_Gene_Index,_ = FM.differential_expression(cluster1, cluster2, selected_number, selected_method='diff pct')
        Marker_Gene_Name = Gene_Name[Marker_Gene_Index]
        GEV1 = GEV(fea_log,Marker_Gene_Index, Marker_Gene_Name)
        GEV1.compare_violin_plot_between_classes(cluster1, cluster2,class_name = [Title_1[ii],'Other'])
    for ii in range(int(len(index_all)/2)):
        cluster1 = index_all[2*ii]
        cluster2 = index_all[2*ii+1]
        Marker_Gene_Index,_ = FM.differential_expression(cluster1, cluster2, selected_number, selected_method='diff pct')
        Marker_Gene_Name = Gene_Name[Marker_Gene_Index]
        GEV1 = GEV(fea_log,Marker_Gene_Index, Marker_Gene_Name)
        GEV1.compare_violin_plot_between_classes(cluster1, cluster2,class_name = [Title_all[2*ii],Title_all[2*ii+1]])

.. image:: _static/MarkerGene_Between_Cluster.png
.. image:: _static/MarkerGene_Between_SubCluster.png



