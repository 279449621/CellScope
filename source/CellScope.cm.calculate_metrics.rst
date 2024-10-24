CellScope.cm.calculate_metrics
====================================================

CellScope.cm. ``calculate_metrics`` (**true_labels** : :class:`np.ndarray`, **pred_labels** : :class:`np.ndarray`)

`CellScope.cm.calculate_metrics` is a function designed to evaluate the consistency between clustering results and true labels. 
By inputting true labels and predicted labels, the function computes a series of standard clustering evaluation metrics, including Accuracy, Normalized Mutual Information (NMI), Adjusted Rand Index (ARI), F1 Score, and Jaccard Index. 

.. rubric:: Parameters

- **true_labels** (:class:`np.ndarray`): 

   A string-based label vector with a length of `num_cell`, representing the true labels for each cell.

- **pred_labels** (:class:`np.ndarray`): 
  
   A string-based label vector with a length of `num_cell`, representing the predicted labels for each cell.

.. rubric:: Return

- **metrics** (:class:`dict`): 
   
   A dictionary containing the clustering accuracy metrics between the `true_labels` and the `pred_labels`:

  - **"Accuracy"** (:class:`float`): 
   The proportion of correctly predicted labels, calculated as the ratio of correct predictions to total samples.

  - **"NMI"** (:class:`float`): 
    
   Normalized Mutual Information, which measures the amount of information shared between true and predicted labels, normalized to the range [0, 1]. Higher values indicate better clustering performance.

  - **"ARI"** (:class:`float`): 
   
   Adjusted Rand Index, assessing the similarity between true and predicted clustering labels. It adjusts for chance, with values ranging from -1 to 1, where 1 indicates perfect agreement.

  - **"F1 Score"** (:class:`float`): 
    
   The weighted harmonic mean of precision and recall, considering the balance between the two metrics. It ranges from 0 to 1, where 1 indicates perfect precision and recall.

  - **"Jaccard"** (:class:`float`): 
    
   Jaccard Index, measuring the intersection-over-union of true and predicted labels, ranging from 0 to 1. Higher values indicate better overlap between the two sets of labels.

