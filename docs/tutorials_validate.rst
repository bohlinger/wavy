Validate the collocated time series
###################################

Having collocated a quick validation can be performed using the validationmod. validation_specs.yaml can be adjusted.

.. code-block:: python3

   >>> val_dict = cco_raw.validate_collocated_values()

   # ---
   Validation stats
   # ---
   Correlation Coefficient: 0.95
   Mean Absolute Difference: 0.22
   Root Mean Squared Difference: 0.27
   Normalized Root Mean Squared Difference: 0.08
   Debiased Root Mean Squared Difference: 0.24
   Bias: -0.13
   Normalized Bias: -0.04
   Scatter Index: 8.05
   Model Activity Ratio: 0.95
   Mean of Model: 3.02
   Mean of Observations: 3.14
   Number of Collocated Values: 72

The entire validation dictionary will then be in val_dict and can be used further in the code.
