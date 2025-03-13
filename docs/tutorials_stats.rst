Validation metrics
##################

**wavy** offers a selection of metrics that are useful for validation. Since these computations are not costly they are all done at once on demand. The available validation metrics are:

* Correlation Coefficient
* Mean Absolute Difference
* Root Mean Squared Difference
* Normalized Root Mean Squared Difference
* Debiased Root Mean Squared Difference
* Bias
* Normalized Bias
* Scatter Index
* Model Activity Ratio
* Mean of Model
* Mean of Observations
* Number of Collocated Values

Variable names and the computation of Bias, MAD, and RMSD follow the notation of `WMO <https://library.wmo.int/doc_num.php?explnum_id=11599>`_. The correlation coefficient and scatter index follows the convention from `Wilks (2011) <https://doi.org/10.1016/B978-0-12-385022-5.00008-7>`_.

* :math:`x_f` is the forecast values of the chosen parameter
* :math:`x_v` is the observed value to evaluate against ("ground truth")
* n is the number of values to be used for verification/the number of collocated values

Bias or mean deviation (MD) is the difference between the average forecasts and observations. Positive bias indicates that the forecast is too large on average, while negative bias means that the forecast is too small on average. Bias is also known as a systematic error but gives no information about the typical magnitude of the forecast error. It is defined with the same physical dimension as the forecast and observations.

Bias = MD = :math:`\frac{1}{n}\sum_{i=1}^{n}(x_{_f}-x_{_v})_{_i}`

Normalized bias or normalized mean deviation (NMD) is the bias normalized by the sum of the observations. This version of the bias is favorable when the data are related but not strictly comparable, for example, subject to seasonal variations. One downside is that the NMD is unitless.

Normalized bias = NMD = :math:`\frac{\sum_{i=1}^{n}(x_{_f}-x_{_v})_{_i}}{\sum_{i=1}^{n}(x_{_v})_{_i}}`

Mean absolute deviation (MAD) is the arithmetic average of the absolute difference between the members of each pair of forecast and observed quantities. The MAD is zero if the difference between each pair is zero, and increases as discrepancies between the pairs become larger. It has the same units as the forecast and observations.

MAD = :math:`\frac{1}{n}\sum_{i=1}^{n}\lvert{x_{_f}-x_{_v}}\rvert_{_i}`

Root mean squared deviation (RMSD) is the square root of the average squared difference between the forecast and observation pairs. Squaring the difference produces positive terms, increasing from zero to larger values similar to the MAD. However, the RMSD will weigh the large errors more than MAD because the errors are squared. Taking the square root of the difference gives RMSD the same physical dimension as the forecast and observations. The MAD and RMSD can be interpreted as a typical magnitude for the forecast error.

RMSD = :math:`\sqrt{\frac{1}{n}\sum_{i=1}^{n}(x_{_f}-x_{_v})^{2}_{_i}}`

Normalized root mean squared deviation (NRMSD) is the RMSD, normalized by the sum of the observations squared. The NRMSD is introduced because the RMSE does not perform well if comparing models fits for different response variables or if the response variable is modified such as standardized or log-transformed. The disadvantage is that the NRMSD will lose the units associated with the forecast and observations.

NRMSD = :math:`\sqrt{\frac{\sum_{i=1}^{n}(x_{_f}-x_{_v})^{2}_{_i}}{\sum_{i=1}^{n}(x_{_f})^{2}_{_i}}}`

Debiased root mean squared deviation (DRMSD) is the square root of the difference between the mean squared deviation and the squared bias.

DRMSD = :math:`\sqrt{\frac{1}{n}{\sum_{i=1}^{n}(x_{_f}-x_{_v})^{2}_{_i}-\left(\frac{1}{n}\sum_{i=1}^{n}(x_{_f}-x_{_v})_{_i}\right)^2}}`

Pearson product-moment correlation coefficient (r) is a dimensionless, single-value measure of the association between forecast and observation. It is defined as the ratio of the covariance between forecast and observation to the product of the standard deviations of each. If the correlation coefficient is zero, there is no correlation between the variables. A correlation coefficient of 1 implies perfect positive correlation, while -1 implies perfect negative correlation.

:math:`r_{_{f,v}} = \frac{cov(x_{_f},x_{_v})}{std(f)\cdot std(v)}`

The scatter index (SI) measures the magnitude of deviation between the forecast and observations relative to the magnitude of the observations. The smaller SI, the better agreement between forecast and observations. SI can be expressed in terms of the standard deviation of the difference or the root mean squared deviation, normalized by the mean of the observations. The SI is multiplied by a hundred to return the value in percentage.

:math:`SI_{std} = \frac{std(x_{_f} - x_{_v})}{mean(x_{_v})}\cdot100`

:math:`SI_{rmsd} = \frac{rmsd}{mean(x_{_v})}\cdot100`

Model activity ratio (MAR) is the dimensionless ratio of the standard deviation of the observation to the standard deviation of the forecast. The MAR expresses the ratio of variability in each product. A value >1 implies that the numerator has a higher variability, whereas a value <1 implies the opposite.

MAR = :math:`\frac{std(x_{_v})}{std(x_{_f})}`
