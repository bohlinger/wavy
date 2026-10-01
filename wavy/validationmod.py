"""
Module to organize the validation procedure
Consists mostly of functions computing validation metrics
"""

import numpy as np
from scipy import stats
import pandas as pd

# define global functions


def calc_model_activity_ratio(a, b):
    """
    computes the model activity ratio of input a (mode) and input b (obs)
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    mar = np.std(a1) / np.std(b1)
    return mar


def calc_rmsd(a, b):
    """
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    n = len(a1)
    diff2 = (a1 - b1) ** 2
    msd = diff2.sum() / n
    rmsd = np.sqrt(msd)
    return msd, rmsd


def calc_nrmsd(a, b):
    """
    Normalized root mean square deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    diff2 = (a1 - b1) ** 2
    msd = diff2.sum() / np.sum(b1**2)
    rmsd = np.sqrt(msd)
    return msd, rmsd


def calc_drmsd(a, b):
    """
    debiased root mean square deviation
    if nans exist the prinziple of marginalization is applied
    """
    a, b = np.array(a), np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    n = len(a1)
    diff2 = (a1 - b1) ** 2
    msd = diff2.sum() / n
    dmsd = msd - calc_bias(a, b) ** 2
    drmsd = np.sqrt(dmsd)
    return dmsd, drmsd


def calc_scatter_index(model, obs):
    """
    Scatter index based on rmse and on std of diff
    """
    _, rmsd = calc_rmsd(obs, model)
    stddiff = np.nanstd(obs - model)
    SIrmse = rmsd / np.nanmean(obs) * 100.0
    SIstd = stddiff / np.nanmean(obs) * 100.0
    return SIrmse, SIstd


def calc_corrcoef(a, b):
    """
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    corr = np.corrcoef(a1, b1)[1, 0]
    return corr


def calc_bias(a, b):
    """
    Bias
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    N = len(a1)
    bias = np.sum(a1 - b1) / N
    return bias


def calc_nbias(a, b):
    """
    Normalized Bias [dimensionless]
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    nbias = np.sum(a1 - b1) / np.sum(b1)
    return nbias


def calc_mad(a, b):
    """
    mean absolute deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1 = a[idx]
    b1 = b[idx]
    N = len(a1)
    mad = np.sum(np.abs(a1 - b1)) / N
    return mad


def disp_validation(valid_dict):
    """
    Print to screen validation scores.
    """
    print("\n")
    print("# ---")
    print("Validation stats")
    print("# ---")
    print("Correlation Coefficient: " + "{:0.2f}".format(valid_dict["corr"]))
    print("Mean Absolute Difference: " + "{:0.2f}".format(valid_dict["mad"]))
    print("Root Mean Squared Difference: " + "{:0.2f}".format(valid_dict["rmsd"]))
    print(
        "Normalized Root Mean Squared Difference: "
        + "{:0.2f}".format(valid_dict["nrmsd"])
    )
    print(
        "Debiased Root Mean Squared Difference: "
        + "{:0.2f}".format(valid_dict["drmsd"])
    )
    print("Bias: " + "{:0.2f}".format(valid_dict["bias"]))
    print("Normalized Bias: " + "{:0.2f}".format(valid_dict["nbias"]))
    print("Scatter Index: " + "{:0.2f}".format(valid_dict["SI"][1]))
    print("Model Activity Ratio: " + "{:0.2f}".format(valid_dict["mar"]))
    print("Mean of Model: " + "{:0.2f}".format(valid_dict["mop"]))
    print("Mean of Observations: " + "{:0.2f}".format(valid_dict["mor"]))
    print("Number of Collocated Values: " + str(valid_dict["nov"]))
    print("\n")
    pass


def validate(results_dict, boot=None):
    import numpy as np

    """
    vars in dict: np.arrays with np.nan for invalids

    produced metrics:
    mean of product --> mop
    mean of reference --> mor
    mean square difference --> msd
    number of data values --> nov
    scatter index --> SI
    """
    # date_matches = results_dict['datetime']
    if isinstance(results_dict["model_values"], list):
        model_matches = np.array(results_dict["model_values"])
    else:
        model_matches = results_dict["model_values"]
    if isinstance(results_dict["obs_values"], list):
        obs_matches = np.array(results_dict["obs_values"])
    else:
        obs_matches = results_dict["obs_values"]
    if boot is None or boot is False:
        mop = np.nanmean(model_matches)
        mor = np.nanmean(obs_matches)
        msd, rmsd = calc_rmsd(model_matches, obs_matches)
        _, nrmsd = calc_nrmsd(model_matches, obs_matches)
        _, drmsd = calc_drmsd(model_matches, obs_matches)
        nov = len(obs_matches)
        mad = calc_mad(model_matches, obs_matches)
        corr = calc_corrcoef(model_matches, obs_matches)
        bias = calc_bias(model_matches, obs_matches)
        nbias = calc_nbias(model_matches, obs_matches)
        SI = calc_scatter_index(model_matches, obs_matches)
        mar = calc_model_activity_ratio(model_matches, obs_matches)
        validation_dict = {
            "mop": mop,
            "mor": mor,
            "msd": msd,
            "nov": nov,
            "rmsd": rmsd,
            "nrmsd": nrmsd,
            "drmsd": drmsd,
            "corr": corr,
            "mad": mad,
            "bias": bias,
            "nbias": nbias,
            "SI": SI,
            "mar": mar,
        }
    elif boot is True:
        from wavy.utils import bootstr, marginalize

        reps = 1000
        newmodel, newobs, _ = marginalize(model_matches, obs_matches)
        obs_boot, boot_idx = bootstr(newobs, reps)
        print(len(obs_boot[np.isnan(obs_boot)]))
        RMSD = np.zeros(reps) * np.nan
        MSD = np.zeros(reps) * np.nan
        BIAS = np.zeros(reps) * np.nan
        CORR = np.zeros(reps) * np.nan
        for i in range(reps):
            results_dict = {
                #'date_matches':date_matches[newidx[boot_idx[:,i]]],
                "model_matches": newmodel[boot_idx[:, i]],
                "sat_matches": newobs[boot_idx[:, i]],
            }
            try:
                RMSD[i] = validate(results_dict)["rmsd"]
                MSD[i] = validate(results_dict)["mad"]
                BIAS[i] = validate(results_dict)["bias"]
                CORR[i] = validate(results_dict)["corr"]
            except IndexError as e:
                print(e)
        validation_dict = {"rmsd": RMSD, "mad": MSD, "bias": BIAS, "corr": CORR}
    return validation_dict


def linreg_evm(x, y, **kwargs):
    #  Linear regression by the maximum likelihood effective variance method.
    #
    #  K.K.Kahma 1991. Iterative solution replaced by explicit solution 1998.
    #  J.-V. Björkqvist 2020. From MATLAB to Python
    #
    #  Reference: Orear,J 1982: Least squares when both variables have
    #             uncertanties J.Am Phys 50(10)

    stdx = kwargs.get("stdx", 1)
    stdy = kwargs.get("stdy", 1)

    x0 = np.mean(x)
    y0 = np.mean(y)

    sx2 = stdx**2
    sy2 = stdy**2
    Sx2 = sum((x - x0) ** 2)
    Sy2 = sum((y - y0) ** 2)
    Sxy = sum((x - x0) * (y - y0))

    if sx2 == 0 or Sxy == 0:
        P = np.array([Sxy / Sx2])
    else:
        P = np.array(
            [
                (
                    sx2 * Sy2
                    - sy2 * Sx2
                    + np.sqrt(
                        (sy2 * Sx2) ** 2
                        - 2 * Sx2 * sy2 * sx2 * Sy2
                        + (sx2 * Sy2) ** 2
                        + 4 * Sxy**2 * sx2 * sy2
                    )
                )
                / (2 * Sxy * sx2)
            ]
        )

    P = np.append(P, y0 - P[0] * x0)
    return P


def linreg_ievm(x, y, **kwargs):
    #  Informed effective variance method.
    #  Extended evm by Patrik Bohlinger for accounting for
    #  non-stationary error variances.

    stdx = [kwargs.get("stdx", 1)]
    stdy = [kwargs.get("stdy", 1)]

    if len(stdx) == len(x) and len(stdy) == len(y):
        # stdx and stdy are vectors assumed derived from a
        # non-stationary error variance function
        pass
    elif len(stdx) == 1 and len(stdy) == 1:
        # create vectors of same length
        stdx = np.ones(len(x)) * stdx[0]
        stdy = np.ones(len(y)) * stdy[0]
    else:
        print("The format of the provided error variances is not correct!")

    x0 = np.mean(x)
    y0 = np.mean(y)

    sx2 = stdx**2
    sy2 = stdy**2

    s2 = sx2 * sy2

    Syx2 = sum(sy2 * (x - x0) ** 2)
    Sxy2 = sum(sx2 * (y - y0) ** 2)

    # prodsum1 = sum(np.sqrt(s2)*(x-x0)**2)*sum(np.sqrt(s2)*(y-y0)**2)
    prodsum2 = Syx2 * Sxy2
    prodsum = prodsum2
    # print(prodsum1)
    # print(prodsum2)

    # Sxxyy = sum(stdx*(x-x0)*stdy*(y-y0))
    Sxxyy = sum(np.sqrt(s2) * (x - x0) * (y - y0))

    Sxxy = sum(sx2 * (x - x0) * (y - y0))

    # Add what to do in case of only one value for error variance and one is 0

    P = np.array(
        [
            (Sxy2 - Syx2 + np.sqrt(Syx2**2 + Sxy2**2 - 2 * prodsum + 4 * Sxxyy**2))
            / (2 * Sxxy)
        ]
    )

    P = np.append(P, y0 - P[0] * x0)
    return P


def linreg_deming(x, y, **kwargs):
    #  Deming regression
    #  implemented by Patrik Bohlinger
    #  can account for non-stationary error variances.

    data = pd.DataFrame({"x": x, "y": y})

    stdx = [kwargs.get("stdx", 1)]
    stdy = [kwargs.get("stdy", 1)]

    if len(stdx) == len(x) and len(stdy) == len(y):
        # stdx and stdy are vectors assumed derived from a
        # non-stationary error variance function
        pass
    elif len(stdx) == 1 and len(stdy) == 1:
        # create vectors of same length
        stdx = np.ones(len(x)) * stdx[0]
        stdy = np.ones(len(y)) * stdy[0]
    else:
        print("The format of the provided error variances is not correct!")

    cov = data.cov()
    mean_x = data["x"].mean()
    mean_y = data["y"].mean()
    s_xx = cov["x"]["x"]
    s_yy = cov["y"]["y"]
    s_xy = cov["x"]["y"]
    delta = stdy**2 / stdx**2

    slope = (
        s_yy - delta * s_xx + np.sqrt((s_yy - delta * s_xx) ** 2 + 4 * delta * s_xy**2)
    ) / (2 * s_xy)
    intercept = mean_y - slope * mean_x

    P = np.append(np.mean(slope), np.mean(intercept))
    return P

def linreg_deming_homoscedastic(x, y, stdx=1.0, stdy=1.0):
    """
    Deming regression with a constant error standard deviation
    ratio stdy/stdx.

    Returns
    -------
    slope, intercept
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")

    if x.size < 2:
        raise ValueError("At least two observations are required.")

    if stdx <= 0 or stdy <= 0:
        raise ValueError("stdx and stdy must be positive.")

    xbar = np.mean(x)
    ybar = np.mean(y)

    Sxx = np.mean((x - xbar)**2)
    Syy = np.mean((y - ybar)**2)
    Sxy = np.mean((x - xbar) * (y - ybar))

    delta = (stdy / stdx)**2

    if np.isclose(Sxy, 0):
        raise ValueError("Cannot perform Deming regression when covariance is zero.")

    slope = (
        Syy - delta * Sxx
        + np.sqrt((Syy - delta * Sxx)**2 + 4 * delta * Sxy**2)
    ) / (2 * Sxy)

    intercept = ybar - slope * xbar

    return slope, intercept

import numpy as np
from scipy.optimize import minimize_scalar


def linreg_deming_heteroscedastic(x, y, stdx=1.0, stdy=1.0):
    """
    Heteroscedastic Deming regression.

    Parameters
    ----------
    x, y : array-like
        Observed x and y values.

    stdx, stdy : float or array-like
        1-sigma uncertainties in x and y.
        Can be scalars or arrays with one value per observation.

    Returns
    -------
    slope : float
        Deming regression slope.

    intercept : float
        Deming regression intercept.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one-dimensional.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    n = len(x)

    # Convert uncertainties to arrays
    stdx = np.asarray(stdx, dtype=float)
    stdy = np.asarray(stdy, dtype=float)

    if stdx.ndim == 0:
        stdx = np.full(n, stdx)

    if stdy.ndim == 0:
        stdy = np.full(n, stdy)

    if len(stdx) != n:
        raise ValueError(
            "stdx must be scalar or have the same length as x."
        )

    if len(stdy) != n:
        raise ValueError(
            "stdy must be scalar or have the same length as y."
        )

    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError(
            "x and y must contain only finite values."
        )

    if np.any(~np.isfinite(stdx)) or np.any(~np.isfinite(stdy)):
        raise ValueError(
            "stdx and stdy must contain only finite values."
        )

    if np.any(stdx <= 0) or np.any(stdy <= 0):
        raise ValueError(
            "All uncertainties must be positive."
        )

    # Convert standard deviations to variances
    varx = stdx**2
    vary = stdy**2

    def objective(slope):
        """
        Objective function after analytically minimizing
        over the intercept.
        """

        variance = vary + slope**2 * varx
        weights = 1.0 / variance

        # Optimal intercept for this slope
        intercept = np.sum(
            weights * (y - slope * x)
        ) / np.sum(weights)

        residual = y - intercept - slope * x

        return np.sum(
            weights * residual**2
        )

    # Ordinary least-squares slope as initial estimate
    slope_ols, _ = np.polyfit(x, y, 1)

    # Use a bounded scalar optimization around a
    # sufficiently broad slope range.
    slope_scale = max(
        10.0 * abs(slope_ols),
        1.0
    )

    result = minimize_scalar(
        objective,
        bounds=(-slope_scale, slope_scale),
        method="bounded",
        options={
            "xatol": 1e-12,
            "maxiter": 10000,
        },
    )

    if not result.success:
        raise RuntimeError(
            f"Deming regression failed: {result.message}"
        )

    slope = result.x

    # Calculate final intercept
    variance = vary + slope**2 * varx
    weights = 1.0 / variance

    intercept = np.sum(
        weights * (y - slope * x)
    ) / np.sum(weights)

    return slope, intercept

def deming_bootstrap(
    x,
    y,
    stdx=1.0,
    stdy=1.0,
    n_boot=2000,
    random_state=None,
):
    """
    Bootstrap uncertainty estimates for heteroscedastic
    Deming regression.

    Parameters
    ----------
    x, y : array-like
        Observed x and y values.

    stdx, stdy : float or array-like
        1-sigma uncertainties for each observation.
        Scalars or arrays of length n.

    n_boot : int
        Number of bootstrap replicates.

    random_state : int or None
        Random seed for reproducibility.

    Returns
    -------
    result : dict
        Dictionary containing the original fit, bootstrap
        distributions, and percentile confidence intervals.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    n = len(x)

    if len(y) != n:
        raise ValueError(
            "x and y must have the same length."
        )

    # Expand scalar uncertainties
    stdx = np.asarray(stdx, dtype=float)
    stdy = np.asarray(stdy, dtype=float)

    if stdx.ndim == 0:
        stdx = np.full(n, stdx)

    if stdy.ndim == 0:
        stdy = np.full(n, stdy)

    if len(stdx) != n or len(stdy) != n:
        raise ValueError(
            "stdx and stdy must be scalar or have length n."
        )

    rng = np.random.default_rng(random_state)

    # Fit original data
    slope, intercept = linreg_deming_heteroscedastic(
        x,
        y,
        stdx=stdx,
        stdy=stdy,
    )

    slopes = np.empty(n_boot)
    intercepts = np.empty(n_boot)

    # Bootstrap
    for i in range(n_boot):

        # Resample observations together with their uncertainties
        idx = rng.integers(0, n, size=n)

        x_boot = x[idx]
        y_boot = y[idx]
        stdx_boot = stdx[idx]
        stdy_boot = stdy[idx]

        try:
            slope_boot, intercept_boot = (
                linreg_deming_heteroscedastic(
                    x_boot,
                    y_boot,
                    stdx=stdx_boot,
                    stdy=stdy_boot,
                )
            )

            slopes[i] = slope_boot
            intercepts[i] = intercept_boot

        except RuntimeError:
            slopes[i] = np.nan
            intercepts[i] = np.nan

    # Remove failed bootstrap fits
    valid = (
        np.isfinite(slopes)
        & np.isfinite(intercepts)
    )

    slopes = slopes[valid]
    intercepts = intercepts[valid]

    if len(slopes) == 0:
        raise RuntimeError(
            "All bootstrap fits failed."
        )

    # 95% percentile confidence intervals
    slope_ci = np.percentile(
        slopes,
        [2.5, 97.5],
    )

    intercept_ci = np.percentile(
        intercepts,
        [2.5, 97.5],
    )

    return {
        "slope": slope,
        "intercept": intercept,
        "slope_bootstrap": slopes,
        "intercept_bootstrap": intercepts,
        "slope_ci": slope_ci,
        "intercept_ci": intercept_ci,
        "n_boot": len(slopes),
    }

def linreg_std(x, y, **kwargs):
    slope, intercept, r, p, std_err = stats.linregress(x, y)
    return {"slope": slope, "intercept": intercept}
