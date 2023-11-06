# Module to organize gridding data

# imports
import numpy as np
import tqdm
from wavy.validationmod import validate

def grid_mean(gco, **kwargs):
    """
    purpose: computes gridded means
    arguments:
        Midx -> index matrix Midx[0]: indices for grid longitude
                                      indices for grid latitude
        glons -> grid longitude
        glats -> grid latitude
        ovals -> observation values

    returns:
        var_grid -> gridded variable
        lon_grid -> longitude grid
        lat_grid -> latitude grid
    """ 

    # collect needed variables
    if gco is None:
        Midx = kwargs.get('Midx')
        glons = kwargs.get('glons')
        glats = kwargs.get('glats')
        ovals = kwargs.get('ovals')
    else:
        Midx = gco.Midx_clean
        glons = gco.glons
        glats = gco.glats
        ovals = gco.ovals_clean

    # initialize grid
    var_grid = np.full((len(glons), len(glats)), 0.0)
    N = np.full((len(glons), len(glats)), 0)
    lat_grid, lon_grid = np.meshgrid(glats, glons)

    assert len(ovals) == len(Midx.T)

    pbarlen = len(ovals)
    pbar = tqdm.tqdm(total=pbarlen)

    for o, idx in zip(ovals, Midx.T):
        iy = idx[0]
        ix = idx[1]

        var_grid[iy, ix] += o
        N[iy, ix] += 1

        pbar.update(1)

    mx = N > 0
    var_grid[mx] = var_grid[mx] / N[mx]
    var_grid[~mx] = np.nan

    return var_grid, lon_grid, lat_grid

def grid_mean_group(gco, **kwargs):
    """
    purpose: computes gridded means with group strategy
    arguments:
        Midx -> index matrix Midx[0]: indices for grid longitude
                                      indices for grid latitude
        glons -> grid longitude
        glats -> grid latitude
        ovals -> observation values

    returns:
        var_grid -> gridded variable
        lon_grid -> longitude grid
        lat_grid -> latitude grid
    """

    # collect needed variables
    if gco is None:
        Midx = kwargs.get('Midx')
        glons = kwargs.get('glons')
        glats = kwargs.get('glats')
        ovals = kwargs.get('ovals')
    else:
        Midx = gco.Midx_clean
        glons = gco.glons
        glats = gco.glats
        ovals = gco.ovals_clean

    val_grid = np.full((len(glons), len(glats)), np.nan)
    lat_grid, lon_grid = np.meshgrid(glats, glons)

    Midx = Midx.astype(int)

    # Only do this once:
    fidx = np.ravel_multi_index(Midx, val_grid.shape)

    isort = fidx.argsort()
    fidx = fidx[isort]
    vals = ovals[isort]

    isplit = np.unique(fidx, return_index=True)[1]

    gfidx = fidx[isplit]
    gvals = np.split(vals, isplit[1:])

    assert len(gfidx) == len(gvals)

    iyy, ixx = np.unravel_index(gfidx, val_grid.shape)

    pbarlen = len(gfidx)
    pbar = tqdm.tqdm(total=pbarlen)

    for iy, ix, vs in zip(iyy, ixx, gvals):
        val_grid[iy, ix] = np.mean(vs)
        pbar.update(1)

    return val_grid, lon_grid, lat_grid

def grid_stats_group(gco, **kwargs):
    """
    purpose: computes gridded rmse with group strategy
    arguments:
        Midx -> index matrix Midx[0]: indices for grid longitude
                                      indices for grid latitude
        glons -> grid longitude
        glats -> grid latitude
        ovals -> observation values
        mvals -> model values

    returns:
        var_grid -> gridded variable
        lon_grid -> longitude grid
        lat_grid -> latitude grid
    """

    # collect needed variables
    if gco is None:
        Midx = kwargs.get('Midx')
        glons = kwargs.get('glons')
        glats = kwargs.get('glats')
        ovals = kwargs.get('ovals')
        mvals = kwargs.get('mvals')
    else:
        Midx = gco.Midx_clean
        glons = gco.glons
        glats = gco.glats
        ovals = gco.ovals_clean
        mvals = gco.mvals_clean

    mop_grid = np.full((len(glons), len(glats)), np.nan)
    mor_grid = np.full((len(glons), len(glats)), np.nan)
    msd_grid = np.full((len(glons), len(glats)), np.nan)
    rmsd_grid = np.full((len(glons), len(glats)), np.nan)
    nrmsd_grid = np.full((len(glons), len(glats)), np.nan)
    drmsd_grid = np.full((len(glons), len(glats)), np.nan)
    corr_grid = np.full((len(glons), len(glats)), np.nan)
    mad_grid = np.full((len(glons), len(glats)), np.nan)
    bias_grid = np.full((len(glons), len(glats)), np.nan)
    nbias_grid = np.full((len(glons), len(glats)), np.nan)
    SI_grid = np.full((len(glons), len(glats)), np.nan)
    nov_grid = np.full((len(glons), len(glats)), np.nan)
    mar_grid = np.full((len(glons), len(glats)), np.nan)

    lat_grid, lon_grid = np.meshgrid(glats, glons)

    Midx = Midx.astype(int)

    # Only do this once:
    fidx = np.ravel_multi_index(Midx, mop_grid.shape)

    isort = fidx.argsort()
    fidx = fidx[isort]
    ovals_sort = ovals[isort]
    if mvals is None:
        mvals = np.zeros(ovals.shape)*np.nan
    mvals_sort = mvals[isort]

    isplit = np.unique(fidx, return_index=True)[1]

    gfidx = fidx[isplit]
    g_ovals_sort = np.split(ovals_sort, isplit[1:])
    g_mvals_sort = np.split(mvals_sort, isplit[1:])

    assert len(gfidx) == len(g_ovals_sort)
    assert len(gfidx) == len(g_mvals_sort)

    iyy, ixx = np.unravel_index(gfidx, mop_grid.shape)

    pbarlen = len(gfidx)
    pbar = tqdm.tqdm(total=pbarlen)

    for iy, ix, ov, mv in zip(iyy, ixx, g_ovals_sort, g_mvals_sort):
        rdict = {
                 'model_values': mv,
                 'obs_values': ov
                 }
        validation_dict = validate(rdict)
        mop_grid[iy, ix] = validation_dict['mop']
        mor_grid[iy, ix] = validation_dict['mor']
        mad_grid[iy, ix] = validation_dict['mad']
        msd_grid[iy, ix] = validation_dict['msd']
        rmsd_grid[iy, ix] = validation_dict['rmsd']
        nrmsd_grid[iy, ix] = validation_dict['nrmsd']
        drmsd_grid[iy, ix] = validation_dict['drmsd']
        corr_grid[iy, ix] = validation_dict['corr']
        bias_grid[iy, ix] = validation_dict['bias']
        nbias_grid[iy, ix] = validation_dict['nbias']
        SI_grid[iy, ix] = validation_dict['SI'][1]
        mar_grid[iy, ix] = validation_dict['mar']
        nov_grid[iy, ix] = validation_dict['nov']

        pbar.update(1)

    val_grid_dict = {
            'mop': mop_grid,
            'mor': mor_grid,
            'mad': mad_grid,
            'msd': msd_grid,
            'rmsd': rmsd_grid,
            'nrmsd': nrmsd_grid,
            'drmsd': drmsd_grid,
            'corr': corr_grid,
            'bias': bias_grid,
            'nbias': nbias_grid,
            'SI': SI_grid,
            'mar': mar_grid,
            'nov': nov_grid
            }

    return val_grid_dict, lon_grid, lat_grid


def grid_rmse(gco, **kwargs):
    """
    purpose: computes gridded rmse
    arguments:
        Midx -> index matrix Midx[0]: indices for grid longitude
                                      indices for grid latitude
        glons -> grid longitude
        glats -> grid latitude
        ovals -> observation values
        mvals -> model values

    returns:
        var_grid -> gridded variable
        lon_grid -> longitude grid
        lat_grid -> latitude grid
    """

    # collect needed variables
    if gco is None:
        Midx = kwargs.get('Midx')
        glons = kwargs.get('glons')
        glats = kwargs.get('glats')
        ovals = kwargs.get('ovals')
        mvals = kwargs.get('mvals')
    else:
        Midx = gco.Midx_clean
        glons = gco.glons
        glats = gco.glats
        ovals = gco.ovals_clean
        mvals = gco.mvals_clean

    # initialize grid
    var_grid = np.full((len(glons), len(glats)), 0.0)
    N = np.full((len(glons), len(glats)), 0)
    lat_grid, lon_grid = np.meshgrid(glats, glons)

    assert len(ovals) == len(Midx.T)
    assert len(mvals) == len(Midx.T)

    pbarlen = len(ovals)
    pbar = tqdm.tqdm(total=pbarlen)

    for o, m, idx in zip(ovals, mvals, Midx.T):
        iy = idx[0]
        ix = idx[1]

        var_grid[iy, ix] += (o-m)**2
        N[iy, ix] += 1

        pbar.update(1)

    mx = N > 0
    var_grid[mx] = np.sqrt(var_grid[mx] / N[mx])
    var_grid[~mx] = np.nan

    return var_grid, lon_grid, lat_grid


def apply_metric(gco=None, **kwargs):
    '''
    dispatch table for various validation metrics
    '''
    metric = kwargs.get('metric', 'all')
    print("Computing gridded metric:", metric, "...")
    dispatch_reader = {
            'mean': grid_mean,
            'mean_group': grid_mean_group,
            'rmse': grid_rmse,
            'all': grid_stats_group,
            }
    var_gridded = dispatch_reader[metric](gco, **kwargs)
    print("Computing gridded metric:", metric, " -> done")
    return var_gridded
