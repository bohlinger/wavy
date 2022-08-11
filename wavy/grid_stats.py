# Module to organize gridding data

# imports
import numpy as np
import tqdm

def grid_mean(gco,**kwargs):
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
    var_grid = np.full((len(glons),len(glats)), 0.0)
    N = np.full((len(glons),len(glats)), 0)
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

    mx = N>0
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

    val_grid = np.full((len(glons),len(glats)), np.nan)
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


def grid_rmse(gco,**kwargs):
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
    var_grid = np.full((len(glons),len(glats)), 0.0)
    N = np.full((len(glons),len(glats)), 0)
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

    mx = N>0
    var_grid[mx] = np.sqrt(var_grid[mx] / N[mx])
    var_grid[~mx] = np.nan

    return var_grid, lon_grid, lat_grid


def apply_metric(gco=None,**kwargs):
    '''
    dispatch table for various validation metrics
    '''
    metric = kwargs.get('metric')
    print("Computing gridded metric",metric,"...")
    dispatch_reader = {
            'mean':grid_mean,
            'mean_group':grid_mean_group,
            'rmse':grid_rmse,
            }
    var_gridded = dispatch_reader[metric](gco,**kwargs)
    print("Computing gridded",metric," done")
    return var_gridded
