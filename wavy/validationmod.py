"""
    Module to organize the validation procedure
    Consists mostly of functions computing validation metrics
"""
import yaml
import numpy as np
import os

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/region_specs.yaml'))
with open(moddir,'r') as stream:
    region_dict=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/quicklook_specs.yaml'))
if os.path.exists(moddir):
    with open(moddir,'r') as stream:
        quicklook_dict=yaml.safe_load(stream)


# define global functions
def calc_rmsd(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    '''
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    n = len(a1)
    diff2 = (a1-b1)**2
    msd = diff2.sum()/n
    rmsd = np.sqrt(msd)
    return msd, rmsd

def calc_drmsd(a,b):
    '''
    debiased root mean square deviation
    if nans exist the prinziple of marginalization is applied
    '''
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    n = len(a1)
    diff2 = (a1-b1)**2
    msd = diff2.sum()/n
    dmsd = msd - calc_bias(a,b)**2
    rmsd = np.sqrt(msd)
    drmsd = np.sqrt(dmsd)
    return dmsd, drmsd

def calc_scatter_index(obs,model):
    '''
    Scatter index based on rmse and on std of diff
    '''
    msd,rmsd = calc_rmsd(obs,model)
    stddiff = np.nanstd(obs-model)
    SIrmse = rmsd/np.nanmean(obs)*100.
    SIstd = stddiff/np.nanmean(obs)*100.
    return SIrmse,SIstd

def calc_corrcoef(a,b):
    '''
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    '''
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    corr = np.corrcoef(a1,b1)[1,0]
    return corr

def calc_bias(a,b):
    """
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    bias = np.sum(a1-b1)/N
    return bias

def calc_mad(a,b):
    """
    mean absolute deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    mad = np.sum(np.abs(a1-b1))/N
    return mad

def disp_validation(valid_dict):
    '''
    Print to screen validation scores.
    '''
    print('\n')
    print('# ---')
    print('Validation stats')
    print('# ---')
    print('Correlation Coefficient: '
            + '{:0.2f}'.format(valid_dict['corr']))
    print('Mean Absolute Difference: ' + '{:0.2f}'.format(valid_dict['mad']))
    print('Root Mean Squared Difference: '
            + '{:0.2f}'.format(valid_dict['rmsd']))
    print('Debiased Root Mean Squared Difference: '
            + '{:0.2f}'.format(valid_dict['drmsd']))
    print('Bias: ' + '{:0.2f}'.format(valid_dict['bias']))
    print('Scatter Index: ' + '{:0.2f}'.format(valid_dict['SI'][1]))
    print('Mean of Model: ' + '{:0.2f}'.format(valid_dict['mop']))
    print('Mean of Observations: ' + '{:0.2f}'.format(valid_dict['mor']))
    print('Number of Collocated Values: ' + str(valid_dict['nov']))
    print('\n')
    pass


class validation_class():


    def __init__(self,date):
        print ('# ----- ')
        print (" ### Initializing validation_class instance ###")
        print ('# ----- ')

def validate(results_dict,boot=None):
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
    #date_matches = results_dict['datetime']
    if isinstance(results_dict['model_values'],list):
        model_matches = np.array(results_dict['model_values'])
    else: model_matches = results_dict['model_values']
    if isinstance(results_dict['obs_values'],list):
        obs_matches = np.array(results_dict['obs_values'])
    else: obs_matches = results_dict['obs_values']
    if (boot is None or boot ==  False):
        mop = np.nanmean(model_matches)
        mor = np.nanmean(obs_matches)
        msd, rmsd = calc_rmsd(model_matches,obs_matches)
        dmsd, drmsd = calc_drmsd(model_matches,obs_matches)
        nov = len(obs_matches)
        mad = calc_mad(model_matches,obs_matches)
        corr = calc_corrcoef(model_matches,obs_matches)
        bias = calc_bias(model_matches,obs_matches)
        SI = calc_scatter_index(model_matches,obs_matches)
        validation_dict = {
            'mop':mop,
            'mor':mor,
            'msd':msd,
            'nov':nov,
            'rmsd':rmsd,
            'drmsd':drmsd,
            'corr':corr,
            'mad':mad,
            'bias':bias,
            'SI':SI}
    elif boot is True:
        from utils import bootstr, marginalize
        reps=1000
        newmodel,newobs,newidx = marginalize(model_matches,obs_matches)
        obs_boot,boot_idx=bootstr(newobs,reps)
        print (len(obs_boot[np.isnan(obs_boot)]))
        RMSD=np.zeros(reps)*np.nan
        MSD=np.zeros(reps)*np.nan
        BIAS=np.zeros(reps)*np.nan
        CORR=np.zeros(reps)*np.nan
        for i in range(reps):
            results_dict = {
                        #'date_matches':date_matches[newidx[boot_idx[:,i]]],
                        'model_matches':newmodel[boot_idx[:,i]],
                        'sat_matches':newobs[boot_idx[:,i]]}
            try:
                RMSD[i]=validate(results_dict)['rmsd']
                MSD[i]=validate(results_dict)['mad']
                BIAS[i]=validate(results_dict)['bias']
                CORR[i]=validate(results_dict)['corr']
            except IndexError as e:
                print (e)
        validation_dict = {'rmsd':RMSD,'mad':MSD,'bias':BIAS,'corr':CORR}
    return validation_dict
