import pytest
from wavy.insitu_module import insitu_class as ic
from wavy.satellite_module import satellite_class as sc
from wavy.model_module import model_class as mc
import xarray as xr
import wavy.triple_collocation as tc
import numpy as np
import pandas as pd


def test_triple_collocation(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {"insitu": ico, "satellite": sco, "model": mco}
    # Apply triple collocation
    ref = "insitu"
    tc_res = tc.triple_collocation(dict_data, ref=ref)

    assert isinstance(tc_res, pd.DataFrame)
    assert tc_res.attrs["ref"] in dict_data.keys()
    assert len(tc_res) == 3
    assert len(list(tc_res)) == 6

    # Simulated data
    n = 1000
    T = [np.sin(0.2 * i) + 1.2 for i in range(n)]
    b_x = 1.1
    b_y = 0.5
    b_z = 1.8

    np.random.seed(1)
    s_x = 0.1
    e_x = np.random.normal(0, s_x, n)

    np.random.seed(5)
    s_y = 0.2
    e_y = np.random.normal(0, s_y, n)

    np.random.seed(11)
    s_z = 0.5
    e_z = np.random.normal(0, s_z, n)

    X = [b_x * T[i] + e_x[i] for i in range(n)]
    Y = [b_y * T[i] + e_y[i] for i in range(n)]
    Z = [b_z * T[i] + e_z[i] for i in range(n)]

    dict_data = {"X": X, "Y": Y, "Z": Z}
    ref = "X"

    tc_res = tc.triple_collocation(dict_data, ref=ref)

    assert isinstance(tc_res, pd.DataFrame)
    assert tc_res.attrs["ref"] in dict_data.keys()
    assert len(tc_res) == 3
    assert len(list(tc_res)) == 6


def test_calibration_triplets_cdf_matching(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {
        "insitu": ico.vars.Hs.values,
        "satellite": sco.vars.Hs.values,
        "model": mco.vars.Hs.values,
    }
    # Apply triple collocation
    ref = "insitu"

    data_cal = tc.calibration_triplets_cdf_matching(dict_data, ref=ref, step=0.04)

    assert list(dict_data.keys()) == list(data_cal.keys())


def test_calibration_triplets_tc(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {
        "insitu": ico.vars.Hs.values,
        "satellite": sco.vars.Hs.values,
        "model": mco.vars.Hs.values,
    }
    # Apply triple collocation
    ref = "insitu"

    data_cal = tc.calibration_triplets_tc(dict_data, ref=ref)

    assert list(dict_data.keys()) == list(data_cal.keys())


def test_least_squares_merging(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {
        "insitu": ico.vars.Hs.values,
        "satellite": sco.vars.Hs.values,
        "model": mco.vars.Hs.values,
    }
    # Apply triple collocation
    ref = "insitu"

    least_squares_merge = tc.least_squares_merging(dict_data)

    assert len(least_squares_merge) == len(dict_data["insitu"])


def test_get_mean_spectra(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))

    spectra = tc.get_mean_spectra(ico.vars, varname="Hs", fs=6, nsample=64)

    assert len(spectra) == 32


def test_integrate_r2(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))

    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))

    ps_ico = tc.get_mean_spectra(ico.vars, varname="Hs", fs=6, nsample=64)
    ps_mod = tc.get_mean_spectra(mco.vars, varname="Hs", fs=6, nsample=64)

    r2 = tc.integrate_r2(ps_ico["spectra"], ps_mod["spectra"], ps_ico["f"])

    assert isinstance(r2, float)
    assert not (np.isnan(r2))


def test_filter_collocation_distance(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {"insitu": ico, "satellite": sco, "model": mco}

    data_filtered = tc.filter_collocation_distance(
        dict_data, dist_max=10, name="satellite"
    )

    assert dict_data.keys() == data_filtered.keys()
    assert len(dict_data["insitu"].vars.time) >= len(data_filtered["insitu"].vars.time)
    assert np.max(data_filtered["satellite"].vars.colloc_dist) <= 10


def test_filter_values(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {
        "insitu": ico.vars.Hs.values,
        "satellite": sco.vars.Hs.values,
        "model": mco.vars.Hs.values,
    }
    # Apply triple collocation
    ref = "insitu"

    data_filtered = tc.filter_values(dict_data, ref_data=ref)

    assert len(data_filtered["insitu"]) <= len(dict_data["insitu"])


def test_filter_dynamic_collocation(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd="2014-01-01", ed="2018-12-31", nID="history_cmems_NRT", name="Norne")
    ico.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_ico.nc"))
    # Import satellite data
    sco = sc(sd="2014-01-01", ed="2018-12-31", nID="CCIv1_L3", name="multi")
    sco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_sco.nc"))
    # Import model data
    mco = mc(sd="2014-01-01", ed="2018-12-31", nID="NORA3_hc_waves")
    mco.vars = xr.open_dataset(str(test_data / "triple_collocation/Norne_mco.nc"))
    # Create dictionary for triple collocation function
    dict_data = {
        "insitu": ico.vars.Hs.values,
        "satellite": sco.vars.Hs.values,
        "model": mco.vars.Hs.values,
    }

    data_filtered = tc.filter_dynamic_collocation(
        dict_data, "insitu", "model", max_rel_diff=0.05
    )

    assert dict_data.keys() == data_filtered.keys()
    assert len(dict_data["insitu"]) >= len(data_filtered["insitu"])
    assert all(
        (data_filtered["insitu"][i] - data_filtered["model"][i])
        / data_filtered["insitu"][i]
        <= 0.05
        for i in range(len(data_filtered["insitu"]))
    )


def test_spatial_variance(test_data):

    varalias = "Hs"  # default
    sd = "2023-7-2 00"
    ed = "2023-7-3 00"
    nID = "MO_Draugen_monthly"
    name = "Draugen"
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
    new = ico.populate(
        path=str(
            test_data / "insitu/monthly/Sulafjorden/AR_TS_MO_A-Sulafjorden_202307.nc"
        )
    )
    new.vars = new.vars.dropna(dim="time")

    variance = tc.spatial_variance(
        new.vars,
        period_meas=60,
        varalias="Hs",
        sd="2023-07-02",
        ed="2023-07-03",
        n_max=2,
        n_min=2,
        time_unit="min",
    )

    assert len(variance) == 1
    assert variance["res"][0] == 2
    assert variance["var"][0] >= 0


def test_merge_variance(test_data):

    list_filenames = ["df_spat_var_sat_s3a.csv", "df_spat_var_sat_s3b.csv"]
    df_var_sat = tc.merge_variance(
        str(test_data / "triple_collocation") + "/",
        list_filenames,
        res_max=None,
        res_factor=1.0,
    )

    assert len(df_var_sat) == 127
    assert df_var_sat["res"].values[0] == 2.0
    assert df_var_sat["res"].values[-1] == 128.0
    assert all(df_var_sat["var"] >= 0)


def test_calculate_r2_spatial_variance(test_data):

    path_files = str(test_data / "triple_collocation/") + "/"
    data_tc = {
        "in-situ": xr.open_dataset(path_files + "Norne_ico.nc").Hs.values,
        "satellite": xr.open_dataset(path_files + "Norne_sco.nc").Hs.values,
        "model": xr.open_dataset(path_files + "Norne_mco.nc").Hs.values,
    }

    list_filenames = ["df_spat_var_sat_s3a.csv", "df_spat_var_sat_s3b.csv"]
    df_var_sat = tc.merge_variance(
        path_files, list_filenames, res_max=None, res_factor=1.0
    )
    list_filenames = ["df_spat_var_mod_s3a.csv", "df_spat_var_mod_s3b.csv"]
    df_var_mod = tc.merge_variance(
        path_files, list_filenames, res_max=None, res_factor=1.0
    )

    df_fit, r2, s_z, cal_cst = tc.calculate_r2_spatial_variance(
        df_var_sat, df_var_mod, data_tc, "in-situ", "satellite", "model", 5
    )

    assert isinstance(df_fit, pd.DataFrame)
    assert isinstance(s_z, float)
    assert isinstance(r2, float)
    assert isinstance(cal_cst, dict)
    assert cal_cst.keys() == data_tc.keys()
    assert df_fit["res"].values[0] == 2.0
    assert round(df_fit["res"].values[-1], 8) == 128.0


def test_power_spectra(test_data):

    varalias = "Hs"  # default
    sd = "2023-7-2 00"
    ed = "2023-7-3 00"
    nID = "MO_Draugen_monthly"
    name = "Draugen"
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
    new = ico.populate(
        path=str(
            test_data / "insitu/monthly/Sulafjorden/AR_TS_MO_A-Sulafjorden_202307.nc"
        )
    )
    new.vars = new.vars.dropna(dim="time")

    df_ps = tc.power_spectra(
        new.vars,
        period_meas=10,
        varalias="Hs",
        time_unit="min",
        fs=6,
        nsample=64,
        sd="2023-07-02",
        ed="2023-07-03",
    )

    assert len(df_ps) == 32
    assert 1 / df_ps["f"].values[-1] == 2 * (1 / 6)
    assert 1 / df_ps["f"].values[0] == 64 * (1 / 6)
    assert all(df_ps["spectra"] >= 0)


def test_merge_spectra(test_data):

    list_filenames = ["df_spectra_sat_s3a.csv", "df_spectra_sat_s3b.csv"]
    df_spec_sat = tc.merge_spectra(
        str(test_data / "triple_collocation") + "/",
        list_filenames,
        res_max=None,
        res_factor=1.0,
    )

    assert len(df_spec_sat) == 32
    assert round(df_spec_sat["res"].values[-1], 8) == 6.3 * 2
    assert round(df_spec_sat["res"].values[0], 8) == 6.3 * 64
    assert all(df_spec_sat["spectra"] >= 0)


def test_calculate_r2_spectra(test_data):

    cal_cst = {"in-situ": 1.0, "satellite": 0.95, "model": 0.91}
    path_files = str(test_data / "triple_collocation") + "/"

    list_filenames = ["df_spectra_sat_s3a.csv", "df_spectra_sat_s3b.csv"]
    df_spec_sat = tc.merge_spectra(
        path_files, list_filenames, res_max=None, res_factor=1.0
    )
    list_filenames = ["df_spectra_mod_s3a.csv", "df_spectra_mod_s3b.csv"]
    df_spec_mod = tc.merge_spectra(
        path_files, list_filenames, res_max=None, res_factor=1.0
    )

    df_r2 = tc.calculate_r2_spectra(
        df_spec_sat,
        df_spec_mod,
        cal_cst,
        "satellite",
        "model",
        np.arange(6.3, 18.9 + 1, 6.3),
    )

    assert len(df_r2) == 3
    assert df_r2["res"].values[0] == 6.3
    assert all(df_r2["r2"] >= 0)


def test_bin_tc(test_data):

    path_files = str(test_data / "triple_collocation") + "/"
    data_tc = {
        "in-situ": xr.open_dataset(path_files + "Norne_ico.nc").Hs.values,
        "satellite": xr.open_dataset(path_files + "Norne_sco.nc").Hs.values,
        "model": xr.open_dataset(path_files + "Norne_mco.nc").Hs.values,
    }

    res_bin = tc.bin_tc(
        data_tc,
        metric="rmse",
        vmin=0.5,
        vmax=4.5,
        step=1.0,
        ref_filter="satellite",
        ref_tc="in-situ",
        transfo_func=None,
        cal=True,
    )

    assert isinstance(res_bin, pd.DataFrame)
    assert len(res_bin) == 4
    assert all(res_bin >= 0)


def test_MARD(test_data):

    path_files = str(test_data / "triple_collocation") + "/"
    data_tc = {
        "in-situ": xr.open_dataset(path_files + "Norne_ico.nc").Hs.values,
        "satellite": xr.open_dataset(path_files + "Norne_sco.nc").Hs.values,
        "model": xr.open_dataset(path_files + "Norne_mco.nc").Hs.values,
    }

    res_bin = tc.bin_tc(
        data_tc,
        metric="rmse",
        vmin=0.5,
        vmax=4.5,
        step=1.0,
        ref_filter="satellite",
        ref_tc="in-situ",
        transfo_func=None,
        cal=True,
    )

    tc_res = tc.triple_collocation(data_tc, ref="in-situ")
    df_MARD = tc.MARD(tc_res, res_bin, metric="rmse")

    assert isinstance(df_MARD, pd.DataFrame)
    assert len(df_MARD) == 1
    assert list(df_MARD.columns) == list(data_tc.keys())
    assert all(df_MARD > 0)
