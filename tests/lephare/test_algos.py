import numpy as np
from rail.estimation.algos.lephare import LephareInformer, LephareEstimator
import numpy as np
import lephare as lp
import os
from rail.core.stage import RailStage
from rail.core.data import TableHandle
import matplotlib.pyplot as plt
import tables_io
import pytest

DS = RailStage.data_store
DS.__class__.allow_overwrite = True


def test_informer_basic():
    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
    )

    inform_lephare.validate()
    assert inform_lephare.name == "LephareInformer"
    assert inform_lephare.config["name"] == "inform_Lephare"
    # Check config zgrid updated to stage param defaults:
    assert inform_lephare.config["lephare_config"]["Z_STEP"] == "0.01,0.0,3.0"


def test_informer_and_estimator(test_data_dir: str):
    trainFile = os.path.join(test_data_dir, "output_table_conv_train.hdf5")
    testFile = os.path.join(test_data_dir, "output_table_conv_test.hdf5")
    traindata_io = DS.read_file("training_data", TableHandle, trainFile)
    testdata_io = DS.read_file("test_data", TableHandle, testFile)
    # Load the test params with a sparse redshift grid
    lephare_config_file = os.path.join(test_data_dir, "lsst.para")
    lephare_config = lp.read_config(lephare_config_file)
    lp.data_retrieval.get_auxiliary_data(
        keymap=lephare_config,
        additional_files=["examples/output.para"],
    )

    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
        lephare_config=lp.keymap_to_string_dict(lephare_config),
        # Use a very sparse redshift grid to speed up test:
        zmin=0,
        zmax=5,
        nzbins=6,
    )

    inform_lephare.inform(traindata_io)

    assert os.path.isfile(f"{lp.dm.LEPHAREWORK}/lib_bin/LSST_GAL_BIN.bin")

    estimate_lephare = LephareEstimator.make_stage(
        name="test_Lephare",
        nondetect_val=np.nan,
        model=inform_lephare.get_handle("model"),
        hdf5_groupname="",
        aliases=dict(input="test_data", output="lephare_estim"),
    )

    lephare_estimated = estimate_lephare.estimate(testdata_io)
    assert np.isclose(
        np.sum(lephare_estimated.data[0].pdf(np.linspace(0, 3, 300))), 99.66482986408954
    )
