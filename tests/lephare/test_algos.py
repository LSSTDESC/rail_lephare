import numpy as np
from rail.estimation.algos.lephare import LephareInformer, LephareEstimator
import numpy as np
import lephare as lp
import os
from rail.core.stage import RailStage
import matplotlib.pyplot as plt
import tables_io

DS = RailStage.data_store
DS.__class__.allow_overwrite = True


def test_informer_basic():
    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
    )

    assert inform_lephare.name == "LephareInformer"
    assert inform_lephare.config["name"] == "inform_Lephare"


@pytest.mark.slow
def test_informer_and_estimator(test_data_dir: str):
    trainFile = os.path.join(test_data_dir, "output_table_conv_train.hdf5")
    testFile = os.path.join(test_data_dir, "output_table_conv_test.hdf5")
    traindata_io = tables_io.read(trainFile)
    testdata_io = tables_io.read(testFile)
    lephare_config_file = os.path.join(test_data_dir, "lsst.para")
    lephare_config = lp.read_config(lephare_config_file)
    lp.data_retrieval.get_auxiliary_data(keymap=lephare_config)

    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
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
        np.sum(lephare_estimated.data[0].pdf(np.linspace(0, 3, 300))), 99.66505408115891
    )
