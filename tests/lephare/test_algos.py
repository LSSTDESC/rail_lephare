from rail.estimation.algos.lephare import LephareInformer, LephareEstimator
import numpy as np


def test_basic_estimator_stage_creation():
    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
    )

    # We need to add some testing data. From the main RAIL repo?
    # inform_lephare.inform(train_data)

    # estimate_lephare = LephareEstimator.make_stage(
    #     name="test_Lephare",
    #     nondetect_val=np.nan,
    #     model=inform_lephare.get_handle("model"),
    #     hdf5_groupname="",
    #     aliases=dict(input="test_data", output="lephare_estim"),
    # )
