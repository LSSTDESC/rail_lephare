import numpy as np
from rail.estimation.algos.lephare import LephareInformer


def test_basic_informer_stage_creation():
    inform_lephare = LephareInformer.make_stage(
        name="inform_Lephare",
        nondetect_val=np.nan,
        model="lephare.pkl",
        hdf5_groupname="",
    )

    assert inform_lephare.name == "LephareInformer"

    # We need to add some testing data. From the main RAIL repo?
    # inform_lephare.inform(train_data)

# def test_basic_estimator_stage_creation():
#     estimate_lephare = LephareEstimator.make_stage(
#         name="test_Lephare",
#         nondetect_val=np.nan,
#         model=inform_lephare.get_handle("model"),
#         hdf5_groupname="",
#         aliases=dict(input="test_data", output="lephare_estim"),
#     )

#     assert estimate_lephare.name == "LephareEstimator"
