from rail.estimation.algos.lephare import LephareEstimator

def test_basic_estimator_stage_creation():

    lephare_dict = {}

    estimation_stage = LephareEstimator.make_stage(
        name="lephare_estimation",
        **lephare_dict
    )

    assert estimation_stage.name == 'LephareEstimator'
