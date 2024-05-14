from rail.lephare import *
from rail.core.stage import RailPipeline

def test_lephare_pipeline():
    """Simple test to ensure that a lephare pipeline can be created and saved."""
    lephare = LephareInformer.make_stage(name="lephare_inform")
    pipe = RailPipeline()
    pipe.add_stage(lephare)
    pipe.save('dummy.yml')
