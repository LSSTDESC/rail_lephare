# import qp
from rail.estimation.estimator import CatEstimator, CatInformer



class LephareInformer(CatInformer):
    """Inform stage for LephareEstimator
    
    This class will set templates and filters required for the photoz estimation
    """
    name = "LephareInformer"
    config_options = CatInformer.config_options.copy()
    config_options.update()
    
    
    def run(self):
        """Compute the best fit prior parameters
        """

class LephareEstimator(CatEstimator):
    """LePhare-base CatEstimator
    """

    name = 'LephareEstimator'
    config_options = CatEstimator.config_options.copy()

    # Add Lephare-specific configuration options here
    config_options.update()

    def __init__(self, args, comm=None):
        CatEstimator.__init__(self, args, comm=comm)

    def _process_chunk(self, start, end, data, first):
        """Placeholder for the estimation calculation
        """

        # some calculations take place here
        # pdfs = lephare.estimate(data)
        # qp_dstn = qp.Ensemble(qp.interp, data=dict(xvals=self.zgrid, yvals=pdfs))

        # write the results of estimation for this chunk of data
        # self._do_chunk_output(qp_dstn, start, end, first)
