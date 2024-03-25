# import qp
from rail.estimation.estimator import CatEstimator, CatInformer
from rail.core.common_params import SHARED_PARAMS
import os
import lephare as lp


class LephareInformer(CatInformer):
    """Inform stage for LephareEstimator

    This class will set templates and filters required for photoz estimation
    """

    name = "LephareInformer"
    config_options = CatInformer.config_options.copy()
    config_options.update(
        zmin=SHARED_PARAMS,
        zmax=SHARED_PARAMS,
        nzbins=SHARED_PARAMS,
        nondetect_val=SHARED_PARAMS,
        mag_limits=SHARED_PARAMS,
        bands=SHARED_PARAMS,
        err_bands=SHARED_PARAMS,
        ref_band=SHARED_PARAMS,
        redshift_col=SHARED_PARAMS,
    )

    def __init__(self, args, comm=None):
        """Init function, init config stuff (COPIED from rail_bpz)"""
        CatInformer.__init__(self, args, comm=comm)
        # Default local parameters
        self.config_file = "lsst.para"
        print("config file: " + self.config_file)
        self.lephare_config = lp.read_config(self.config_file)

    def _set_config(self, lephare_config):
        """Update the lephare config

        Parameters
        ==========
        lepahre_config : `dict`
            A dictionary of the lephare config keywords.
        """
        self.lephare_config = lephare_config

    def _create_filter_library(self):
        """Make the filter library files in lephare format"""
        # load filters from config file
        filterLib = lp.FilterSvc.from_config(self.config_file)
        # Get location to store filter files
        filter_output = os.path.join(
            os.environ["LEPHAREWORK"], "filt", self.lephare_config["FILTER_FILE"].value
        )
        # Write filter files
        lp.write_output_filter(
            filter_output + ".dat", filter_output + ".doc", filterLib
        )

    def _create_sed_library(self):
        """Make the SED binary library files in lephare format.

        We separately create the star, quasar and galaxy libraries.
        """
        sedlib = lp.Sedtolib(config_keymap=self.lephare_config)
        sedlib.run(typ="STAR", star_sed="$LEPHAREDIR/examples/STAR_MOD_ALL.list")
        sedlib.run(typ="QSO", qso_sed="$LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list")
        sedlib.run(typ="GAL", gal_sed="$LEPHAREDIR/examples/COSMOS_MOD.list")

    def _create_mag_library(self):
        """Make the magnitudes library file in lephare format.

        We separately create the star, quasar and galaxy libraries.

        TODO: replace hardcoded config options with class config options.
        """
        maglib = lp.MagGal(config_keymap=self.lephare_config)
        maglib.run(typ="STAR", lib_ascii="YES")
        maglib.run(
            typ="QSO",
            lib_ascii="YES",
            mod_extinc="0,1000",
            eb_v="0.,0.1,0.2,0.3",
            extinc_law="SB_calzetti.dat",
        )
        maglib.run(
            typ="GAL",
            lib_ascii="YES",
            mod_extinc="18,26,26,33,26,33,26,33",
            extinc_law=(
                "SMC_prevot.dat,SB_calzetti.dat,"
                + "SB_calzetti_bump1.dat,SB_calzetti_bump2.dat"
            ),
            em_lines="EMP_UV",
            em_dispersion="0.5,0.75,1.,1.5,2.",
        )

    def run(self):
        """Run rail_lephare inform stage.

        This is the basic informer which takes the config and templates and
        makes the inputs required for the run.

        In addition to the three lephare stages making the filter, sed, and
        magnitude libraries we also do some tasks required by all rail inform
        stages.
        """
        # Set training data required for all informers?
        if self.config.hdf5_groupname:
            training_data = self.get_data("input")[self.config.hdf5_groupname]
        else:  # pragma: no cover
            training_data = self.get_data("input")

        # Get number of sources
        ngal = len(training_data[self.config.ref_band])

        # The three main lephare specific inform tasks
        self._create_filter_library()
        self._create_sed_library()
        self._create_mag_library()

        # Spectroscopic redshifts
        self.szs = training_data[self.config.redshift_col]

        # Give principle inform config 'model' to instance.
        self.model = dict(config_file=self.config_file)
        self.add_data("model", self.model)


class LephareEstimator(CatEstimator):
    """LePhare-base CatEstimator"""

    name = "LephareEstimator"
    config_options = CatEstimator.config_options.copy()

    # Add Lephare-specific configuration options here
    config_options.update(
        zmin=SHARED_PARAMS,
        zmax=SHARED_PARAMS,
        nzbins=SHARED_PARAMS,
        nondetect_val=SHARED_PARAMS,
        mag_limits=SHARED_PARAMS,
        bands=SHARED_PARAMS,
        ref_band=SHARED_PARAMS,
        err_bands=SHARED_PARAMS,
        redshift_col=SHARED_PARAMS,
    )

    def __init__(self, args, comm=None):
        CatEstimator.__init__(self, args, comm=comm)

    def _process_chunk(self, start, end, data, first):
        """Placeholder for the estimation calculation"""

        # some calculations take place here
        # pdfs = lephare.estimate(data)
        # qp_dstn = qp.Ensemble(qp.interp, data=dict(xvals=self.zgrid, yvals=pdfs))

        # write the results of estimation for this chunk of data
        # self._do_chunk_output(qp_dstn, start, end, first)
