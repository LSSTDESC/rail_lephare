from rail.estimation.estimator import CatEstimator, CatInformer
from rail.core.common_params import SHARED_PARAMS
from ceci.config import StageParameter as Param
import os
import lephare as lp
import numpy as np
from astropy.table import Table
import qp


class LephareInformer(CatInformer):
    """Inform stage for LephareEstimator

    This class will set templates and filters required for photoz estimation.
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
        lephare_config=Param(
            dict,
            lp.read_config(
                "{}/{}".format(os.path.dirname(os.path.abspath(__file__)), "lsst.para")
            ),
            msg="The lephare config keymap.",
        ),
        star_config=Param(
            dict,
            dict(LIB_ASCII=lp.keyword("LIB_ASCII", "YES")),
            msg="Star config overrides.",
        ),
        gal_config=Param(
            dict,
            dict(
                LIB_ASCII=lp.keyword("LIB_ASCII", "YES"),
                MOD_EXTINC=lp.keyword("MOD_EXTINC", "18,26,26,33,26,33,26,33"),
                EXTINC_LAW=lp.keyword(
                    "EXTINC_LAW",
                    "SMC_prevot.dat,SB_calzetti.dat,"
                    "SB_calzetti_bump1.dat,SB_calzetti_bump2.dat",
                ),
                EM_LINES=lp.keyword("EM_LINES", "EMP_UV"),
                EM_DISPERSION=lp.keyword("EM_DISPERSION", "0.5,0.75,1.,1.5,2."),
            ),
            msg="Galaxy config overrides.",
        ),
        qso_config=Param(
            dict,
            dict(
                LIB_ASCII=lp.keyword("LIB_ASCII", "YES"),
                MOD_EXTINC=lp.keyword("MOD_EXTINC", "0,1000"),
                EB_V=lp.keyword("EB_V", "0.,0.1,0.2,0.3"),
                EXTINC_LAW=lp.keyword("EXTINC_LAW", "SB_calzetti.dat"),
            ),
            msg="QSO config overrides.",
        ),
    )

    def __init__(self, args, comm=None):
        """Init function, init config stuff (COPIED from rail_bpz)"""
        CatInformer.__init__(self, args, comm=comm)
        self.lephare_config = self.config["lephare_config"]

    def _set_config(self, lephare_config):
        """Update the lephare config

        Parameters
        ==========
        lepahre_config : `dict`
            A dictionary of the lephare config keywords.
        """
        self.lephare_config = lephare_config

    def run(self):
        """Run rail_lephare inform stage.

        This informer takes the config and templates and makes the inputs
        required for the run.

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
        lp.prepare(
            self.lephare_config,
            star_config=self.config["star_config"],
            gal_config=self.config["gal_config"],
            qso_config=self.config["qso_config"],
        )

        # Spectroscopic redshifts
        self.szs = training_data[self.config.redshift_col]

        # Run auto adapt on training sample

        # We must make a string dictionary to allow pickling and saving
        config_text_dict = dict()
        for k in self.config["lephare_config"]:
            config_text_dict[k] = self.config["lephare_config"][k].value
        # Give principle inform config 'model' to instance.
        self.model = dict(lephare_config=config_text_dict)
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
        lephare_config=Param(
            dict,
            lp.read_config(
                "{}/{}".format(os.path.dirname(os.path.abspath(__file__)), "lsst.para")
            ),
            msg="The lephare config keymap.",
        ),
        output_keys=Param(
            list,
            ["Z_BEST", "ZQ_BEST", "MOD_STAR"],
            msg="The output keys to add to ancil. These must be in the output para file.",
        ),
        offsets=Param(
            list,
            [None, None],
            msg=(
                "The offsets to apply to photometry. If set to None "
                "autoadapt will be run if that key is set in the config."
            ),
        ),
    )

    def __init__(self, args, comm=None):
        CatEstimator.__init__(self, args, comm=comm)
        self.lephare_config = self.config["lephare_config"]
        self.photz = lp.PhotoZ(self.lephare_config)

    def _estimate_pdf(self, onesource):
        """Return the pdf of a single source.

        Do we want to resample on RAIL z grid?
        """
        # Check this is the best way to access pdf
        pdf = onesource.pdfmap[11]  # 11 = Bayesian galaxy redshift
        # return the PDF as an array alongside lephare native zgrid
        return np.array(pdf.vPDF), np.array(pdf.xaxis)

    def _rail_to_lephare_input(self, data):
        """Take the rail data input and convert it to that expected by lephare

        Parameters
        ==========
        data : pandas
            The RAIL input data chunk

        Returns
        =======
        input : astropy.table.Table
            The lephare input


        """
        ng = data["redshift"].shape[0]
        # Make input catalogue in standard lephare format
        input = Table()
        try:
            input["id"] = data["id"]
        except KeyError:
            input["id"] = np.arange(ng)
        # Add all available magnitudes
        for k in data.keys():
            if k.startswith("mag_err"):
                input[k.replace("mag_err", "mag")] = data[k.replace("mag_err", "mag")].T
                input[k] = data[k].T
        # Set context to use all bands
        input["context"] = np.sum([2**n for n in np.arange(ng)])
        input["zspec"] = data["redshift"]
        input["string_data"] = " "
        return input

    def _process_chunk(self, start, end, data, first):
        """Process an individual chunk of sources using lephare

        Run the equivalent of zphota and get the PDF for every source.
        """
        # write the results of estimation for this chunk of data
        self.zgrid = np.linspace(self.config.zmin, self.config.zmax, self.config.nzbins)
        input = self._rail_to_lephare_input(data)
        if self.config["offsets"][0] is None:
            offsets = None
        else:
            offsets = self.config["offsets"]
        output, pdfs, zgrid = lp.process(self.lephare_config, input, offsets=offsets)
        self.zgrid = zgrid

        ng = data["redshift"].shape[0]
        zmode = np.zeros(ng)
        zmean = np.zeros(ng)

        for i in range(ng):
            # Take median in case multiple probability densities are equal
            zmode[i] = np.median(
                zgrid[np.where(pdfs[i] == np.max(pdfs[i]))[0].astype(int)]
            )
            zmean[i] = (zgrid * pdfs[i]).sum() / pdfs[i].sum()

        qp_dstn = qp.Ensemble(qp.interp, data=dict(xvals=zgrid, yvals=np.array(pdfs)))
        ancil = dict(zmode=zmode, zmean=zmean)
        # Add the requested outputs.
        for c in self.config["output_keys"]:
            ancil[c] = np.array(output[c])
        qp_dstn.set_ancil(ancil)
        self._do_chunk_output(qp_dstn, start, end, first)
