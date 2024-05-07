from rail.estimation.estimator import CatEstimator, CatInformer
from rail.core.common_params import SHARED_PARAMS
from ceci.config import StageParameter as Param
import os
import lephare as lp
import numpy as np
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
        self.lephare_config = lp.read_config(self.config["lephare_config"])

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

        lp.data_retrieval.get_auxiliary_data(keymap=self.lephare_config)

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

        # Give principle inform config 'model' to instance.
        self.model = dict(lephare_config=self.config["lephare_config"])
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
        lephare_config_file=Param(
            str,
            "{}/{}".format(os.path.dirname(os.path.abspath(__file__)), "lsst.para"),
            msg="Path to the lephare config in .para format",
        ),
    )

    def __init__(self, args, comm=None):
        CatEstimator.__init__(self, args, comm=comm)
        self.lephare_config = lp.read_config(self.config["lephare_config_file"])
        self.photz = lp.PhotoZ(self.lephare_config)

    def _estimate_pdf(self, onesource):
        """Return the pdf of a single source.

        Do we want to resample on RAIL z grid?
        """
        # Check this is the best way to access pdf
        pdf = onesource.pdfmap[11]  # 11 = Bayesian galaxy redshift
        # return the PDF as an array alongside lephare native zgrid
        return np.array(pdf.vPDF), np.array(pdf.xaxis)

    def _process_chunk(self, start, end, data, first):
        """Process an individual chunk of sources using lephare

        Run the equivalent of zphota and get the PDF for every source.
        """
        # write the results of estimation for this chunk of data
        self.zgrid = np.linspace(self.config.zmin, self.config.zmax, self.config.nzbins)

        nz = len(self.zgrid)
        ng = data["redshift"].shape[0]
        flux = np.array([data["mag_{}_lsst".format(b)] for b in "ugrizy"]).T
        flux_err = np.array([data["mag_err_{}_lsst".format(b)] for b in "ugrizy"]).T
        zspec = data["redshift"]

        pdfs = []  # np.zeros((ng, nz))
        zmode = np.zeros(ng)
        zmean = np.zeros(ng)
        zgrid = self.zgrid

        # Loop over all ng galaxies!
        srclist = []
        for i in range(ng):
            oneObj = lp.onesource(i, self.photz.gridz)
            oneObj.readsource(str(i), flux[i], flux_err[i], 63, zspec[i], " ")
            self.photz.prep_data(oneObj)
            srclist.append(oneObj)

        # Run autoadapt to improve zero points
        a0, a1 = self.photz.run_autoadapt(srclist)
        offsets = ",".join(np.array(a0).astype(str))
        offsets = "# Offsets from auto-adapt: " + offsets + "\n"
        print(offsets)

        photozlist = []
        for i in range(ng):
            oneObj = lp.onesource(i, self.photz.gridz)
            oneObj.readsource(str(i), flux[i], flux_err[i], 63, zspec[i], " ")
            self.photz.prep_data(oneObj)
            photozlist.append(oneObj)

        self.photz.run_photoz(photozlist, a0, a1)

        for i in range(ng):
            pdf, zgrid = self._estimate_pdf(photozlist[i])
            pdfs.append(pdf)
            # Take median in case multiple probability densities are equal
            zmode[i] = np.median(
                zgrid[np.where(pdfs[i] == np.max(pdfs[i]))[0].astype(int)]
            )
            zmean[i] = (zgrid * pdfs[i]).sum() / pdfs[i].sum()

        qp_dstn = qp.Ensemble(qp.interp, data=dict(xvals=zgrid, yvals=np.array(pdfs)))
        qp_dstn.set_ancil(dict(zmode=zmode, zmean=zmean))
        self._do_chunk_output(qp_dstn, start, end, first)
