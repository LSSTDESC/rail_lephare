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
        lephare_config_file=Param(
            str,
            "{}/{}".format(os.path.dirname(os.path.abspath(__file__)), "lsst.para"),
            msg="Path to the lephare config in .para format",
        ),
        star_sed=Param(
            str,
            "$LEPHAREDIR/examples/STAR_MOD_ALL.list",
            msg="Path to text file containing list of star SED templates",
        ),
        qso_sed=Param(
            str,
            "$LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list",
            msg="Path to text file containing list of galaxy SED templates",
        ),
        gal_sed=Param(
            str,
            "$LEPHAREDIR/examples/COSMOS_MOD.list",
            msg="Path to text file containing list of quasar SED templates",
        ),
        star_mag_dict=Param(
            dict,
            dict(
                lib_ascii="YES",
            ),
            msg="Dictionary of values sent to MagGal for stars",
        ),
        gal_mag_dict=Param(
            dict,
            dict(
                lib_ascii="YES",
                mod_extinc="18,26,26,33,26,33,26,33",
                extinc_law=(
                    "SMC_prevot.dat,SB_calzetti.dat,"
                    "SB_calzetti_bump1.dat,SB_calzetti_bump2.dat"
                ),
                em_lines="EMP_UV",
                em_dispersion="0.5,0.75,1.,1.5,2.",
            ),
            msg="Dictionary of values sent to MagGal for galaxies",
        ),
        qso_mag_dict=Param(
            dict,
            dict(
                lib_ascii="YES",
                mod_extinc="0,1000",
                eb_v="0.,0.1,0.2,0.3",
                extinc_law="SB_calzetti.dat",
            ),
            msg="Dictionary of values sent to MagGal for quasars",
        ),
    )

    def __init__(self, args, comm=None):
        """Init function, init config stuff (COPIED from rail_bpz)"""
        CatInformer.__init__(self, args, comm=comm)
        self.lephare_config = lp.read_config(self.config["lephare_config_file"])

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
        filterLib = lp.FilterSvc.from_config(self.config["lephare_config_file"])
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
        sedlib.run(typ="STAR", star_sed=self.config["star_sed"])
        sedlib.run(typ="GAL", gal_sed=self.config["gal_sed"])
        sedlib.run(typ="QSO", qso_sed=self.config["qso_sed"])

    def _create_mag_library(self):
        """Make the magnitudes library file in lephare format.

        We separately create the star, quasar and galaxy libraries.

        TODO: replace hardcoded config options with class config options.
        """
        maglib = lp.MagGal(config_keymap=self.lephare_config)
        maglib.run(typ="STAR", **self.config["star_mag_dict"])
        maglib.run(typ="GAL", **self.config["gal_mag_dict"])
        maglib.run(typ="QSO", **self.config["qso_mag_dict"])

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
        self._create_filter_library()
        self._create_sed_library()
        self._create_mag_library()

        # Spectroscopic redshifts
        self.szs = training_data[self.config.redshift_col]

        # Give principle inform config 'model' to instance.
        self.model = dict(lephare_config_file=self.config["lephare_config_file"])
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

        # Default local parameters
        self.config_file = "{}/{}".format(
            os.path.dirname(os.path.abspath(__file__)), "lsst.para"
        )
        self.lephare_config = lp.read_config(self.config_file)
        self.photz = lp.PhotoZ(self.lephare_config)
        print("init")

    # def open_model(self, **kwargs):
    #     CatEstimator.open_model(self, **kwargs)
    #     self.modeldict = self.model

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
