import importlib
import os
from typing import Any

import lephare as lp
import numpy as np
import qp
from astropy.table import Table
from ceci.config import StageParameter as Param
from ceci.config import StageConfig
from rail.core.common_params import SHARED_PARAMS
from rail.estimation.estimator import CatEstimator, CatInformer

# We start with the COSMOS default and override with LSST specific values.
lsst_default_config = lp.default_cosmos_config.copy()
lsst_default_config.update(
    {
        "VERBOSE": "NO",
        "CAT_IN": "undefined",
        "ERR_SCALE": "0.02",
        "FILTER_CALIB": "0",
        "FILTER_FILE": "filter_lsst",
        "FILTER_LIST": "lsst/total_u.pb,lsst/total_g.pb,lsst/total_r.pb,lsst/total_i.pb,lsst/total_z.pb,lsst/total_y.pb",
        "GAL_LIB": "LSST_GAL_BIN",
        "GAL_LIB_IN": "LSST_GAL_BIN",
        "GAL_LIB_OUT": "LSST_GAL_MAG",
        "GLB_CONTEXT": "0",
        "INP_TYPE": "M",
        "MABS_CONTEXT": "0",
        "MABS_REF": "1",
        "MAG_REF": "2",
        "MAG_ABS_QSO": "-30,-20.5",
        "QSO_LIB": "LSST_QSO_BIN",
        "QSO_LIB_IN": "LSST_QSO_BIN",
        "QSO_LIB_OUT": "LSST_QSO_MAG",
        "STAR_LIB": "LSST_STAR_BIN",
        "STAR_LIB_IN": "LSST_STAR_BIN",
        "STAR_LIB_OUT": "LSST_STAR_MAG",
        "ZPHOTLIB": "LSST_STAR_MAG,LSST_GAL_MAG,LSST_QSO_MAG",
        "Z_STEP": "0.02,0.,6.",
    }
)

lsst_phys_config = lsst_default_config.copy()

lsst_phys_config.update(
    {
        # In order to get the physical parameters you need to use
        # Composite Stellar Population synthesis models. Here Bruzual & Charlot (2003).
        # This can be done only for galaxies.
        "GAL_SED": "$LEPHAREDIR/sed/GAL/BC03_CHAB/BC03COMB_MOD.list",
        "GAL_LIB": "LIB_BC03",
        "SEL_AGE": "$LEPHAREDIR/sed/GAL/BC03_CHAB/AGE_BC03COMB.dat",
        # Library
        "GAL_LIB_IN":"LIB_BC03",
        "GAL_LIB_OUT":"BC03_LSST",
        "MOD_EXTINC": "0,12,0,12",
        "EXTINC_LAW": "SB_calzetti.dat,SMC_prevot.dat",
        "EM_LINES": "PHYS",
        "EM_DISPERSION": "1.",
        # FIT
        # We set the redshift to the spec-z value
        "ZPHOTLIB":"BC03_LSST",
        "ZFIX": "YES",
        "SPEC_OUT": "spec_phys",  # We would like to see the output
    }
)

star_default_config = dict(LIB_ASCII="YES")

gal_default_config = dict(
    LIB_ASCII="YES",
    MOD_EXTINC="18,26,26,33,26,33,26,33",
    EXTINC_LAW="SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat",
    EM_LINES="EMP_UV",
    EM_DISPERSION="0.5,0.75,1.,1.5,2.",
)

qso_default_config = dict(
    LIB_ASCII="YES",
    MOD_EXTINC="0,1000",
    EB_V="0.,0.1,0.2,0.3",
    EXTINC_LAW="SB_calzetti.dat",
)


def _add_sub_config(
    config: dict[str, Any],
    sub_config: dict[str, Any],
    prefix: str,
) -> None:
    """Add all sub-config parameters to the stage config with
    the requested prefix.  This will make the correct parameter types
    and defaults
    """
    for key, val in sub_config.items():
        dtype = type(val)
        default = val
        param = Param(dtype=dtype, default=default)
        config[f"{prefix}{key}"] = param


def _get_sub_config(config: StageConfig, prefix: str) -> dict[str, Any]:
    """Extract all config parameters that start with a
    particular prefix into a dict"""
    out_dict = {
        key[len(prefix) :]: val for key, val in config.items() if key.find(prefix) == 0
    }
    return out_dict


class LephareInformer(CatInformer):
    """Inform stage for LephareEstimator

    This class will set templates and filters required for photoz estimation.
    """

    name = "LephareInformer"
    entrypoint_function = "inform"  # the user-facing science function for this class
    interactive_function = "lephare_informer"
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
    _add_sub_config(config_options, lsst_default_config, "lephare.")
    _add_sub_config(config_options, star_default_config, "star.")
    _add_sub_config(config_options, gal_default_config, "gal.")
    _add_sub_config(config_options, qso_default_config, "qso.")

    def __init__(self, args, **kwargs):
        """Init function, init config stuff (COPIED from rail_bpz)"""
        super().__init__(args, **kwargs)

    def validate(self):
        self.lephare_config = _get_sub_config(self.config, "lephare.")

        # Put something in place to allow for not rerunning the prepare stage
        try:
            self.do_prepare = self.config["do_prepare"]
            if self.do_prepare.__class__ is not bool:
                raise RuntimeError("do_prepare argument must be a bool")
        except KeyError:
            self.do_prepare = True

        # We need to ensure the requested redshift grid is propagated
        self.zmin = self.config["zmin"]
        self.zmax = self.config["zmax"]
        self.nzbins = self.config["nzbins"]
        self.dz = (self.zmax - self.zmin) / (self.nzbins - 1)
        Z_STEP = f"{self.dz},{self.zmin},{self.zmax}"
        print(
            f"rail_lephare is setting the Z_STEP config to {Z_STEP} based on the informer params."
        )
        self.config["lephare.Z_STEP"] = Z_STEP
        # We create a run directory with the informer name
        self.run_dir = _set_run_dir(self.config["name"])

    def _set_config(self, lephare_config):
        """Update the lephare config

        Parameters
        ==========
        lephare_config : `dict`
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

        star_config = _get_sub_config(self.config, "star.")
        gal_config = _get_sub_config(self.config, "gal.")
        qso_config = _get_sub_config(self.config, "qso.")

        # The three main lephare specific inform tasks
        if self.do_prepare:
            lp.prepare(
                self.lephare_config,
                star_config=star_config,
                gal_config=gal_config,
                qso_config=qso_config,
            )
        else:
            print(
                f"do_prepare set to False, using precomputed files in {self.run_dir}."
            )

        # Spectroscopic redshifts
        self.szs = training_data[self.config.redshift_col]

        # Run auto adapt on training sample
        input_table = _rail_to_lephare_input(
            training_data, self.config.bands, self.config.err_bands
        )
        # This will return zeros if AUTO_ADAPT is NO
        offsets = lp.calculate_offsets_from_input(self.lephare_config, input_table)
        # Give principle inform config 'model' to instance.
        self.model = dict(
            lephare_version=lp.__version__,
            lephare_config=self.lephare_config,
            offsets=offsets,
            run_dir=self.run_dir,
            star_config=star_config,
            gal_config=gal_config,
            qso_config=qso_config,
        )
        self.add_data("model", self.model)


class LephareEstimator(CatEstimator):
    """LePhare-base CatEstimator"""

    name = "LephareEstimator"
    entrypoint_function = "estimate"  # the user-facing science function for this class
    interactive_function = "lephare_estimator"
    config_options = CatEstimator.config_options.copy()

    # Add Lephare-specific configuration options here
    config_options.update(
        nondetect_val=SHARED_PARAMS,
        mag_limits=SHARED_PARAMS,
        bands=SHARED_PARAMS,
        ref_band=SHARED_PARAMS,
        err_bands=SHARED_PARAMS,
        redshift_col=SHARED_PARAMS,
        lephare_config_from_model=Param(
            bool,
            True,
            "Load lephare config keymap from model",
        ),
        use_inform_offsets=Param(
            bool,
            True,
            msg="Use the zero point offsets computed in the inform stage.",
        ),
        posterior_output=Param(
            int,
            11,
            msg=(
                "Which posterior distribution to output."
                "MASS: 0"
                "SFR: 1"
                "SSFR: 2"
                "LDUST: 3"
                "LIR: 4"
                "AGE: 5"
                "COL1: 6"
                "COL2: 7"
                "MREF: 8"
                "MIN_ZG: 9"
                "MIN_ZQ: 10"
                "BAY_ZG: 11"
                "BAY_ZQ: 12"
            ),
        ),
        output_keys=Param(
            list,
            ["Z_BEST", "CHI_BEST", "ZQ_BEST", "CHI_QSO", "MOD_STAR", "CHI_STAR"],
            msg=(
                "The output keys to add to ancil. These must be in the "
                "output para file. By default we include the best galaxy "
                "and QSO redshift and best star alongside their respective "
                "chi squared."
            ),
        ),
        run_dir=Param(
            str,
            "None",
            msg=(
                "Override for the LEPHAREWORK directory. If None we load it "
                "from the model which is set during the inform stage. This "
                "is to facilitate manually moving intermediate files."
            ),
        ),
        write_outputs=Param(
            bool,
            False,
            msg="Whether to write the output files.",
        ),
    )
    _add_sub_config(config_options, lsst_default_config, "lephare.")

    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)
        self.lephare_config: dict = {}
        self.zmin: float | None = None
        self.zmax: float | None = None
        self.nzbins: int | None = None

    def open_model(self, **kwargs):
        CatEstimator.open_model(self, **kwargs)
        if self.config["lephare_config_from_model"]:
            self.lephare_config = self.model["lephare_config"]
        else:
            self.lephare_config = _get_sub_config(self.config, "lephare.")

        Z_STEP = self.model["lephare_config"]["Z_STEP"]
        self.lephare_config["Z_STEP"] = Z_STEP
        self.dz = float(Z_STEP.split(",")[0])
        self.zmin = float(Z_STEP.split(",")[1])
        self.zmax = float(Z_STEP.split(",")[2])
        self.nzbins = int((self.zmax - self.zmin) / self.dz)
        if self.config["run_dir"] == "None":
            self.run_dir = self.model["run_dir"]
        else:
            self.run_dir = self.config["run_dir"]
        _update_lephare_env(None, self.run_dir)

    def _process_chunk(self, start, end, data, first):
        """Process an individual chunk of sources using lephare.

        Run the equivalent of zphota and get the PDF for every source.
        """
        # Create the lephare input table
        input_table = _rail_to_lephare_input(
            data, self.config.bands, self.config.err_bands
        )
        # Set the desired offsets estimate config overide lephare config overide inform offsets
        if self.config["use_inform_offsets"] and self.model["offsets"] is not None:
            offsets = self.model["offsets"]
            self.lephare_config["APPLY_SYSSHIFT"] = ",".join([str(o) for o in offsets])

        output, photozlist = lp.process(
            self.lephare_config,
            input_table,
            write_outputs=self.config["write_outputs"],
        )

        ng = data[self.config.bands[0]].shape[0]
        # Unpack the pdfs for galaxies
        pdfs = []
        posterior_output = self.config["posterior_output"]
        for i in range(ng):
            pdf = np.array(photozlist[i].pdfmap[posterior_output].vPDF)
            pdfs.append(pdf)
        zgrid = np.array(photozlist[i].pdfmap[posterior_output].xaxis)
        pdfs = np.array(pdfs)
        self.zgrid = zgrid
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
            if output[c].dtype.kind in ["U", "O"]:
                output[c] = output[c].astype("S")
            ancil[c] = np.array(output[c])
        qp_dstn.set_ancil(ancil)
        self._do_chunk_output(qp_dstn, start, end, first, data=data)


def _rail_to_lephare_input(data, mag_cols, mag_err_cols):
    """Take the rail data input and convert it to that expected by lephare

    Parameters
    ==========
    data : OrderedDict
        The RAIL input data chunk being and ordered dictionary of magnitude arrays.
    mag_cols : list
        The names of the magnitude or flux columns.
    mag_err_cols : list
        The names of the magnitude or flux error columns.

    Returns
    =======
    input_table : astropy.table.Table
        The lephare input


    """
    ng = data[mag_cols[0]].shape[0]
    # Make input catalogue in standard lephare format
    input_table = Table()
    try:
        input_table["id"] = data["id"]
    except KeyError:
        input_table["id"] = np.arange(ng)
    # Add all available magnitudes

    context = np.full(ng, 0)
    for n in np.arange(len(mag_cols)):
        input_table[mag_cols[n]] = data[mag_cols[n]].T
        input_table[mag_err_cols[n]] = data[mag_err_cols[n]].T
        # Shall we allow negative fluxes?
        mask = input_table[mag_cols[n]] > 0
        mask &= ~np.isnan(input_table[mag_cols[n]])
        mask &= input_table[mag_err_cols[n]] > 0
        mask &= ~np.isnan(input_table[mag_err_cols[n]])
        context += mask * 2**n
    # Set context to data value if set or else exclude all negative and nan values
    try:
        input_table["context"] = data["context"]
    except KeyError:
        input_table["context"] = context
    try:
        input_table["zspec"] = data["redshift"]
    except KeyError:
        input_table["zspec"] = np.full(ng, -99.0)
    input_table["string_data"] = " "
    return input_table


def _update_lephare_env(lepharedir, lepharework):
    """Update the environment variables and reset the lephare package.

    We may be using the same Python session to run inform with different
    settings. These produce intermediate files which are distinct and we must
    use different runs.

    This simple function updates the environment variables and reloads lephare
    to ensure they are properly used.
    """
    if lepharedir is not None:
        os.environ["LEPHAREDIR"] = lepharedir
    if lepharework is not None:
        os.environ["LEPHAREWORK"] = lepharework
    importlib.reload(lp)


def _set_run_dir(name=None, full_path=None):
    """Create a named run if it doesn't exist otherwise set it to existing.

    lephare has the functionality to set a timed or named run. In general we
    must ensure that each inform run has a distinct run to ensure that
    intermediate files are not overwritten.

    Parameters
    ==========
    name : str
        The name to set the run. If not set we use a default timestamped run.
    full_path : str
        If set we create a run directory wherever the user sets it.
    """
    if name is None:
        run_directory = lp.dm.create_new_run()
    elif full_path:
        run_directory = full_path
    else:
        run_directory = os.path.realpath(f"{lp.dm.lephare_work_dir}/../{name}")
    _update_lephare_env(lp.LEPHAREDIR, run_directory)
    return run_directory
