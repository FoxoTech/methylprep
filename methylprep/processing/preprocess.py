# Normal-exponential using out-of-band probes
# normex: negative control probes
# noob: ‘out-of-band’ Infinium I probes

# Lib
import logging
import numpy as np
from statsmodels import robust
from scipy.stats import norm
# App
from ..models import ControlType


__all__ = ['preprocess_noob']


LOGGER = logging.getLogger(__name__)


class BackgroundCorrectionParams():
    __slots__ = (
        'bg_mean',
        'bg_mad',
        'mean_signal',
        'offset',
    )

    def __init__(self, bg_mean, bg_mad, mean_signal, offset=15):
        self.bg_mean = bg_mean
        self.bg_mad = bg_mad
        self.mean_signal = mean_signal
        self.offset = offset


def preprocess_noob(data_container):
    """ the main preprocessing function. Applies background-subtraction and
    NOOB. Sets data_container.methylated and unmethylated values for sample."""
    #LOGGER.info('NOOB: %s', data_container.sample)

    bg_correct_green, params_green = normexp_bg_corrected(data_container.fg_green, data_container.oob_green)
    bg_correct_red, params_red = normexp_bg_corrected(data_container.fg_red, data_container.oob_red)

    data_container.methylated.set_bg_corrected(bg_correct_green, bg_correct_red)
    data_container.unmethylated.set_bg_corrected(bg_correct_green, bg_correct_red)

    ctrl_green = normexp_bg_correct_control(data_container.ctrl_green, params_green)
    ctrl_red = normexp_bg_correct_control(data_container.ctrl_red, params_red)

    mask_green = ctrl_green['Control_Type'].isin(ControlType.normalization_green())
    mask_red = ctrl_red['Control_Type'].isin(ControlType.normalization_red())

    avg_green = ctrl_green[mask_green]['bg_corrected'].mean()
    avg_red = ctrl_red[mask_red]['bg_corrected'].mean()

    rg_ratios = avg_red / avg_green

    red_factor = 1 / rg_ratios

    data_container.methylated.set_noob(red_factor)
    data_container.unmethylated.set_noob(red_factor)


def normexp_bg_corrected(fg_probes, ctrl_probes):
    fg_means = fg_probes['mean_value']
    fg_mean, _fg_mad = huber(fg_means)
    bg_mean, bg_mad = huber(ctrl_probes['mean_value'])
    mean_signal = np.maximum(fg_mean - bg_mean, 10)

    params = BackgroundCorrectionParams(bg_mean, bg_mad, mean_signal)

    corrected_signals = apply_bg_correction(fg_means, params)
    fg_probes['bg_corrected'] = corrected_signals
    return fg_probes, params


def normexp_bg_correct_control(control_probes, params):
    """Function for getting xcs controls for preprocessNoob"""
    control_means = control_probes['mean_value']
    corrected_signals = apply_bg_correction(control_means, params)
    control_probes['bg_corrected'] = corrected_signals
    return control_probes


def apply_bg_correction(mean_values, params):
    """ this function won't work with float16 in practice. limits use to float32 """
    if not isinstance(params, BackgroundCorrectionParams):
        raise ValueError('params is not a BackgroundCorrectionParams instance')

    bg_mean = params.bg_mean
    bg_mad = params.bg_mad
    mean_signal = params.mean_signal
    offset = params.offset

    mu_sf = mean_values - bg_mean - (bg_mad ** 2) / mean_signal

    signal = mu_sf + (bg_mad ** 2) * \
        np.exp(norm(mu_sf, bg_mad).logpdf(0) - norm(mu_sf, bg_mad).logsf(0))

    signal = np.maximum(signal, 1e-6)
    true_signal = signal + offset

    return true_signal


def huber(vector):
    """Huber function. Designed to mirror MASS huber function in R

    Parameters
    ----------
    vector: list
        list of float values

    Returns
    -------
    local_median: float
        calculated mu value
    mad_scale: float
        calculated s value
    """
    num_values = len(vector)
    positive_factor = 1.5
    convergence_tol = 1.0e-6
    mad_scale = robust.mad(vector)
    local_median = np.median(vector)
    init_local_median = local_median

    if not (local_median or mad_scale):
        return local_median, mad_scale

    while True:
        yy = np.minimum(
            np.maximum(
                local_median - positive_factor * mad_scale,
                vector,
            ),
            local_median + positive_factor * mad_scale,
        )

        init_local_median = sum(yy) / num_values

        if abs(local_median - init_local_median) < convergence_tol * mad_scale:
            return local_median, mad_scale

        local_median = init_local_median
