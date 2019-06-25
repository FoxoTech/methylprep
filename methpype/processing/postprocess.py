# Lib
import numpy as np


__all__ = ['calculate_beta_value', 'calculate_m_value']


def calculate_beta_value(methylated_noob, unmethylated_noob, offset=100):
    methylated = max(methylated_noob, 0)
    unmethylated = max(unmethylated_noob, 0)

    total_intensity = methylated + unmethylated + offset
    intensity_ratio = methylated / total_intensity
    return intensity_ratio



def calculate_m_value(methylated_noob, unmethylated_noob, offset=1):
    methylated = methylated_noob + offset
    unmethylated = unmethylated_noob + offset

    intensity_ratio = np.true_divide(methylated, unmethylated)
    return np.log2(intensity_ratio)


def calculate_copy_number(methylated_noob, unmethylated_noob):
    total_intensity = methylated_noob + unmethylated_noob
    copy_number = np.log2(total_intensity)
    return copy_number
