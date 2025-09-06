#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Physical constants for atmospheric science.

This module defines a collection of widely used physical constants
for thermodynamics and atmospheric physics. Each constant is exposed
as a simple Python variable, documented with its value, units and
scientific meaning.
All constants are of :class:`float` type.
"""

# =============================================================================
# Fundamental constants
# =============================================================================

g = 9.80665
"""Standard gravity acceleration (m s-2).

.. math:: g = 9.80665 \; \\text{m s}^{-2}
"""

boltzman = 1.380658e-23
"""Boltzmann constant (J K-1).

.. math:: k_B = 1.380658 \\times 10^{-23} \; \\text{J K}^{-1}
"""

avogadro = 6.0221367e23
"""Avogadro's number (mol-1).

.. math:: N_A = 6.0221367 \\times 10^{23} \; \\text{mol}^{-1}
"""

R = avogadro * boltzman
"""Universal gas constant (J mol-1 K-1).

.. math:: R = N_A \\times k_B \\approx 8.31451 \; \\text{J mol}^{-1} \\text{ K}^{-1}
"""
sigma = 5.67e-8
"""Stefan–Boltzmann constant (W m-2 K-4).

.. math:: \sigma = 5.67 \\times 10^{-8} \; \\text{W m}^{-2} \\text{ K}^{-4}
"""

# =============================================================================
# Molecular masses
# =============================================================================

Md = 28.9644
"""Dry air molar mass (g mol-1).

.. math:: M_d = 28.9644 \; \\text{g mol}^{-1}
"""

Mv = 18.0153
"""Water vapor molar mass (g mol-1).

.. math:: M_v = 18.0153 \; \\text{g mol}^{-1}
"""

Mo3 = 47.9982
"""Ozone molar mass (g mol-1).

.. math:: M_{O_3} = 47.9982 \; \\text{g mol}^{-1}
"""

# =============================================================================
# Specific gas constants
# =============================================================================

Rd = 1000.0 * R / Md
"""Gas constant for dry air (J kg-1 K-1).

.. math:: R_d = 1000. \\times \\frac{R}{M_d} \\approx 287.1 \; \\text{g mol}^{-1} \\text{ K}^{-1}
"""

Rv = 1000.0 * R / Mv
"""Gas constant for water vapor (J kg-1 K-1).

.. math:: R_v = 1000. \\times \\frac{R}{M_v} \\approx 461.5 \; \\text{g mol}^{-1} \\text{ K}^{-1}
"""

# =============================================================================
# Specific heats at constant pressure
# =============================================================================

Cpd = 3.5 * Rd
"""Dry air specific heat at constant pressure (J kg-1 K-1).

.. math:: C_{pd} = 3.5 \\times R_d \\approx 1004.7 \; \\text{J kg}^{-1} \\text{ K}^{-1}
"""

Cpv = 4.0 * Rv
"""Water vapor specific heat at constant pressure (J kg-1 K-1).

.. math:: C_{pv} = 4 \\times R_v \\approx 1846.1 \; \\text{J kg}^{-1} \\text{ K}^{-1}
"""

# =============================================================================
# Latent heats (at 0 °C)
# =============================================================================

Lv = 2.5008e6
"""Latent heat of water vaporization (J kg-1).

.. math:: L_v = 2.5008 \\times 10^6 \; \\text{J kg}^{-1}
"""

Ls = 2.8345e6
"""Latent heat of water sublimation (J kg-1).

.. math:: L_s = 2.8345 \\times 10^6 \; \\text{J kg}^{-1}
"""

Lf = Ls - Lv
"""Latent heat of water fusion (J kg-1).

.. math:: L_f = L_s - L_v = 0.3337 \\times 10^6 \; \\text{J kg}^{-1}
"""

# =============================================================================
# Derived constants
# =============================================================================

kappa = Rd / Cpd
"""Adiabatic dry air constant (dimensionless).

.. math:: \kappa = \\frac{R_d}{C_{pd}} \\approx 0.286
"""

eps = Rd / Rv
"""Ratio between dry air and water vapor gas constants (dimensionless).

.. math:: \epsilon = \\frac{R_d}{R_v} \\approx 0.622
"""

# =============================================================================
# Reference values
# =============================================================================

p0 = 100000.0
"""Reference pressure (Pa).

.. math:: p_0 = 1000 \; \\text{hPa}
"""


