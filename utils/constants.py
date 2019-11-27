#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

# Gravity acceleration (m s-2)
g=9.80665

# Boltzman constant (J K-1)
boltzman=1.380658e-23

# Avogadro number (mol-1)
avogadro=6.0221367e+23

# Gaz constant (J mol-1 K-1)
R=avogadro*boltzman

# Dry air molar mass (g mol-1)
Md=28.9644

# Water vapor molar mass (g mol-1)
Mv=18.0153

# Dry air gaz constant (J kg-1 K-1)
Rd=1000.*R/Md

# Water vapor gaz constant (J kg-1 K-1)
Rv=1000.*R/Mv

# Dry air specific capacity at constant pressure (J kg-1 K-1)
Cpd=3.5*Rd

# Water vapor specific capacity at constant pressure (J kg-1 K-1)
Cpv=4.*Rv

# Derived constant
kappa=Rd/Cpd
eps=Rv/Rd-1.0

# Reference pressure
p0 = 100000.

