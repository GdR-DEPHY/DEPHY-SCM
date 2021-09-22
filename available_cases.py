#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 03 January 2020

@author: Romain Roehrig

Define the list of available case and a dictionary of associated subcases
"""

####################################
#### Some initialization
# List of cases
cases = []
# Dictionnary (case, list of subcases)
subcases = {}

####################################
#### CASES
####################################

####################################
#### Stable boundary-layer cases

# GABLS1 Case
case = 'GABLS1'
cases.append(case)
subcases[case] = ['REF','MESONH']

# GABLS4 Case
case = 'GABLS4'
cases.append(case)
subcases[case] = ['STAGE3','STAGE3-SHORT']

####################################
#### Dry convection cases

# AYOTTE Cases
case = 'AYOTTE'
cases.append(case)
subcases[case] = ['00SC','00WC','03SC','05SC','05WC','24SC']

# IHOP Case
case = 'IHOP'
cases.append(case)
subcases[case] = ['REF']

####################################
#### Shallow convection cases

# RICO Case
case = 'RICO'
cases.append(case)
subcases[case] = ['SHORT','MESONH']

# ARMCU Case
case = 'ARMCU'
cases.append(case)
subcases[case] = ['REF','E3SM','MESONH']

# BOMEX Case
case = 'BOMEX'
cases.append(case)
subcases[case] = ['REF']

# SCMS Case
case = 'SCMS'
cases.append(case)
subcases[case] = ['REF']

# MPACE Case
case = 'MPACE'
cases.append(case)
subcases[case] = ['REF']

####################################
#### Stratocumulus cases

# FIRE case
case = 'FIRE'
cases.append(case)
subcases[case] = ['REF']

# SANDU composite cases
case = 'SANDU'
cases.append(case)
subcases[case] = ['REF','FAST','SLOW']

####################################
#### Deep convection cases

case = 'DYNAMO'
cases.append(case)
subcases[case] = ['NSA3A','NSA3A_D1','NSA3A_MJO1']

####################################

def available(case=None):
    """
    List available cases/subcases.
    """

    if case is None:
        print('-'*30, 'Available cases')
        for cc in sorted(cases):
            tmp = '{0:>10}: '.format(cc)
            for ss in sorted(subcases[cc]):
                tmp += ss + ' '
            print(tmp)
        print('-'*60)
    else:
        print('-'*30, 'Available subcase for case =', case)
        for ss in sorted(subcases[case]):
            print (ss)
        print('-'*60)

if __name__ == "__main__":
    available()
