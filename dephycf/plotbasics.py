#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide basic plotting function based on matplotlib

Created on 27 November 2019
@author: Romain Roehrig
"""

import numpy as np
import matplotlib.pyplot as plt
plt.set_loglevel('error')

def plot(x,y,x2=None,y2=None,xlim=None,ylim=None,xlabel=None,ylabel=None,title=None,rep_images=None,name=None,label="",label2="",yunits=None):
    plt.plot(x,y,'k',label=label)
    if not(xlim is None): plt.xlim(xlim)
    if not(ylim is None): plt.ylim(ylim)
    if yunits in ['hPa','Pa']: plt.gca().invert_yaxis()
    if not(xlabel is None): plt.xlabel(xlabel)
    if not(ylabel is None): plt.ylabel(ylabel)
    if not(title is None): plt.title(title)

    if not(x2 is None) and not(y2 is None):
        plt.plot(x2,y2,'r',label=label2)
        plt.legend(loc='best')

    if rep_images is None:
      plt.savefig(name)
    else:
      plt.savefig(rep_images + name)
    plt.close()

def plot2D(x,y,z,xlim=None,ylim=None,xlabel=None,ylabel=None,title=None,rep_images=None,name=None,yunits=None):
    nt,nz = z.shape
    X = np.tile(x,(nz,1))
    Y = np.tile(y,(nt,1))
    plt.contourf(X,np.transpose(Y),np.transpose(z))
    if not(xlim is None): plt.xlim(xlim)
    if not(ylim is None): plt.ylim(ylim)
    if yunits in ['hPa','Pa']: plt.gca().invert_yaxis()
    if not(xlabel is None): plt.xlabel(xlabel)
    if not(ylabel is None): plt.ylabel(ylabel)
    if not(title is None): plt.title(title)
    plt.colorbar()
    if rep_images is None:
      plt.savefig(name)
    else:
      plt.savefig(rep_images + name)
    plt.close()
