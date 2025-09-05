#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic plotting functions using matplotlib.

This module provides basic functions to generate simple
1D and 2D plots and save them as image files.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Reduce matplotlib logging noise
plt.set_loglevel("error")

__all__ = ["plot", "plot2D"]


def plot(
    x,
    y,
    x2=None,
    y2=None,
    xlim=None,
    ylim=None,
    xlabel=None,
    ylabel=None,
    title=None,
    rep_images=None,
    name=None,
    label="",
    label2="",
    yunits=None,
):
    """
    Create and save a 1D line plot.

    Parameters
    ----------
    x : array-like
        X-axis values for the first curve.
    y : array-like
        Y-axis values for the first curve.
    x2 : array-like, optional
        X-axis values for the second curve (if any).
    y2 : array-like, optional
        Y-axis values for the second curve (if any).
    xlim : :class:`tuple`, optional
        Limits for the x-axis, e.g., (xmin, xmax).
    ylim : :class:`tuple`, optional
        Limits for the y-axis, e.g., (ymin, ymax).
    xlabel : :class:`str`, optional
        Label for the x-axis.
    ylabel : :class:`str`, optional
        Label for the y-axis.
    title : :class:`str`, optional
        Title of the plot.
    rep_images : :class:`str`, optional
        Directory where the image should be saved.
        Default is `None`, and in this case, the current working directory is used.
    name : :class:`str`
        File name of the output image (e.g., "plot.png").
    label : :class:`str`, optional
        Legend label for the first curve defined by x and y..
    label2 : :class:`str`, optional
        Legend label for the second curve (if any), as defined by x2 and y2.
    yunits : :class:`str`, optional
        Units for the y-axis. If 'hPa' or 'Pa', the y-axis is inverted.

    Returns
    -------
    None
        The plot is saved to disk as an image file.
    """
    plt.plot(x, y, "k", label=label)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if yunits in ["hPa", "Pa"]:
        plt.gca().invert_yaxis()
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)

    if x2 is not None and y2 is not None:
        plt.plot(x2, y2, "r", label=label2)
        plt.legend(loc="best")

    # Save figure
    filepath = name if rep_images is None else os.path.join(rep_images, name)
    plt.savefig(filepath)
    plt.close()


def plot2D(
    x,
    y,
    z,
    xlim=None,
    ylim=None,
    xlabel=None,
    ylabel=None,
    title=None,
    rep_images=None,
    name=None,
    yunits=None,
):
    """
    Create and save a 2D contour plot.

    Parameters
    ----------
    x : array-like
        Values along the x-axis (e.g., time).
    y : array-like
        Values along the y-axis (e.g., height).
    z : 2D array-like
        Values to plot (shape must be (nt, nz)).
    xlim : :class:`tuple`, optional
        Limits for the x-axis, e.g., (xmin, xmax).
    ylim : :class:`tuple`, optional
        Limits for the y-axis, e.g., (ymin, ymax).
    xlabel : :class:`str`, optional
        Label for the x-axis.
    ylabel : :class:`str`, optional
        Label for the y-axis.
    title : :class:`str`, optional
        Title of the plot.
    rep_images : :class:`str`, optional
        Directory where the image should be saved.
        Default is None, and in this case the current working directory is used.
    name : :class:`str`
        File name of the output image (e.g., "plot2d.png").
    yunits : :class:`str`, optional
        Units for the y-axis. If 'hPa' or 'Pa', the y-axis is inverted.

    Returns
    -------
    None
        The plot is saved to disk as an image file.
    """
    nt, nz = z.shape
    X = np.tile(x, (nz, 1))
    Y = np.tile(y, (nt, 1))

    plt.contourf(X, Y.T, z.T)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if yunits in ["hPa", "Pa"]:
        plt.gca().invert_yaxis()
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)

    plt.colorbar()

    # Save figure
    filepath = name if rep_images is None else os.path.join(rep_images, name)
    plt.savefig(filepath)
    plt.close()
