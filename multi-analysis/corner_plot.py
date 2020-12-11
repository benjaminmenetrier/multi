#!/usr/bin/env python

##############################################################################################
# CornerPlot.py                                                 
#						        
# created by N.Baillot
##############################################################################################

# imported pachages:
import numpy as np
import pickle
import corner
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import corner_plot

##############################################################################################
def make_corner_plot(data_set,labels_set,quant,smoothing):
    """Make a corner plot:
    """
    figure = corner.corner(data_set,
                           #range=ranges*len(chain_pars),
                           title_fmt=".5f",
                           smooth = smoothing,
                           labels=labels_set,
                           show_titles=True,
                           bins=100,
                           plot_datapoints = False,
                           plot_density=False,
                           no_fill_contours=False,
                           use_math_text=True,
                           fill_contours=True,
                           weights = [1.]*len(data_set),
                           #truths=expect_pars,
                           #quantiles=(0.00135,0.02275,0.16, 0.84,0.97725,0.99865),
                           #levels=(1-np.exp(-0.5*1**2),1-np.exp(-0.5*2**2),1-np.exp(-0.5*3**2))
                           quantiles=quant,
                           levels=(1-np.exp(-0.5*1**2),1-np.exp(-0.5*2**2),1-np.exp(-0.5*3**2))
                           )
    return figure
##############################################################################################
