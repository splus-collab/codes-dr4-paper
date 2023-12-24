# Read data files using astropy.

import os
import warnings
# warnings.filterwarnings("ignore")
from pathlib import Path
import sys
if str(Path.cwd().parent) not in sys.path:
    sys.path.append(str(Path.cwd().parent))

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, join, unique, vstack, Column
import pandas as pd
import seaborn as sns
import matplotlib as mpl

pd.set_option('display.max_columns', None)
from matplotlib import rc
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["figure.figsize"] = (15,11)
plt.rcParams["font.size"] = 18
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.labelsize"] = 13
plt.rcParams["ytick.labelsize"] = 13
mpl.rcParams.update({'font.size': 16, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
mpl.rcParams['axes.linewidth'] = 1.0

colors = {
    'i': (193,5,0),
    'r': (254,95,3),
    'g': (5,195,255),
    'z': (151,0,8),
    'u': (219,3,174),
    'J0378': (195,2,188),
    'J0395': (164,8,199),
    'J0410': (142,2,231),
    'J0430': (61,1,232),
    'J0515': (32,240,17),
    'J0660': (246,1,8),
    'J0861': (178,0,4)
}


# read tables from data folder

psf = Table.read('../data/detections_in_bands_psf.csv')
single = Table.read('../data/detections_in_bands_single.csv')
dual = Table.read('../data/detections_in_bands.csv')

# make 4x3 grid using gridspec

fig = plt.figure(constrained_layout=True, figsize=(16, 9), dpi=300)
gs = fig.add_gridspec(nrows=3, ncols=4, wspace=0.1, hspace=0.1) #, height_ratios=[1], width_ratios=[2,2,2])

# gs.update(left=0.08, right=0.98, bottom=0.10, top=0.98, wspace=0.4, hspace=0.03)

# creating subplots for the 12 bands
ax_list = []
# add loop for plots in gs
for row in range(3):
    for col in range(4):
        ax_list.append(fig.add_subplot(gs[row, col]))

for i,band in enumerate(['u','J0378','J0395','J0410','J0430','g','J0515','r','J0660','i','J0861','z']):
    band2r_psf = []
    band2r_single = []
    band2r_dual = []
    # band2r_field = []
    for limit in range(14,24):
        band2r_psf.append(psf[(psf['limit_range']==limit) & (psf['band']==band)]['count']/psf[(psf['limit_range']==limit) & (psf['band']=='r')]['count'])
        band2r_single.append(single[(single['limit_range']==limit) & (single['band']==band)]['count']/single[(single['limit_range']==limit) & (single['band']=='r')]['count'])
        band2r_dual.append(dual[(dual['limit_range']==limit) & (dual['band']==band)]['count']/dual[(dual['limit_range']==limit) & (dual['band']=='r')]['count'])
        # band2r_field.append(field[(field['limit_range']==limit) & (field['band']==band)]['count']/field[(field['limit_range']==limit) & (field['band']=='r')]['count'])
    
    ax_list[i].plot(range(14,24), band2r_dual, label='dual', color=np.array(colors[band])/255, lw=2, ls='-', marker='o', markersize=6)
    ax_list[i].plot(range(14,24), band2r_single, label="single", color=np.array(colors[band])/255, ls=':', lw=2, marker='^', markersize=6)
    ax_list[i].plot(range(14,24), band2r_psf, label="psf", color=np.array(colors[band])/255, ls='--', lw=2, marker='s', markersize=6)

        # ax_list[i].plot(range(14,24), band2r_field, label=[band+"_dual", band+"_single", band+"_psf"], color=np.array(colors[band])/255, lw=3, ls='-.', marker='s', markersize=8)

    # ax_list[i].set_ylabel('completeness')
    # ax_list[i].set_xlabel('R')
    
    if i % 4 == 0:
        ax_list[i].set_ylabel('completeness', fontsize=22)
    
    if i >= 8:
        ax_list[i].set_xlabel('r (auto)', fontsize=22)
    
    ax_list[i].text(14.5,0.5,band, fontsize=18)
    ax_list[i].legend(loc='lower left', ncol=1, frameon=False, fontsize=14)
    ax_list[i].set_xlim(13.5,23.5)
    ax_list[i].set_ylim(-0.1,1.1)
    ax_list[i].set_xticks(np.arange(14,24,2))
    ax_list[i].set_yticks(np.arange(0,1.1,0.2))
    ax_list[i].grid(ls='--', alpha=0.5)
    ax_list[i].set_axisbelow(True)
    ax_list[i].tick_params(axis='x', labelsize=14)  # Set x-ticks font size
    ax_list[i].tick_params(axis='y', labelsize=14)  # Set y-ticks font size


plt.savefig('../plots/completeness.png', dpi=100)
plt.close()
