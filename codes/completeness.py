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
plt.rcParams["figure.figsize"] = (13,11)
plt.rcParams["font.size"] = 18
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.labelsize"] = 13
plt.rcParams["ytick.labelsize"] = 13
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
mpl.rcParams['axes.linewidth'] = 2.0

colors = {
    'i': (97, 0, 0),
    'r': (255, 99, 0),
    'g': (0, 195, 255),
    'z': (97, 0, 0),
    'u': (97, 0, 97),
    'j0378': (97, 0, 97),
    'j0395': (127, 0, 157),
    'j0410': (126, 0, 217),
    'j0430': (64, 0, 255),
    'j0515': (22, 255, 0),
    'j0660': (232, 0, 0),
    'j0861': (97, 0, 0)
}

# make 4x3 grid using gridspec

fig = plt.figure(constrained_layout=True, figsize=(11, 5), dpi=100)
gs = fig.add_gridspec(nrows=3, ncols=4) #, height_ratios=[1], width_ratios=[2,2,2])
# gs.update(left=0.08, right=0.98, bottom=0.10, top=0.98, wspace=0.4, hspace=0.03)

# add loop for plots in gs
for row in range(3):
    for col in range(4):
        ax = fig.add_subplot(gs[row, col])
        # label = 'Width: {}\nHeight: {}'.format(widths[col], heights[row])
        # ax.annotate(label, (0.1, 0.5), xycoords='axes fraction', va='center')

plt.savefig('plots/completeness.png', dpi=100)
plt.close()