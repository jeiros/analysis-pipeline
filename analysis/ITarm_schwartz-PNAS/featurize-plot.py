"""Plot RMSD results

msmbuilder autogenerated template version 2
created 2017-07-03T15:49:19.699572
please cite msmbuilder in any publications
"""

import numpy as np
from matplotlib import pyplot as plt
from msmexplorer import plot_trace
from msmbuilder.io import load_trajs
from matplotlib.ticker import FuncFormatter
from plot_utils import cleanup_top_right_axes
from cycler import cycler
from random import shuffle
import pandas as pd
from plot_utils import figure_dims
import matplotlib.gridspec as gridspec
import seaborn.apionly as sns
import matplotlib.patches as mpatches
sns.set_style('ticks')
colors = ["#cfa0ce",
          "#63db56",
          "#7941cf",
          "#b4e741",
          "#ce51c1",
          "#cecf3c",
          "#4f2b7f",
          "#629f3b",
          "#7879d7",
          "#d7a835",
          "#36203a",
          "#b4d97f",
          "#da487a",
          "#69dca1",
          "#e14733",
          "#79d6d9",
          "#9c392b",
          "#76a3c9",
          "#d57b38",
          "#555e85",
          "#d6c37b",
          "#8c3964",
          "#47682a",
          "#d38580",
          "#53987d",
          "#603728",
          "#cdd0b9",
          "#2d4334",
          "#917333",
          "#8d8271"]


def wrapAngle(x):
    """Wraps an angle between -180 and 180 degrees"""
    x = (x + 180) % 360
    if x < 0:
        x += 360
    return x - 180


def constrainAngle(x):
    """Constrains an angle between 0 and 360 degrees"""
    x = x % 360
    if x < 0:
        x += 360
    return x


def plot_angle(data, N=50, title=None, ax1=None, ax2=None, color=None, wrap=True):
    if ax1 is None or ax2 is None:
        gs = gridspec.GridSpec(2, 6)
        ax1 = plt.subplot(gs[:1, :2], polar=True)
        ax2 = plt.subplot(gs[:1, 2:])

    if wrap:
        vf = np.vectorize(wrapAngle)
    else:
        vf = np.vectorize(constrainAngle)
    x = vf(data)

    sns.distplot(x, bins=N, ax=ax2, color=color, kde=True)
    radii, theta = np.histogram(x, bins=N, normed=True)
    ax1.set_yticklabels([])

    if wrap:
        ax1ticks = [0, 45, 90, 135, 180, -135, -90, -45]
        ax2ticks = list(range(-180, 180 + 45, 45))
        ax1.set_xticklabels(['{}°'.format(x) for x in ax1ticks])
        ax2.set_xlim(-180, 180)
        ax2.set_xticks(ax2ticks)
        ax2.set_xticklabels(['{}°'.format(x) for x in ax2ticks])

    else:
        ax2ticks = list(range(0, 360 + 45, 45))
        ax2.set_xlim(0, 360)
        ax2.set_xticks(ax2ticks)
        ax2.set_xticklabels(['{}°'.format(x) for x in ax2ticks])

    ax2.set_yticks([])
    ax2.set(xlabel='Angle', ylabel='Density')

    sns.despine(ax=ax2)
    width = (2 * np.pi) / N

    ax1.bar(np.deg2rad(theta[1:]), radii, width=width, color=color, alpha=.5)

    if title is not None:
        plt.suptitle(title)

    plt.tight_layout()

    f = plt.gcf()
    return f, (ax1, ax2)


if __name__ == '__main__':
    # Load
    meta, dtrajs = load_trajs('dtrajs')

    dihed_names = [
        'cTnI dihedral',
        'cTnT dihedral'
    ]
    data_dict = {}  # This will have keys, one for each sim type
    for t in meta['type'].unique():
        indexes = meta[meta['type'] == t].index  # The simulation indexes
        dihed_arr = np.concatenate([dtrajs.get(key) for key in indexes])
        data_dict[t] = dihed_arr

    n_diheds = list(data_dict.values())[0].shape[1]

    for i in range(n_diheds):
        n_types = 0
        patch_list = []
        for k, v in data_dict.items():
            data = np.rad2deg(v[:, i])
            p = mpatches.Patch(label=k, color=colors[n_types])
            patch_list.append(p)
            if n_types == 0:
                f, (ax1, ax2) = plot_angle(data, color=colors[n_types])
            else:
                f, (ax1, ax2) = plot_angle(data, ax1=ax1, ax2=ax2, color=colors[n_types])

            n_types += 1
            plt.legend(handles=patch_list, loc='best')
        f.set_size_inches(figure_dims(600, 0.9))
        plt.suptitle(dihed_names[i])
        f.savefig('dihed-no{}.pdf'.format(i))
