#!/usr/bin/env python3

import argparse
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import dask.dataframe as dd
import itertools
import multiprocessing as mp

os.environ["DASK_MEMORY_LIMIT"] = "32G"


def load_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=" ".join([
                                     'Plot the photometry of a given field.']))
    parser.add_argument('--work_dir', type=str,
                        help='Working directory. Default: current directory.',
                        default=os.getcwd())
    parser.add_argument('--fields_list', type=str,
                        help='List of fields to get the photometry.',
                        default='DR4_pointings.csv')
    parser.add_argument('--nprocs', type=int,
                        help='Number of processes to use. Default: 1.',
                        default=1)
    parser.add_argument('--npoints', type=int,
                        help='Number of points to plot. Use -1 to plot the entire table. Default: 10000',
                        default=10000)
    parser.add_argument('--nbins', type=int,
                        help='Number of bins for shadow regions. Default: 20.',
                        default=20)
    parser.add_argument('--plot', action='store_true',
                        help='Plot the photometry.')
    parser.add_argument('--savefig', action='store_true',
                        help='Save the plot.')
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite existing files.')

    return parser.parse_args()


def load_csv_data(csv_file):
    """Load the data from a csv file."""
    field = csv_file.split('/')[-1].split('_')[0]
    try:
        data = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f'File {csv_file} not found.')
        return pd.DataFrame(columns=['RA', 'DEC', 'mag', 'e_mag', 'Field'])
    data['Field'] = [field]*len(data)
    data.columns = ['RA', 'DEC', 'mag', 'e_mag', 'Field']

    return dd.from_pandas(data, npartitions=1)


def process_csvfiles(workdir, modes, filters, nprocs, args):
    for mode, f in itertools.product(modes, filters):
        outputcsv = os.path.join(workdir, f'{mode}/all_{f}.csv')
        if os.path.exists(outputcsv) and not args.clobber:
            print(f'File {outputcsv} exists. Skipping.')
            continue
        print(f'Processing mode {mode} and filter {f}...')
        list_csvs = glob.glob(os.path.join(workdir, mode, f'*_{f}.csv'))

        pool = mp.Pool(nprocs)
        ddfs = pool.map(load_csv_data, list_csvs)
        pool.close()
        pool.join()

        alldf = dd.concat(ddfs, ignore_index=True)

        print(f'Saving file {outputcsv}.')
        alldf.to_csv(outputcsv, index=False)


def plot_photometry(workdir, modes, filters, args):
    """Plot the photometry."""
    colours = {'u': 'indigo', 'j0378': 'darkviolet', 'j0395': 'navy',
               'j0410': 'b', 'j0430': 'dodgerblue',  'j0515': 'lime',
               'g': 'turquoise', 'r': 'limegreen', 'j0660': 'y',
               'i': 'darkorange', 'j0861': 'orangered', 'z': 'darkred'}
    mode_names = {'dual_auto': 'Dual', 'dual_PStotal': 'PStotal',
                  'single_auto': 'Single', 'psf_psf': 'PSF'}

    fig, ax = plt.subplots(12, 4, figsize=(10, 12))
    ax = ax.ravel()
    iter_modes_filters = itertools.product(filters, modes)
    for i, (f, m) in enumerate(iter_modes_filters):
        mode_dir = os.path.join(workdir, f'{m}')
        allcsv = os.path.join(mode_dir, f'all_{f}.csv')
        print(i, f, m, allcsv)
        if os.path.exists(allcsv):
            if args.npoints == -1:
                try:
                    df = dd.read_csv(allcsv, blocksize=25e6)
                except ValueError:
                    print(
                        f'An error occurred while reading the file {allcsv}.')
                    continue
            else:
                try:
                    df = dd.read_csv(allcsv, blocksize=25e6).head(
                        n=args.npoints)
                except ValueError:
                    print(
                        f'An error occurred while reading the file {allcsv}.')
                    continue
        else:
            print(f'File {allcsv} not found.')
            continue
        df['mag'] = df['mag'].astype(float)
        df['e_mag'] = df['e_mag'].astype(float)
        df['mag'] = df['mag'].round(3)
        df['e_mag'] = df['e_mag'].round(3)

        mask = (df['mag'] > 10) & (df['mag'] < 30)
        mask &= (df['e_mag'] > 0) & (df['e_mag'] < 10.9)
        completeness = len(df[mask]) / len(df)
        xlims = [10.1, 26.9]
        ylims = [-0.5, 10.9]
        xticks = [12, 16, 20, 24]
        yticks = [1, 3, 5, 7, 9]

        df = df[mask]
        df = df.sort_values(by='mag')
        df['mag_bin'] = pd.qcut(df['mag'], q=25, labels=False)
        percentiles = [.16, .50, .84]
        mag_bins = df.groupby('mag_bin')[['mag', 'e_mag']].quantile(
            percentiles).unstack()
        mag_intersect = mag_bins['mag'][percentiles[1]][
            mag_bins['e_mag'][percentiles[1]] < 0.3619].max()

        ax[i].scatter(df['mag'], df['e_mag'], s=1, alpha=0.4, color='gray')
        print(min(df['mag']), max(df['mag']), min(df['e_mag']),
              max(df['e_mag']))
        ax[i].plot(mag_bins['mag'][percentiles[0]],
                   mag_bins['e_mag'][percentiles[0]],
                   'o-', lw=2, color=colours[f], ms=2)
        ax[i].fill_between(mag_bins['mag'][percentiles[0]],
                           mag_bins['e_mag'][percentiles[0]],
                           mag_bins['e_mag'][percentiles[2]],
                           color=colours[f], alpha=0.3)
        if mag_intersect <= max(df['mag']):
            ax[i].axvline(mag_intersect, color='k', ls='--', lw=1)
            if i in [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47]:
                ax[i].text(mag_intersect+0.1, 0.05, r'$\mathrm{%.2f}$' %
                           mag_intersect, fontsize=6)
            else:
                ax[i].text(mag_intersect+0.1, 9, r'$\mathrm{%.2f}$' %
                           mag_intersect, fontsize=6)

        if i in [44, 45, 46, 47]:
            ax[i].set_xlabel(r'$\mathrm{Mag}$', fontsize=12)
        else:
            ax[i].set_xticklabels(ax[i].get_xticklabels(), visible=False)
        if i in [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44]:
            ax[i].set_ylabel(r'$\mathrm{e_{Mag}}$', fontsize=12)
        if i in [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47]:
            ax2 = ax[i].twinx()
            if f in ['u', 'g', 'r', 'i', 'z']:
                fname = f'{f}'
            else:
                fname = f'{f}'.upper()
            ax2.set_ylabel(r"$%s$" % fname, rotation=90,
                           labelpad=10, fontsize=12)
            ax2.set_yticklabels(ax2.get_yticklabels(), visible=False)
            ax2.set_yticks([])
            ax[i].set_xlim(10.1, 23)
            ax[i].set_xticks(xticks)
            ax[i].set_ylim(-0.02, 0.51)
        else:
            ax[i].set_xlim(xlims)
            ax[i].set_ylim(ylims)
            ax[i].set_xticks(xticks)
            ax[i].set_yticks(yticks)

        if i in [0, 1, 2, 3]:
            ax[i].set_title(f'{mode_names[m]}', fontsize=12)

        ax[i].grid(alpha=0.5)
        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        ax[i].text(0.05, 0.9, r'%i (%.3f)' % (len(df), completeness),
                   transform=ax[i].transAxes,
                   fontsize=6, verticalalignment='top', bbox=bbox_props)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.02)
    if args.savefig:
        figpath = os.path.join(workdir, 'dr4-photxerrors.png')
        print(f'Saving figure {figpath}.')
        plt.savefig(figpath, format='png', dpi=300)
    plt.show()
    plt.close()


def main(args):

    workdir = args.work_dir
    modes = ['dual_auto', 'dual_PStotal', 'single_auto', 'psf_psf']
    filters = ['u', 'j0378', 'j0395', 'j0410', 'j0430', 'g',  'j0515', 'r',
               'j0660', 'i', 'j0861', 'z']
    nprocs = args.nprocs

    process_csvfiles(workdir, modes, filters, nprocs, args)

    if args.plot:
        plot_photometry(workdir, modes, filters, args)


if __name__ == '__main__':
    args = load_args()

    main(args)
