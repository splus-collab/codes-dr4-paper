#!/usr/bin/env python3

import argparse
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import itertools


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

    return parser.parse_args()


def main(args):

    workdir = args.work_dir
    fields_list = args.fields_list
    modes = ['dual_auto', 'dual_PStotal', 'single_auto', 'psf_psf']
    filters = ['u', 'g', 'r', 'i', 'z', 'j0378', 'j0395', 'j0410', 'j0430',
               'j0515', 'j0660', 'j0861']
    fields = pd.read_csv(os.path.join(workdir, fields_list))

    fig, ax = plt.subplots(12, 4, figsize=(10, 10))
    ax = ax.ravel()
    iter_modes_filters = itertools.product(filters, modes)
    all_data = {}
    for i, (f, m) in enumerate(iter_modes_filters):
        print(i, f, m)
        mode_dir = os.path.join(workdir, f'{m}')
        all_data[m] = pd.DataFrame(columns=['RA', 'DEC', f"{f}_{m.split('_')[1]}",
                                            f"e_{f}_{m.split('_')[1]}", 'Field'])
        for field in fields['Field'][:5]:
            print('reading file', os.path.join(mode_dir, f'{field}_{f}.csv'))
            try:
                df = pd.read_csv(os.path.join(mode_dir, f'{field}_{f}.csv'))
            except FileNotFoundError:
                print(f'File {field}_{f}.csv not found.')
                continue
            df['Field'] = [field]*len(df)
            df.columns = all_data[m].columns
            all_data[m] = pd.concat([all_data[m], df], ignore_index=True)
        # all_data = all_data.drop_duplicates(subset=['RA', 'DEC'])
        mask = all_data[m][f"e_{f}_{m.split('_')[1]}"] < 30
        mask &= all_data[m][f"e_{f}_{m.split('_')[1]}"] > -30
        mask &= all_data[m][f"{f}_{m.split('_')[1]}"] < 30
        mask &= all_data[m][f"{f}_{m.split('_')[1]}"] > 0

        ax[i].scatter(all_data[m][f"{f}_{m.split('_')[1]}"][mask],
                      all_data[m][f"e_{f}_{m.split('_')[1]}"][mask],
                      s=0.5, alpha=0.5)
        # ax[i].set_xlabel(f"{f}_{m.split('_')[1]}")
        # ax[i].set_ylabel(f"e_{f}_{m.split('_')[1]}")
        # ax[i].set_title(f'{f}_{m.split("_")[1]}')
    import pdb
    pdb.set_trace()

    plt.tight_layout()
    # plt.savefig(os.path.join(workdir, 'photxerror.png'), dpi=300)
    plt.show()
    plt.close()


if __name__ == '__main__':
    args = load_args()

    main(args)
