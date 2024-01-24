#!/usr/bin/env python3

import argparse
import os
import json
from splusdata import splusdata
import pandas as pd
import multiprocessing as mp


def load_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=" ".join([
        'Get S-PLUS DR4 photometry.']))
    parser.add_argument('--work_dir', type=str,
                        help='Working directory. Default: current directory.',
                        default=os.getcwd())
    parser.add_argument('--fields_list', type=str,
                        help='List of fields to get the photometry.',
                        default='DR4_pointings.csv')
    parser.add_argument('--ncores', type=int,
                        help='Number of cores to use. Default: 1.',
                        default=1)
    parser.add_argument('--clobber', action='store_true',
                        help='Overwrite existing files.')

    args = parser.parse_args()
    return args


def get_dr4_photometry(field, sfilter='r', mode='dual_auto', save=True, output='output.csv'):
    """
    Get the photometry of a given field and filter from the DR4.
    Parameters
    ----------
    field : str
        Name of the field.
    filter : str
        Name of the filter.
    save : bool
        If True, save the output as a csv file.
    Returns
    -------
    df : pandas.DataFrame
        Pandas dataframe with the photometry of the field.
    """
    user, pswd = get_user_credentials()
    try:
        conn = splusdata.Core(user, pswd)
    except Exception as e:
        print(e)
        print('Error connecting to the database. Using anonymous connection.')
        conn = splusdata.Core()

    mode = mode.split('_')
    if mode[0] in 'dual':
        columns = 'RA,DEC'
    elif mode[0] in ['single', 'psf']:
        columns = f"RA_{sfilter},DEC_{sfilter}"

    if mode[0] in ['pdf']:
        query = f"SELECT {columns},{sfilter}_{mode[1]},e_{sfilter}_{mode[1]} FROM idr4_{mode[0]}.idr4_{mode[0]}_{sfilter} WHERE Field='{field}'"
    else:
        query = f"SELECT {columns},{sfilter}_{mode[1]},e_{sfilter}_{mode[1]},SEX_FLAGS_{sfilter} FROM idr4_{mode[0]}.idr4_{mode[0]}_{sfilter} WHERE Field='{field}'"
    print('Querying the database...')
    try:
        data = conn.query(query)
    except Exception as e:
        print(e)
        print(f'Error querying the database for {field} in {sfilter} filter.')
        data = pd.DataFrame()

    if data is None:
        print(f'No data for {field} in {sfilter} filter.')
    else:
        if save:
            if len(data) > 0:
                try:
                    data.write(output, format='csv', overwrite=True)
                    print(f"{output} successfully saved.")
                except Exception as e:
                    print(e)
                    print(f'Error saving {output}')
    return


def get_user_credentials():
    """
    Get the user credentials from the environment variables.
    Returns
    -------
    user : str
        Username.
    pswd : str
        Password.
    """
    print('Getting user credentials...')
    json_file = os.path.join(os.environ.get('HOME'), '.splus/credentials.json')
    if os.path.exists(json_file):
        json_data = json.load(open(json_file))
        user = json_data['splus.cloud']['user']
        pswd = json_data['splus.cloud']['password']
    else:
        user = os.environ.get('SPLUS_USER')
        pswd = os.environ.get('SPLUS_PSWD')

    return user, pswd


def process_field(args):
    field, mode, workdir, filters, save, clobber = args
    for f in filters:
        output = os.path.join(workdir, mode, f'{field}_{f}.csv')
        if os.path.exists(output) and not clobber:
            print(f'{output} already exists.')
            continue
        else:
            print(f'Getting {field} in {f} filter.')
            get_dr4_photometry(
                field, sfilter=f, mode=mode, save=save, output=output)


def main_run(fields, filters, workdir='./', modes=['dual_auto'], save=True, ncores=1, clobber=False):
    """
    Get the photometry of a list of fields and filters from the DR4.
    Parameters
    ----------
    fields : pandas.DataFrame
        Pandas dataframe with the list of fields.
    filters : list
        List of filters.
    save : bool
        If True, save the output as a csv file.
    Returns
    -------
    df : pandas.DataFrame
        Pandas dataframe with the photometry of the fields.
    """
    for mode in modes:
        if not os.path.exists(os.path.join(workdir, mode)):
            os.makedirs(os.path.join(workdir, mode))

    pool = mp.Pool(processes=ncores)
    field_args = [(field, mode, workdir, filters, save, clobber)
                  for field in fields['Field'] for mode in modes]
    pool.map(process_field, field_args)
    pool.close()
    pool.join()


if __name__ == '__main__':
    args = load_args()
    workdir = args.work_dir
    fields_list = args.fields_list
    ncores = args.ncores
    clobber = args.clobber
    fields = pd.read_csv(os.path.join(workdir, fields_list))
    filters = ['u', 'j0378', 'j0395', 'j0410', 'j0430',
               'g', 'j0515', 'r', 'j0660', 'i', 'j0861', 'z']
    # modes = ['single_auto', 'dual_auto', 'dual_PStotal', 'psf_psf']
    modes = ['single_auto', 'dual_auto', 'dual_PStotal']
    main_run(fields, filters, workdir=workdir,
             modes=modes, save=True, ncores=ncores, clobber=clobber)
