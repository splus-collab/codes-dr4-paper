#!/bin/python3.10

import argparse
import os
import sys
import json
from splusdata import splusdata as sd
import pandas as pd


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=" ".join([
        'Calculate the transmission curve or a given filtre.',
        'Estimate the central lambda from the FWHM of that filter.']))
    parser.add_argument('--work_dir', type=str, help='Working directory. Default: current directory.',
                        default=os.getcwd())
    parser.add_argument('--fields_list', type=str, help='List of fields to get the photometry.',
                        default='DR4_pointings.csv')

    # if len(sys.argv) <= 1:
    #     parser.print_help()
    #     sys.exit(0)

    args = parser.parse_args()
    return args


def get_dr4_photometry(field, sfilter='r', mode='dual_auto'):
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
        conn = sd.Core(user, pswd)
    except Exception as e:
        print(e)
        print('Error connecting to the database. Using anonymous connection.')
        conn = sd.Core()

    mode = mode.split('_')
    if mode[0] in 'dual':
        columns = 'RA,DEC'
    elif mode[0] in ['single', 'psf']:
        columns = f"RA_{sfilter},DEC_{sfilter}"

    query = f"SELECT {columns},{sfilter}_{mode[1]},e_{sfilter}_{mode[1]} FROM idr4_{mode[0]}.idr4_{mode[0]}_{sfilter} WHERE Field='{field}'"
    print('Querying the database...')
    try:
        data = conn.query(query)
    except Exception as e:
        print(e)
        print(f'Error querying the database for {field} in {sfilter} filter.')
        data = pd.DataFrame()
        import pdb
        pdb.set_trace()

    if data is None:
        import pdb
        pdb.set_trace()
    return data


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


def main_run(fields, filters, workdir='./', modes=['dual_auto'], save=True):
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

        for field in fields['Field']:
            for f in filters:
                output = os.path.join(workdir, mode, f'{field}_{f}.csv')
                if os.path.exists(output):
                    print(f'{output} already exists.')
                    continue
                else:
                    print(f'Getting {field} in {f} filter.')
                    df = get_dr4_photometry(field, sfilter=f, mode=mode)
                    if save:
                        if df is None:
                            print(f'No data for {field} in {f} filter.')
                            import pdb
                            pdb.set_trace()
                        if len(df) > 0:
                            try:
                                df.write(output, format='csv', overwrite=True)
                                print(f"{output} successfully saved.")
                            except Exception as e:
                                print(e)
                                print(f'Error saving {output}')


if __name__ == '__main__':
    args = get_args()
    # fields = pd.read_csv(os.path.join(args.work_dir, args.fields_list))
    workdir = '/storage/splus/DR4_phot/'
    fields_list = 'DR4_pointings.csv'
    fields = pd.read_csv(os.path.join(workdir, fields_list))
    filters = ['u', 'j0378', 'j0395', 'j0410', 'j0430',
               'g', 'j0515', 'r', 'j0660', 'i', 'j0861', 'z']
    modes = ['single_auto', 'dual_auto', 'dual_PStotal', 'psf_psf']
    main_run(fields, filters, workdir=workdir, modes=modes, save=True)
