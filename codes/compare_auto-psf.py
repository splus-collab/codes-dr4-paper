#!/usr/bim/env python3

import os
import json
from splusdata import splusdata
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser(description='Get S-PLUS stamp')
    parser.add_argument('--workdir', type=str,
                        help='Working directory', default=os.getcwd())
    parser.add_argument('--getstamp', action='store_true',
                        help='Get the stamp from S-PLUS')
    parser.add_argument('--ra', type=str,
                        help='RA of the object', default='21 33 27.02')
    parser.add_argument('--dec', type=str,
                        help='DEC of the object', default='-00 49 23.7')
    parser.add_argument('-size', '--size', type=float,
                        help='Size of the stamp')
    parser.add_argument('-filter', '--filter', type=str,
                        help='Filter of the stamp')
    parser.add_argument('--stamp_name', type=str,
                        help='Name of the stamp', default='stamp.fits')
    parser.add_argument('--field', type=str,
                        help='Photometry file', default='photometry.csv')

    args = parser.parse_args()
    return args


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


def get_splus_stamp(args):
    if not os.path.exists(os.path.join(args.workdir, args.stamp_name)):
        print('Downloading stamp from S-PLUS...')
        user, pswd = get_user_credentials()
        try:
            conn = splusdata.Core(user, pswd)
        except Exception as e:
            print(e)
            print('Error connecting to the database. Using anonymous connection.')
            conn = splusdata.Core()

        coords = SkyCoord(ra=args.ra, dec=args.dec, unit=('hourangle', 'deg'))

        hdus_stamp = conn.stamp(coords.ra.value, coords.dec.value, 30, "r")

        hdus_stamp.writeto(os.path.join(args.workdir, args.stamp_name))

    else:
        print('Reading stamp from {}'.format(
            os.path.join(args.workdir, args.stamp_name)))
        hdus_stamp = fits.open(os.path.join(args.workdir, args.stamp_name))

    return hdus_stamp


def get_auto_photometry(args):
    dfautor = pd.read_csv(os.path.join(
        args.workdir, f'dual_PStotal/{args.field}_r.csv'))
    dfautoi = pd.read_csv(os.path.join(
        args.workdir, f'dual_PStotal/{args.field}_i.csv'))

    mask = (dfautor['r_PStotal'] > 0) & (dfautoi['i_PStotal'] > 0)
    mask &= (dfautor['r_PStotal'] < 30) & (dfautoi['i_PStotal'] < 30)

    dfauto = pd.concat([dfautor[mask].reset_index(drop=True),
                        dfautoi[mask].reset_index(drop=True)],
                       axis=1)

    return dfauto


def get_psf_photometry(args):
    dfpsfr = pd.read_csv(os.path.join(
        args.workdir, f'psf_psf/{args.field}_r.csv'))
    dfpsfi = pd.read_csv(os.path.join(
        args.workdir, f'psf_psf/{args.field}_i.csv'))

    coords_r = SkyCoord(ra=dfpsfr['RA_r'],
                        dec=dfpsfr['DEC_r'], unit=('deg', 'deg'))
    coords_i = SkyCoord(ra=dfpsfi['RA_i'],
                        dec=dfpsfi['DEC_i'], unit=('deg', 'deg'))

    idx, d2d, _ = coords_r.match_to_catalog_3d(coords_i)
    sep_constraint = d2d < 1 * u.arcsec

    matched_dfpsfr = dfpsfr[sep_constraint]
    matched_dfpsfi = dfpsfi.iloc[idx][sep_constraint]

    mask = (matched_dfpsfr['r_psf'] > 0).values & (
        matched_dfpsfi['i_psf'] > 0).values
    mask &= (matched_dfpsfr['r_psf'] < 30).values & (
        matched_dfpsfi['i_psf'] < 30).values

    dfpsf = pd.concat([matched_dfpsfr[mask].reset_index(drop=True),
                      matched_dfpsfi[mask].reset_index(drop=True)], axis=1)

    return dfpsf


def plot_photometry(args):

    dfauto = get_auto_photometry(args)
    dfpsf = get_psf_photometry(args)
    import pdb
    pdb.set_trace()

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax.ravel()
    ax[0, 0].scatter(dfauto['RA'], dfauto['DEC'], s=1, c='k')
    ax[0, 1].scatter(dfpsf['RA_r'], dfpsf['DEC_r'], s=1, c='k')
    ax[1, 0].scatter(dfauto['r_PStotal'] - dfauto['i_PStotal'],
                     dfauto['r_PStotal'], s=1, c='k')
    ax[1, 1].scatter(dfpsf['r_psf'] - dfpsf['i_psf'],
                     dfpsf['r_psf'], s=1, c='k')
    plt.show()


def main(args):
    if args.getstamp:
        hdus_stamp = get_splus_stamp(args)
        print(hdus_stamp)

    plot_photometry(args)


if __name__ == '__main__':
    args = get_args()
    if not os.path.exists(args.workdir):
        os.makedirs(args.workdir)
    main(args)
