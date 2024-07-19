#!/usr/bim/env python3

import os
import json
from splusdata import splusdata
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astroquery.vizier import Vizier


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
                        help='Field to query', default='STRIPE82-0119')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save plot to file')

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
    _r_auto_file = os.path.join(args.workdir, f'{args.field}_r_auto.csv')
    _i_auto_file = os.path.join(args.workdir, f'{args.field}_i_auto.csv')
    if not os.path.isfile(_r_auto_file) or not os.path.isfile(_i_auto_file):
        user, pswd = get_user_credentials()
        conn = splusdata.Core(user, pswd)
    if os.path.isfile(_r_auto_file):
        print(f'Reading {args.field} r-band photometry from {_r_auto_file}')
        dfautor = pd.read_csv(_r_auto_file)
    else:
        print('Querying S-PLUS database for filter r...')
        dfautor = conn.query(
            f"SELECT * FROM idr4_dual.idr4_dual_r WHERE Field='{args.field}'")
        print(f"Saving {len(dfautor)} objects to {args.workdir}/{_r_auto_file}")
        dfautor = dfautor.to_pandas()
        dfautor.to_csv(_r_auto_file, index=False)
    if os.path.isfile(_i_auto_file):
        print(f'Reading {args.field} i-band photometry from {_i_auto_file}')
        dfautoi = pd.read_csv(_i_auto_file)
    else:
        print('Querying S-PLUS database for filter i...')
        user, pswd = get_user_credentials()
        conn = splusdata.Core(user, pswd)
        dfautoi = conn.query(
            f"SELECT * FROM idr4_dual.idr4_dual_i WHERE Field='{args.field}'")
        print(f"Saving {len(dfautoi)} objects to {args.workdir}/{_i_auto_file}")
        dfautoi = dfautoi.to_pandas()
        dfautoi.to_csv(_i_auto_file, index=False)

    mask = (dfautor['r_PStotal'] > 0) & (dfautoi['i_PStotal'] > 0)
    mask &= (dfautor['r_PStotal'] < 30) & (dfautoi['i_PStotal'] < 30)

    dfauto = pd.concat([dfautor[mask].reset_index(drop=True),
                        dfautoi[mask].reset_index(drop=True)],
                       axis=1)

    return dfauto


def get_psf_photometry(args):
    _r_psf_file = os.path.join(args.workdir, f'{args.field}_r_psf.csv')
    _i_psf_file = os.path.join(args.workdir, f'{args.field}_i_psf.csv')
    if not os.path.isfile(_r_psf_file) or not os.path.isfile(_i_psf_file):
        user, pswd = get_user_credentials()
        conn = splusdata.Core(user, pswd)
    if os.path.isfile(_r_psf_file):
        print(f'Reading {args.field} r-band photometry from {_r_psf_file}')
        dfpsfr = pd.read_csv(_r_psf_file)
    else:
        print('Querying S-PLUS database for filter r_psf...')
        dfpsfr = conn.query(
            f"SELECT * FROM idr4_psf.idr4_psf_r WHERE Field='{args.field}'")
        print(f"Saving {len(dfpsfr)} objects to {_r_psf_file}")
        dfpsfr = dfpsfr.to_pandas()
        dfpsfr.to_csv(_r_psf_file, index=False)
    if os.path.isfile(_i_psf_file):
        print(f'Reading {args.field} i-band photometry from {_i_psf_file}')
        dfpsfi = pd.read_csv(_i_psf_file)
    else:
        print('Querying S-PLUS database for filter i_psf...')
        dfpsfi = conn.query(
            f"SELECT * FROM idr4_psf.idr4_psf_i WHERE Field='{args.field}'")
        print(f"Saving {len(dfpsfi)} objects to {_i_psf_file}")
        dfpsfi = dfpsfi.to_pandas()
        dfpsfi.to_csv(_i_psf_file, index=False)

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


def get_m2_phot():
    """
    Query Vizier for NGC7089 (M2) photometry
    """
    # Will use E. Hartmann catalogue made using SPLUS data
    # the catalogue is 515/4191/NGC7089
    # https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/MNRAS/515/4191/ngc7089
    # will query for all objects, but only columns RAJ2000, DEJ2000, rmag, e_rmag, imag, e_imag

    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'rmag', 'e_rmag', 'imag', 'e_imag'],
               catalog='J/MNRAS/515/4191/ngc7089')
    v.ROW_LIMIT = 20000
    cache_path = os.path.join(args.workdir, 'ngc7089.vizier')
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)

    v.cache_location = cache_path

    m2data = v.get_catalogs('J/MNRAS/515/4191/ngc7089')

    return m2data


def match_splus_m2(splus, m2cat, gaia_data=None):
    """
    Match SPLUS photometry with M2 photometry
    """
    splus_coords = SkyCoord(
        ra=splus['RA'], dec=splus['DEC'], unit=('deg', 'deg'))
    m2_coords = SkyCoord(ra=m2cat['RAJ2000'],
                         dec=m2cat['DEJ2000'], unit=('hourangle', 'deg'))

    idx, d2d, _ = splus_coords.match_to_catalog_3d(m2_coords)
    sep_constraint = d2d < 1 * u.arcsec

    matched_m2 = pd.concat([m2cat[idx][sep_constraint].to_pandas().reset_index(drop=True),
                           splus[sep_constraint].reset_index(drop=True)],
                           axis=1)

    if gaia_data is None:
        return matched_m2

    # add proper motion
    gaia_coords = SkyCoord(ra=gaia_data[0]['RA_ICRS'],
                           dec=gaia_data[0]['DE_ICRS'], unit=('deg', 'deg'))
    matched_m2_coords = SkyCoord(
        ra=matched_m2['RA'], dec=matched_m2['DEC'], unit=('deg', 'deg'))
    idx, d2d, _ = matched_m2_coords.match_to_catalog_3d(gaia_coords)
    sep_constraint = d2d < 1 * u.arcsec
    matched_m2 = pd.concat([matched_m2[sep_constraint].reset_index(drop=True),
                           gaia_data[0][idx][sep_constraint].to_pandas().reset_index(drop=True)],
                           axis=1)

    return matched_m2


def get_proper_motion(args):
    """
    Get proper motion from Gaia DR3 up to 0.5 deg around the central coords
    """
    # get Gaia Dr3 data from Vizier (I/355/gaiadr3)
    # citation: https://ui.adsabs.harvard.edu/abs/2016A%26A...595A...1G/abstract
    # and
    # https://ui.adsabs.harvard.edu/abs/2023A%26A...674A...1G/abstract
    central_coords = SkyCoord(args.ra, args.dec, unit=('hourangle', 'deg'))
    radius = 0.5 * u.deg
    gaia = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE'],
                  catalog='I/355/gaiadr3')
    gaia.ROW_LIMIT = 999999
    cache_path = os.path.join(args.workdir, 'gaiaedr3.vizier')
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)

    gaia.cache_location = cache_path

    gaia_data = gaia.query_region(central_coords, radius=radius)

    return gaia_data


def plot_photometry(args):

    central_coords = SkyCoord(args.ra, args.dec, unit=('hourangle', 'deg'))

    dfauto = get_auto_photometry(args)
    # drop repeated columns
    dfauto = dfauto.loc[:, ~dfauto.columns.duplicated()]
    dfpsf = get_psf_photometry(args)
    # rename RA and DEC columns if RA or DEC is part of the name
    dfpsf.rename(columns={'RA_r': 'RA', 'DEC_r': 'DEC'}, inplace=True)
    gaia_data = get_proper_motion(args)
    # gaia_data = None

    m2data = get_m2_phot()[0]
    matched_m2_auto = match_splus_m2(dfauto, m2data, gaia_data=gaia_data)
    matched_m2_psf = match_splus_m2(dfpsf, m2data, gaia_data=gaia_data)
    # define a radial mask around central coords
    if gaia_data is not None:
        radii = 1.0  # deg
        mask_auto = central_coords.separation(SkyCoord(
            ra=matched_m2_auto['RA'], dec=matched_m2_auto['DEC'], unit=('deg', 'deg'))) < radii * u.deg
        mask_psf = central_coords.separation(SkyCoord(
            ra=matched_m2_psf['RA'], dec=matched_m2_psf['DEC'], unit=('deg', 'deg'))) < radii * u.deg
        mask_proper_motion = matched_m2_auto['pmRA'] >= (
            3.51 - 0.5)  # mas/year
        # mean proper motions were taken from
        # https://academic.oup.com/mnras/article/482/4/5138/5173087?login=true
        mask_proper_motion &= matched_m2_auto['pmRA'] <= (3.51 + 0.5)
        mask_proper_motion &= matched_m2_auto['pmDE'] >= (-2.16 - 0.5)
        mask_proper_motion &= matched_m2_auto['pmDE'] <= (-2.16 + 0.5)
        mask_auto &= mask_proper_motion
        mask_proper_motion = matched_m2_psf['pmRA'] >= (
            3.51 - 0.5)
        mask_proper_motion &= matched_m2_psf['pmRA'] <= (3.51 + 0.5)
        mask_proper_motion &= matched_m2_psf['pmDE'] >= (-2.16 - 0.5)
        mask_proper_motion &= matched_m2_psf['pmDE'] <= (-2.16 + 0.5)
        mask_psf &= mask_proper_motion
    else:
        radii = 0.25
        mask_auto = central_coords.separation(SkyCoord(
            ra=matched_m2_auto['RA'], dec=matched_m2_auto['DEC'], unit=('deg', 'deg'))) < radii * u.deg
        mask_psf = central_coords.separation(SkyCoord(
            ra=matched_m2_psf['RA'], dec=matched_m2_psf['DEC'], unit=('deg', 'deg'))) < radii * u.deg

    field_lims = {'x': [central_coords.ra.value - 0.5,
                        central_coords.ra.value + 0.5],
                  'y': [central_coords.dec.value - 0.5,
                        central_coords.dec.value + 0.5]}
    colourmag_lims = {'x': [-0.52, 0.9], 'y': [21, 12]}
    labelsize = 16
    ticksizes = 14
    maskflags = (dfauto['SEX_FLAGS_r'] <= 3)
    maskflags &= (dfauto['RA'] >= field_lims['x'][0])
    maskflags &= (dfauto['RA'] < field_lims['x'][1])
    maskflags &= (dfauto['DEC'] >= field_lims['y'][0])
    maskflags &= (dfauto['DEC'] < field_lims['y'][1])
    maskpsf = dfpsf['RA'] >= field_lims['x'][0]
    maskpsf &= dfpsf['RA'] <= field_lims['x'][1]
    maskpsf &= dfpsf['DEC'] >= field_lims['y'][0]
    maskpsf &= dfpsf['DEC'] <= field_lims['y'][1]

    fig = plt.figure(figsize=(8, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 3])
    ax1 = plt.subplot(gs[0])
    ax1.scatter(dfauto['RA'][maskflags], dfauto['DEC'][maskflags], s=1,
                c='grey', alpha=0.25, label='Field objects')
    ax1.scatter(matched_m2_auto['RA'][mask_auto], matched_m2_auto['DEC'][mask_auto],
                marker='*', s=8, c='b', label='M2')
    ax1.set_title(r'$\mathrm{Dual\ mode\ [PStotal]}$', fontsize=labelsize)
    ax1.set_xlim(field_lims['x'])
    ax1.set_ylim(field_lims['y'])
    ax1.set_xlabel(r'$\mathrm{RA\ [deg]}$', fontsize=labelsize)
    ax1.set_ylabel(r'$\mathrm{DEC\ [deg]}$', fontsize=labelsize)
    ax1.tick_params(labelsize=ticksizes)
    ax1.legend(loc='upper left')

    ax2 = plt.subplot(gs[1])
    ax2.scatter(dfpsf['RA'], dfpsf['DEC'], s=1, c='grey',
                alpha=0.25, label='Field objects')
    ax2.scatter(matched_m2_psf['RA'][mask_psf], matched_m2_psf['DEC'][mask_psf],
                marker='*', s=8, c='m', label='M2')
    ax2.set_title(r'$\mathrm{PSF\ photometry}$', fontsize=labelsize)
    ax2.set_xlim(field_lims['x'])
    ax2.set_ylim(field_lims['y'])
    ax2.set_xlabel(r'$\mathrm{RA\ [deg]}$', fontsize=labelsize)
    ax2.set_yticklabels(ax2.get_yticklabels(), visible=False)
    ax2.tick_params(labelsize=ticksizes)
    ax2.legend(loc='upper left')

    ax3 = plt.subplot(gs[2])
    # ax3.grid()
    ax3.scatter(dfauto['r_PStotal'][maskflags] - dfauto['i_PStotal'][maskflags],
                dfauto['r_PStotal'][maskflags], s=1, c='gray', alpha=0.25)
    ax3.scatter(matched_m2_auto['r_PStotal'][mask_auto] - matched_m2_auto['i_PStotal'][mask_auto],
                matched_m2_auto['r_PStotal'][mask_auto], marker='*', s=8, c='b')
    ax3.set_xlim(colourmag_lims['x'])
    ax3.set_ylim(colourmag_lims['y'])
    ax3.set_xlabel(r'$r - i$', fontsize=labelsize)
    ax3.set_ylabel(r'$r$', fontsize=labelsize)
    ax3.tick_params(labelsize=ticksizes)

    ax4 = plt.subplot(gs[3])
    # ax4.grid()
    ax4.scatter(dfpsf['r_psf'][maskpsf] - dfpsf['i_psf'][maskpsf],
                dfpsf['r_psf'][maskpsf], s=1, c='grey', alpha=0.25)
    ax4.scatter(matched_m2_psf['r_psf'][mask_psf] - matched_m2_psf['i_psf'][mask_psf],
                matched_m2_psf['r_psf'][mask_psf], marker='*', s=8, c='m')
    ax4.set_xlim(colourmag_lims['x'])
    ax4.set_ylim(colourmag_lims['y'])
    ax4.set_xlabel(r'$r - i$', fontsize=labelsize)
    ax4.set_yticklabels(ax4.get_yticklabels(), visible=False)
    ax4.tick_params(labelsize=ticksizes)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02, hspace=0.15)

    if args.saveplot:
        print('Saving plot to {}'.format(os.path.join(
            args.workdir, 'pstotalxpsf_photometry.png')))
        plt.savefig(os.path.join(args.workdir, 'pstotalxpsf_photometry.png'))
        plt.close()
    else:
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
