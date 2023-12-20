#!/usr/bim/env python3

import os
import json
from splusdata import splusdata
from astropy.coordinates import SkyCoord
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='Get S-PLUS stamp')
    parser.add_argument('--workdir', type=str,
                        help='Working directory', default=os.getcwd())
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
        import pdb
        pdb.set_trace()

        hdus_stamp = conn.stamp(coords.ra.value, coords.dec.value, 30, "r")

        hdus_stamp.writeto(os.path.join(args.workdir, args.stamp_name))

    else:
        print('Reading stamp from {}'.format(
            os.path.join(args.workdir, args.stamp_name)))
        hdus_stamp = fits.open(os.path.join(args.workdir, args.stamp_name))

    return hdus_stamp


def main(args):
    hdus_stamp = get_splus_stamp(args)
    print(hdus_stamp)


if __name__ == '__main__':
    args = get_args()
    if not os.path.exists(args.workdir):
        os.makedirs(args.workdir)
    main(args)
