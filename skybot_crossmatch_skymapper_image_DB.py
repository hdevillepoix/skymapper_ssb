#!/bin/bash
"""
Retrieve possible asteroids seen by SkyMapper
Using Skybot (http://vo.imcce.fr/webservices/skybot/?conesearch)
"""

import os
from pathlib import Path
import numpy as np

from astroquery.imcce import Skybot

import astropy.units as u
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord


__author__ = "Hadrien Devillepoix"
__license__ = "MIT"
__version__ = "1.0"


MAG_CUT_SAFETY_TERM = 3 * u.mag

SKYMAPPER_FOV = 1.8 * u.degree # fov radius
SKYMAPPER_MPC_CODE = Q61 # technically wrong: 260 = Siding Spring DSS, Q61 = PROMPT, Siding Spring

# image pointing catalogue
SKYMAPPER_IMAGE_CATALOGUE = 'data/skymapper_DR3_image_table_prerelease.csv'

OUTPUT_DEST = 'data/output'

OVERWRITE = False



def main(ifile=SKYMAPPER_IMAGE_CATALOGUE):

    
    print('loading input image catalogue {}'.format(ifile))
    skymapper_images = Table.read(ifile, format='ascii.csv')

    print(f'{len(skymapper_images)} images loaded')
    print(skymapper_images.colnames)
    


    # convert some key fields to astropy ecosystem objects
    skymapper_images['MJD'] = Time(skymapper_images['date'], format='mjd') + TimeDelta(skymapper_images["exp_time"]*u.second)/2
    skymapper_images['skycoord'] = SkyCoord(ra=skymapper_images['ra']*u.rad,
                                        dec=skymapper_images['dcl']*u.rad)
    

    for image in skymapper_images:
        if OVERWRITE or not os.path.isfile(res_fname(image["image_id"])):
            process_image(image)


def res_fname(image_id):
    """
    generate an output filename for a given image_id
    """
    return os.path.join(OUTPUT_DEST, "skymapper_dr3_image_" + str(image_id) + "_potential_solar_system_bodies.csv")
    
def process_image(image):
    
    # hacky way to determine lim mag
    if image["exp_time"] > 50:
        survey_type = 'deep'
        lim_mag = 21 * u.mag
    else:
        survey_type = 'shallow'
        lim_mag = 19 * u.mag
        
    try:
        skybot_cs = cone_search_skybot(image)
    except RuntimeError as e:
        if 'No solar system object was found in the requested FOV' in str(e):
            print('image {} : {}'.format(image["image_id"], str(e)))
            # create empty file
            Path(res_fname(image["image_id"])).touch()
            return
        else:
            print('Server error for image {} :'.format(image["image_id"]), file=sys.stderr)
            print(e, file=sys.stderr)
            return
    
    # C/2002 O4, this is the only object that messes things up
    problematic_comet_mask = (skybot_cs['Name'] == 'C/2002 O4')
    if np.sum(problematic_comet_mask) > 0:
        print('image {} contained problematic comet'.format(image["image_id"]))
        skybot_cs = skybot_cs[~problematic_comet_mask]
        if len(skybot_cs) < 1:
            print('problematic comet was the only results for image {}'.format(image["image_id"]))
            Path(res_fname(image["image_id"])).touch()
            return

    
    all_asts_num = len(skybot_cs)
    
    # activate mag cut
    mag_cut = False
    # only select asteroids that would be at worse 3 magnitudes fainter in the V band
    if mag_cut:
        skybot_cs = skybot_cs[skybot_cs["V"] < lim_mag + MAG_CUT_SAFETY_TERM]
    
    print(f'found {len(skybot_cs)} / {all_asts_num} potential solar system bodies in SkyMapper image {image["image_id"]} ({survey_type})')
    

    
    keep_cols = ['Number',
                'Name',
                'RA',
                'DEC',
                'Type',
                'V',
                'posunc',
                'RA_rate',
                'DEC_rate',
                'geodist',
                'heliodist',
                'alpha',
                'elong']
    
    # clean up columns
    remove_cols = list(set(skybot_cs.colnames) - set(keep_cols))

    
    skybot_cs.remove_columns(remove_cols)
    skybot_cs_reduced_cols = skybot_cs
    
    skybot_cs_reduced_cols["skymapper_image_id"] = image["image_id"]
    
    
    if len(skybot_cs) > 0:
        try:
            save_res(skybot_cs_reduced_cols, image["image_id"])
        except TypeError as e:
            print('Problem with return result from {} :'.format(image["image_id"]), file=sys.stderr)
            print(e, file=sys.stderr)
            for c in skybot_cs.colnames:
                print("{} : {}".format(c, skybot_cs[c]), file=sys.stderr)
    else:
        print('Image {} had matches ({}), but none made the magnitude cut'.format(image["image_id"]), file=sys.stderr)
        
        
    


def save_res(skybot_cs, image_id):
    ofname = res_fname(image_id)
    
    skybot_cs.write(ofname, format='ascii.csv', overwrite=True)
    print(f'results written to {ofname}')


def cone_search_skybot(image):
    try:
        res = Skybot.cone_search(image['skycoord'],
                            rad=SKYMAPPER_FOV,
                            epoch=image['MJD'],
                            location=SKYMAPPER_MPC_CODE,
                            position_error=15*u.arcsec)
    except RuntimeError as e:
        raise e
    
    return res
    
if __name__ == '__main__':
    import sys
    ifile = sys.argv[1]
    
    main(ifile)
