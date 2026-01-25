from flask import Flask, render_template, request

import requests
import json
import os

import time
import random
import math

from astropy.io import fits 
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import ascii

from matplotlib.patches import Circle
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
import numpy as np
import pyvo
import pandas as pd

from photutils.detection import DAOStarFinder
from test_4 import query_apass_to_csv




def login_to_astrometry(api_key: str) -> str:
	url = 'http://nova.astrometry.net/api/login'
	payload = {'request-json': json.dumps({"apikey": api_key})}
	response = requests.post(url, data=payload)
	try:
		session_key = response.json().get('session')
	except Exception:
		session_key = None
	return session_key


def upload_fits_file(file_path, session_token, url="http://nova.astrometry.net/api/upload"):

    request_json = {
        "publicly_visible": "y",
        "allow_modifications": "d",
        "session": session_token,
        "allow_commercial_use": "d"
    }

    with open(file_path, "rb") as f:
        files = {
            "request-json": (
                None,
                json.dumps(request_json),
                "text/plain"
            ),
            "file": (
                file_path,
                f,
                "application/octet-stream"
            )
        }
        response = requests.post(url, files=files)

    # print("Status code:", response.status_code)
    # print("Response:", response.text)
    try:
        subid = response.json().get("subid")
    except Exception:
        subid = None
    return subid

def wait_for_job(sub_id, timeout=180):
    url = f"http://nova.astrometry.net/api/submissions/{sub_id}"
    for i in range(timeout):
        r = requests.get(url).json()
        jobs = r.get("jobs", [])
        if jobs and jobs[0] is not None:
            print("Job found:", jobs[0])
            return jobs[0]
        print("Waiting for job...", i)
        time.sleep(2)
    raise TimeoutError("Job did not appear in time.")
         

def query_apass_to_csv(ra_center, dec_center, radius_deg, output_csv="apass_subset.csv"):
    """
    Query APASS DR10 catalog using GAVO TAP and save results to a CSV file.
    Args:
        ra_center (float): Right Ascension of center (degrees)
        dec_center (float): Declination of center (degrees)
        radius_deg (float): Search radius (degrees)
        output_csv (str): Output CSV filename
    """
    service = pyvo.dal.TAPService("https://dc.g-vo.org/tap")
    query = f"""
    SELECT * 
    FROM apass.dr10
    WHERE 1 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra_center}, {dec_center}, {radius_deg})
    )
    """
    result = service.search(query)
    df = result.to_table().to_pandas()
    df.to_csv(output_csv, index=False)
    print(f"Saved {output_csv}")



def wait_for_calibration(job_id, timeout=180):
    url = f"http://nova.astrometry.net/api/jobs/{job_id}/calibration/"
    
    for _ in range(timeout):
        r = requests.get(url).json()

        # Calibration is ready
        if "ra" in r:
            return r

        time.sleep(2)

    return None
   






    """
    Query APASS DR10 catalog using GAVO TAP and save results to a CSV file.
    Args:
        ra_center (float): Right Ascension of center (degrees)
        dec_center (float): Declination of center (degrees)
        radius_deg (float): Search radius (degrees)
        output_csv (str): Output CSV filename
    """
    service = pyvo.dal.TAPService("https://dc.g-vo.org/tap")
    query = f"""
    SELECT * 
    FROM apass.dr10
    WHERE 1 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra_center}, {dec_center}, {radius_deg})
    )
    """
    result = service.search(query)
    df = result.to_table().to_pandas()
    df.to_csv(output_csv, index=False)
    print(f"Saved {output_csv}")


def apply_calibration_to_fits(input_fits, output_fits, cal):
    """
    Apply Astrometry.net calibration JSON to a FITS file.
    """

    # Extract calibration values
    ra = cal["ra"]
    dec = cal["dec"]
    pixscale = cal["pixscale"] / 3600.0   # arcsec → deg
    orient = np.deg2rad(cal["orientation"])
    parity = cal["parity"]

    # Load FITS
    hdul = fits.open(input_fits)
    hdr = hdul[0].header
    data = hdul[0].data
    ny, nx = data.shape

    # Compute CRPIX
    crpix1 = nx / 2.0
    crpix2 = ny / 2.0

    # CD matrix
    sgn = -1 if parity < 0 else 1
    cd11 = -pixscale * np.cos(orient)
    cd12 =  pixscale * np.sin(orient)
    cd21 =  pixscale * np.sin(orient) * sgn
    cd22 =  pixscale * np.cos(orient) * sgn

    # Build WCS
    w = WCS(naxis=2)
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [crpix1, crpix2]
    w.wcs.cd = np.array([[cd11, cd12],
                         [cd21, cd22]])
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    # Update FITS header
    hdr.update(w.to_header())

    # Save new FITS
    hdul.writeto(output_fits, overwrite=True)
    hdul.close()

    print(f"WCS written to {output_fits}")


def flux(x, y, radius, image):
    hdu = fits.open(image)
    data = hdu[0].data
    flux = np.array([])

    for i in range(len(x)):
        x_low = x[i] - radius - 1 
        x_high = x[i] + radius
        y_low = y[i] - radius - 1
        y_high = y[i] + radius 
        
        ring = data[y_low-1:y_high+1, x_low-1:x_high+1]
        
        

        if (
            x_low - 1 < 0 or x_high + 1 >= data.shape[1] or
            y_low - 1 < 0 or y_high + 1 >= data.shape[0]
        ):
            hdul = fits.open(image)
            wcs = WCS(hdul[0])
            ra, dec = wcs.all_pix2world(x[i], y[i], 0)
            print(f"Skipping star at ({ra}, {dec}): out of bounds")
            continue
        
        
        #plt.imshow(ring, cmap='gray_r')
        #plt.show()

        #print (ring)
        
        background = []
        
        for i in range(len(ring[0])):
            background.append(ring[0][i])
            background.append(ring[len(ring[0]) - 1][i])
            
            if i != 0 and i != len(ring[0] - 1):
                background.append(ring[i][0])
                background.append(ring[i][len(ring[0])-1])

        # print(background)
                
        median = np.median(background)
        
        star = data[y_low:y_high, x_low:x_high]
        
        
        flux_num = 0
        
        for i in range(len(star)):
            for j in range(len(star[0])):
                flux_num += (star[i][j] - median)
        if flux_num < 0:
            flux_num = 0
        
        flux = np.append(flux, flux_num)
    return flux

def magnitudes(csv_file, green_image, red_image, n):
    data = ascii.read(csv_file, format='csv')
    hdul_g = fits.open(green_image)
    hdul_r = fits.open(red_image)
    wcs_g_h = WCS(hdul_g[0].header)
    wcs_r_h = WCS(hdul_r[0].header)

    fits_file_path = green_image
    with fits.open(fits_file_path) as hdul:
        image_data = hdul[0].data

    w_g = WCS(green_image)
    w_r = WCS(red_image)

    col_length = len(data['ra'])
    calibration_num = random.sample(range(0, col_length - 1), n)
    calibration_num = sorted(calibration_num)

    ra_list = np.array([])
    dec_list = np.array([])

    x_pixel_g = np.array([])
    y_pixel_g = np.array([])
    x_pixel_r = np.array([])
    y_pixel_r = np.array([])
    mag_g = np.array([])
    mag_r = np.array([])
    g = np.array([])
    r = np.array([])

    j = 0

    for i in range(col_length):
        if j < len(calibration_num) and i == calibration_num[j]:
            ra = data['ra'][i]
            dec = data['dec'][i]
            x_g, y_g = w_g.all_world2pix(ra, dec, 1)
            x_r, y_r = w_r.all_world2pix(ra, dec, 1)
            ra_list = np.append(ra_list, ra)
            dec_list = np.append(dec_list, dec)
            x_pixel_g = np.append(x_pixel_g, x_g)
            y_pixel_g = np.append(y_pixel_g, y_g)
            x_pixel_r = np.append(x_pixel_r, x_r)
            y_pixel_r = np.append(y_pixel_r, y_r)
            g = np.append(g, data['mag_g'][i])
            r = np.append(r, data['mag_r'][i])

            j += 1
            
    x_pixel_g = x_pixel_g.astype(int)
    y_pixel_g = y_pixel_g.astype(int)
    x_pixel_r = x_pixel_r.astype(int)
    y_pixel_r = y_pixel_r.astype(int)
    
    green_flux = flux(x_pixel_g, y_pixel_g, 8, green_image)
    red_flux = flux(x_pixel_r, y_pixel_r, 8, red_image)

    return green_flux, red_flux, ra_list, dec_list, g, r



def full_calibration_with_subid(image, wcs_image_name, subid_key):

    # 1. Login to Astrometry.net
    api_key = os.environ.get("ASTRO_LOGIN")
    session_key = login_to_astrometry(api_key)
    print("Session:", session_key)

    # 2. Upload image if no subid_key provided
    if subid_key is None:
        subid_key = upload_fits_file(image, session_key)
        print("Submission ID:", subid_key)
    else:
        subid_key = subid_key
        print("Submission ID:", subid_key)

    # 3. Wait for job to finish
    job_id = wait_for_job(subid_key)
    print("Job ID:", job_id)

    # 4. Retrieve calibration results
    astro_results = wait_for_calibration(job_id)
    print("Astrometry Results:", astro_results)

    # 5. Apply WCS to the green image
    image = apply_calibration_to_fits(image, wcs_image_name, astro_results)

    # 6. Query APASS around the solved coordinates
    ra = astro_results["ra"]
    dec = astro_results["dec"]
    radius = round(astro_results["radius"] * 0.5, 1)

    query_apass_to_csv(ra, dec, radius, "apass_subset.csv")

def full_calibration(image, wcs_image_name):

    # 1. Login to Astrometry.net
    api_key = os.environ.get("ASTRO_LOGIN")
    session_key = login_to_astrometry(api_key)
    print("Session:", session_key)

    # 2. Upload image if no subid_key provided
    subid_key = upload_fits_file(image, session_key)
    print("Submission ID:", subid_key)

    # 3. Wait for job to finish
    job_id = wait_for_job(subid_key)
    print("Job ID:", job_id)

    # 4. Retrieve calibration results
    astro_results = wait_for_calibration(job_id)
    print("Astrometry Results:", astro_results)

    # 5. Apply WCS to the green image
    image = apply_calibration_to_fits(image, wcs_image_name, astro_results)

    # 6. Query APASS around the solved coordinates
    ra = astro_results["ra"]
    dec = astro_results["dec"]
    radius = round(astro_results["radius"] * 0.5, 1)

    query_apass_to_csv(ra, dec, radius, "apass_subset.csv")




app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def home():
    user_text = None

    if request.method == "POST":
        user_text = request.form.get("user_input")

    return render_template("index.html", user_text=user_text)

@app.route("/object_calibration", methods=["GET", "POST"])
def object_calibration():
    user_text = None
    g_text = None
    r_text = None

    green_flux = None 
    red_flux = None 
    ra_list = None 
    dec_list = None 
    g = None 
    r = None

    if request.method == "POST":
        # Single input case
        user_text = request.form.get("calibration_input")

        # Both-input case
        g_text = request.form.get("g_link_both")
        r_text = request.form.get("r_link_both")



    if g_text and r_text:

        #full_calibration(g_text, "wcs_green_solution.fits")
        #full_calibration(r_text, "wcs_red_solution.fits")
        full_calibration_with_subid(g_text, "wcs_green_solution.fits", os.environ.get("GREEN_SUBID"))
        full_calibration_with_subid(r_text, "wcs_red_solution.fits", os.environ.get("RED_SUBID"))
        green_flux, red_flux, ra_list, dec_list, g, r = magnitudes("apass_subset.csv", "wcs_green_solution.fits", "wcs_red_solution.fits", 10)


    return render_template(
        "object_calibration.html",
        user_text=user_text,
        g_text=g_text,
        r_text=r_text,
        green_flux=green_flux, 
        red_flux=red_flux, 
        ra_list=ra_list, 
        dec_list=dec_list, 
        g=g, 
        r=r
    )






if __name__ == "__main__":
    app.run(debug=True)


