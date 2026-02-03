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
import astropy.units as u

from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground


from matplotlib.patches import Circle
import matplotlib.pyplot as plt 
import matplotlib.lines as mlines
import numpy as np
import pyvo
import pandas as pd

from photutils.detection import DAOStarFinder

global status_message
status_message = "Waiting for user input..."


def lsrl(x, y):

    N1 = len(x)

    sum_xsq1 = 0
    for terms in x:
        sum_xsq1 += terms ** 2

    sum_xy1 = 0
    for i in range(len(x)):
        sum_xy1 += x[i] * y[i]

    sum_x1 = 0
    for terms in x:
        sum_x1 += terms
        
    sum_y1 = 0
    for terms in y:
        sum_y1 += terms 
    
    q1 = np.array([[sum_y1], [sum_xy1]])
    p1 = np.array([[N1, sum_x1], [sum_x1, sum_xsq1]])
    matrix1 = np.dot(np.linalg.inv(p1), q1)
    
    m1 = matrix1[1]
    b1 = matrix1[0]

    sd1 = 0
    
    for i in range(len(x)):
        sd1 += ((y[i] - (m1 * x[i] + b1)) ** 2)
        
    sd1 /= len(x)
    sd1 = sd1 ** 0.5
    
    x2 = []
    y2 = []
    
    for i in range(len(x)):
        if (((y[i] - (m1 * x[i] + b1)) ** 2) < 4 * (sd1 ** 2)):
            x2.append(x[i])
            y2.append(y[i])
            
            
    N2 = len(x2)
    
    sum_xsq2 = 0
    for terms in x2:
        sum_xsq2 += terms ** 2

    sum_xy2 = 0
    for i in range(len(x2)):
        sum_xy2 += x2[i] * y2[i]

    sum_x2 = 0
    for terms in x2:
        sum_x2 += terms
        
    sum_y2 = 0
    for terms in y2:
        sum_y2 += terms 
    
    q2 = np.array([[sum_y2], [sum_xy2]])
    p2 = np.array([[N2, sum_x2], [sum_x2, sum_xsq2]])
    matrix2 = np.dot(np.linalg.inv(p2), q2)
    
    m2 = matrix2[1]
    b2 = matrix2[0]

    sd2 = 0
    
    for i in range(len(x2)):
        sd2 += ((y2[i] - (m2 * x2[i] + b2)) ** 2)
        
    sd2 /= len(x2)
    sd2 = sd2 ** 0.5
        
    #slope, y-int, error
    
    return (m2, b2, sd2)

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
    global status_message
    url = f"http://nova.astrometry.net/api/submissions/{sub_id}"
    for i in range(timeout):
        r = requests.get(url).json()
        jobs = r.get("jobs", [])
        if jobs and jobs[0] is not None:
            print("Job found:", jobs[0])

            
            #status_message = f"Job found: {jobs[0]}"
            status_message = "Job ID achieved...please wait for astrometry results"

            return jobs[0]
        print("Waiting for job...", i)        
        status_message = f"Waiting to achieve Job ID... {i}"

        time.sleep(2)
    raise TimeoutError("Job did not appear in time.")
         

def query_apass_to_csv(ra_center, dec_center, radius_deg, output_csv="apass_subset.csv"):
    global status_message

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
    df = df.dropna(subset=["ra", "dec", "mag_g", "mag_r"])
    df.to_csv(output_csv, index=False)
    print(f"Saved {output_csv}")
    num_rows = len(df)
    status_message = f"Saved csv file as {output_csv}..."
    status_message = f"Done!"

    return num_rows - 2



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

    global status_message

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

    status_message = f"WCS written to {output_fits}..."



def target_flux(x, y, radius, image):
    """
    Same as flux(), but for a single target star.
    """
    data = fits.getdata(image)

    position = np.array([[x, y]])

    aperture = CircularAperture(position, r=radius)
    annulus = CircularAnnulus(position, r_in=radius+3, r_out=radius+6)

    aper_phot = aperture_photometry(data, aperture)
    ann_phot = aperture_photometry(data, annulus)

    bkg_mean = ann_phot['aperture_sum'] / annulus.area
    final_flux = aper_phot['aperture_sum'] - bkg_mean * aperture.area

    return np.array([max(final_flux[0], 0)])


def flux(x, y, radius, image):
    """
    Perform circular aperture photometry with an annulus background.
    x, y: arrays of centroid positions (floats)
    radius: aperture radius in pixels
    image: FITS filename
    """
    data = fits.getdata(image)

    # Convert to Nx2 array of positions
    positions = np.column_stack((x, y))

    # Circular aperture
    aperture = CircularAperture(positions, r=radius)

    # Background annulus (3 px wide)
    annulus = CircularAnnulus(positions, r_in=radius+3, r_out=radius+6)

    # Perform photometry
    aper_phot = aperture_photometry(data, aperture)
    ann_phot = aperture_photometry(data, annulus)

    # Compute background per pixel
    bkg_mean = ann_phot['aperture_sum'] / annulus.area

    # Subtract background from aperture flux
    final_flux = aper_phot['aperture_sum'] - bkg_mean * aperture.area

    # Replace negative flux with zero
    final_flux = np.where(final_flux > 0, final_flux, 0)

    return np.array(final_flux)

def detect_stars(image, fwhm=3.0, sigma=3.0, threshold=5.0):
    data = fits.getdata(image)
    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std)
    sources = daofind(data - median)
    return sources

def full_calibration_with_subid(image, wcs_image_name, subid_key):
    global status_message

    # 1. Login to Astrometry.net
    api_key = os.environ.get("ASTRO_LOGIN")
    session_key = login_to_astrometry(api_key)
    print("Session:", session_key)

    
    #status_message = f"Session: {session_key}"
    status_message = f"Session key achieved...please wait for submission ID"


    # 2. Upload image if no subid_key provided
    if subid_key is None:
        subid_key = upload_fits_file(image, session_key)
        print("Submission ID:", subid_key)
        #status_message = f"Submission ID: {subid_key}"
        status_message = f"Submission ID achieved...please wait for job ID"
    else:
        subid_key = subid_key
        print("Submission ID:", subid_key)
        #status_message = f"Submission ID: {subid_key}"
        status_message = f"Submission ID achieved...please wait for job ID"

    # 3. Wait for job to finish
    job_id = wait_for_job(subid_key)
    print("Job ID:", job_id)
    #status_message = f"Job ID: {job_id}"
    status_message = f"Job ID achieved...please wait for astrometry results"

    # 4. Retrieve calibration results
    astro_results = wait_for_calibration(job_id)
    print("Astrometry Results:", astro_results)
    #status_message = f"Astrometry Results: {astro_results}"
    status_message = f"Astrometry Results achieved...please wait for WCS application"

    # 5. Apply WCS to the green image
    image = apply_calibration_to_fits(image, wcs_image_name, astro_results)

    # 6. Query APASS around the solved coordinates
    ra = astro_results["ra"]
    dec = astro_results["dec"]
    radius = round(astro_results["radius"] * 0.6, 1)

    num_rows = query_apass_to_csv(ra, dec, radius, "apass_subset.csv")
    status_message = f"Done!"
    
    return num_rows

def full_calibration(image, wcs_image_name):
    global status_message

    # 1. Login to Astrometry.net
    api_key = os.environ.get("ASTRO_LOGIN")
    session_key = login_to_astrometry(api_key)
    print("Session:", session_key)

    #status_message = f"Session: {session_key}"
    status_message = f"Session key achieved...please wait for submission ID...please be patient as this will be the longest step!"

    # 2. Upload image if no subid_key provided
    subid_key = upload_fits_file(image, session_key)
    print("Submission ID:", subid_key)

    #status_message = f"Submission ID: {subid_key}"
    status_message = f"Submission ID achieved...please wait for job ID"

    # 3. Wait for job to finish
    job_id = wait_for_job(subid_key)
    print("Job ID:", job_id)

    #status_message = f"Job ID: {job_id}"
    status_message = f"Job ID achieved...please wait for astrometry results"

    # 4. Retrieve calibration results
    astro_results = wait_for_calibration(job_id)
    print("Astrometry Results:", astro_results)
    #status_message = f"Astrometry Results: {astro_results}"
    status_message = f"Astrometry Results achieved...please wait for WCS application"

    

    # 5. Apply WCS to the green image
    image = apply_calibration_to_fits(image, wcs_image_name, astro_results)

    # 6. Query APASS around the solved coordinates
    ra = astro_results["ra"]
    dec = astro_results["dec"]
    radius = round(astro_results["radius"] * 0.5, 1)

    num_rows = query_apass_to_csv(ra, dec, radius, "apass_subset.csv")
    status_message = f"Done!"

    return num_rows
    
    

def magnitudes(csv_file, green_image, red_image, n, RA, DEC):
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

    '''
    x_pixel_g = x_pixel_g.astype(int)
    y_pixel_g = y_pixel_g.astype(int)
    x_pixel_r = x_pixel_r.astype(int)
    y_pixel_r = y_pixel_r.astype(int)
    '''
   # === 1. Remove APASS stars outside the image ===
    ny, nx = image_data.shape

    inside = (
        (x_pixel_g >= 0) & (x_pixel_g < nx) &
        (y_pixel_g >= 0) & (y_pixel_g < ny)
    )

    x_pixel_g = x_pixel_g[inside]
    y_pixel_g = y_pixel_g[inside]
    x_pixel_r = x_pixel_r[inside]
    y_pixel_r = y_pixel_r[inside]
    g = g[inside]
    r = r[inside]

    print("APASS stars inside image:", len(x_pixel_g))


    # === 2. Detect stars in the green image ===
    sources = detect_stars(green_image, threshold=3.0)
    det_x = np.array(sources['xcentroid'])
    det_y = np.array(sources['ycentroid'])

    print("Detected stars:", len(det_x))

    # === 3. Convert detected centroids to RA/Dec ===
    det_ra, det_dec = w_g.all_pix2world(det_x, det_y, 1)
    det_coords = SkyCoord(det_ra * u.deg, det_dec * u.deg, frame="icrs")

    # === 4. Build SkyCoord for APASS catalog ===
    cat_ra = data['ra']
    cat_dec = data['dec']
    cat_coords = SkyCoord(cat_ra * u.deg, cat_dec * u.deg, frame="icrs")

    # === 5. Match each detected star to nearest APASS star ===
    idx, sep2d, _ = det_coords.match_to_catalog_sky(cat_coords)

    # VERY strict matching: only keep stars within 1 arcsec
    max_sep = 4.5 * u.arcsec
    good = sep2d < max_sep

    print("Matched APASS stars:", np.sum(good))

    # Keep only matched detections
    matched_det_x = det_x[good]
    matched_det_y = det_y[good]
    matched_cat_ix = idx[good]

    # Pull standard magnitudes from APASS
    g = data['mag_g'][matched_cat_ix]
    r = data['mag_r'][matched_cat_ix]



    # === 5. Use detected centroids for flux extraction ===
    x_flux = matched_det_x
    y_flux = matched_det_y

    # === 6. Filter out faint APASS stars (critical!) ===
    bright = (g < 16.5)   # or 15.5 if you want even cleaner stars
    x_flux = x_flux[bright]
    y_flux = y_flux[bright]
    g = g[bright]
    r = r[bright]

    print("Final calibration stars after brightness filter:", len(g))

    
    green_flux = flux(x_flux, y_flux, 8, green_image)
    red_flux = flux(x_flux, y_flux, 8, red_image)


    valid_flux_g = []
    valid_flux_r = []


    mask = (
    (green_flux > 0) &
    np.isfinite(green_flux) &
    (red_flux > 0) &
    np.isfinite(red_flux)
    )   

    valid_flux_g = green_flux[mask]
    valid_flux_r = red_flux[mask]
    g = g[mask]
    r = r[mask]
    
    valid_flux_g = np.array(valid_flux_g)
    valid_flux_r = np.array(valid_flux_r)

    inst_g = -2.5 * np.log10(valid_flux_g)
    inst_r = -2.5 * np.log10(valid_flux_r)

    #X1
    inst_g_r = inst_g - inst_r
    
    # Y1 & X2
    st_g_r = g - r
    
    # Y2
    g_offset = g - inst_g

    ##Specifically solving for the target object now.

    target_pixel_g_x, target_pixel_g_y = w_g.all_world2pix(RA, DEC, 1)
    target_pixel_r_x, target_pixel_r_y = w_r.all_world2pix(RA, DEC, 1)


    '''
    target_pixel_g_x = target_pixel_g_x.astype(int)
    target_pixel_g_y = target_pixel_g_y.astype(int)
    target_pixel_r_x = target_pixel_r_x.astype(int)
    target_pixel_r_y = target_pixel_r_y.astype(int)
    '''
    target_flux_g = target_flux(target_pixel_g_x, target_pixel_g_y, 9, green_image)
    target_flux_r = target_flux(target_pixel_r_x, target_pixel_r_y, 9, red_image)

    target_g_inst_mag = -2.5 * np.log10(target_flux_g)
    target_r_inst_mag = -2.5 * np.log10(target_flux_r)

    m1_b1 = lsrl(inst_g_r, st_g_r)
    m2_b2 = lsrl(st_g_r, g_offset)

    new_std = m1_b1[0]*(inst_g_r) + m1_b1[1]
    new_std_inst = m2_b2[0]*new_std + m2_b2[1] 
    label_text = f'Tgr = {round(m1_b1[0][0], 4)} \n Cgr = {round(m1_b1[1][0], 4)}'
    legend_entry = mlines.Line2D([], [], color='none', label=label_text)
    
    standard_g_r_target = m1_b1[0]*(target_g_inst_mag - target_r_inst_mag) + m1_b1[1]
    standard_g_target = m2_b2[0]*standard_g_r_target + m2_b2[1] + target_g_inst_mag
    
    



    standard_r_target = standard_g_target - standard_g_r_target

    error_g, error_r = m1_b1[2], m2_b2[2]
    
    return (standard_g_target, standard_r_target, error_g, error_r, calibration_num) 



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
    g_file = None
    r_file = None
    g_path = None 
    r_path = None

    standard_g_target = None 
    standard_r_target = None 
    error_g = None 
    error_r = None


    ra_decimal = None
    dec_decimal = None
    ra_hms = None
    dec_dms = None

    upload_folder = "uploads" 
    os.makedirs(upload_folder, exist_ok=True)

    if request.method == "POST":
        # Single input case
        user_text = request.form.get("calibration_input")

        

        # Both-input case
        g_text = request.form.get("g_link_both")
        r_text = request.form.get("r_link_both")

        g_file = request.files.get("g_file")
        r_file = request.files.get("r_file")

        ra_decimal = request.form.get("ra_decimal")
        dec_decimal = request.form.get("dec_decimal")
        ra_hms = request.form.get("ra_hms")
        dec_dms = request.form.get("dec_dms")

        # Option A: decimal degrees
        if ra_decimal and dec_decimal:
            try:
                ra_deg = float(ra_decimal)
                dec_deg = float(dec_decimal)
            except ValueError:
                pass

        # Option B: HMS/DMS
        elif ra_hms and dec_dms:
            try:
                coord = SkyCoord(ra_hms, dec_dms, unit=(u.hourangle, u.deg))
                ra_deg = coord.ra.deg
                dec_deg = coord.dec.deg
            except Exception:
                pass

        # If neither option was provided
        if ra_deg is None or dec_deg is None:
            return render_template(
                "object_calibration.html",
                error_message="Please enter RA/Dec in either decimal or HMS/DMS format."
            )


        

    if g_file and g_file.filename.strip():
        g_path = os.path.join(upload_folder, g_file.filename)
        g_file.save(g_path)

    if r_file and r_file.filename.strip():
        r_path = os.path.join(upload_folder, r_file.filename)
        r_file.save(r_path)


    # If user uploaded files, use those
    if g_path and r_path:
        #full_calibration_with_subid(g_path, "wcs_green_solution.fits", os.environ.get("GREEN_SUBID"))
        #num_rows = full_calibration_with_subid(r_path, "wcs_red_solution.fits", os.environ.get("RED_SUBID"))
        full_calibration(g_path, "wcs_green_solution.fits")
        full_calibration(r_path, "wcs_red_solution.fits")

        standard_g_target, standard_r_target, error_g, error_r, calibration_num = magnitudes(
            "apass_subset.csv",
            "wcs_green_solution.fits",
            "wcs_red_solution.fits",
            num_rows,
            ra_deg,
            dec_deg
        )

    # If user typed paths instead, use those
    elif g_text and r_text:
        #full_calibration_with_subid(g_text, "wcs_green_solution.fits", os.environ.get("GREEN_SUBID"))
        #num_rows = full_calibration_with_subid(r_text, "wcs_red_solution.fits", os.environ.get("RED_SUBID"))
        full_calibration(g_text, "wcs_green_solution.fits")
        num_rows = full_calibration(r_text, "wcs_red_solution.fits")

        standard_g_target, standard_r_target, error_g, error_r, calibration_num = magnitudes(
            "apass_subset.csv",
            "wcs_green_solution.fits",
            "wcs_red_solution.fits",
            num_rows,
            ra_deg,
            dec_deg
        )

    return render_template(
        "object_calibration.html",
        g_text=g_text,
        r_text=r_text,
        standard_g_target=standard_g_target,
        standard_r_target=standard_r_target,
        error_g=error_g,
        error_r=error_r,
    )



@app.route("/status")
def status():
    return status_message




if __name__ == "__main__":
    app.run(debug=True)


