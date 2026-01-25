import requests
import json
import os
from astropy.io import fits 
import matplotlib.pyplot as plt 
import numpy as np
import time
import requests
import pyvo
import pandas as pd
from astropy.wcs import WCS
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from matplotlib.patches import Circle


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

def get_job_results(sub_id):
    url=f"http://nova.astrometry.net/api/submissions/{sub_id}"
    
    for attempt in range(5):
        try:
            response = requests.get(url)
            break  # Success!
        except requests.exceptions.ConnectionError as e:
             print(f"Connection error: {e}, retrying in 5 seconds...")
             time.sleep(5)
    else:
         print("Failed after several attempts.")

    try:
         jobs = response.json().get('jobs')
    except Exception:
         jobs = None
    print("Jobs:", jobs)
    return jobs
         

def get_astro_results(job_id):
    url = f"http://nova.astrometry.net/api/jobs/{job_id}/calibration/"
    response = requests.get(url)
    try:
        result = response.json()
        # ra = result.get("ra")
        # dec = result.get("dec")
        # radius = result.get("radius")
        # return [round(ra, 3), round(dec, 3), round(radius * 0.4, 3)]
        return result

    except Exception as e:
        print("Error parsing calibration results:", e)
        return None
    

def wait_for_calibration(job_id, timeout=180):
    url = f"http://nova.astrometry.net/api/jobs/{job_id}/calibration/"
    
    for _ in range(timeout):
        r = requests.get(url).json()

        # Calibration is ready
        if "ra" in r:
            return r

        time.sleep(2)

    return None





# def query_apass_to_csv(ra_center, dec_center, radius_deg, output_csv="apass_subset.csv"):
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


# Usage:
if __name__ == "__main__":
    api_key = os.environ.get("ASTRO_LOGIN")
    session_key = login_to_astrometry(api_key)
    print(session_key)
    subid_key = upload_fits_file("july25_mvn_red40.fits", session_key)
    # subid_key = os.environ.get("GREEN_SUBID")
    print("subdid key: ", subid_key)
    #job_id = get_job_results(subid_key)[0]
    # print(job_id)

    job_id = wait_for_job(subid_key) 
    print("Job ID:", job_id)

    astro_results = wait_for_calibration(job_id)
    print(astro_results)

    apply_calibration_to_fits("july25_mvn_red40.fits", "wcs_solution.fits", astro_results)



# Load your WCS‑corrected FITS file
hdul = fits.open("wcs_solution.fits")
data = hdul[0].data
w = WCS(hdul[0].header)

vmin, vmax = np.percentile(data, [5, 99])

fig = plt.figure(figsize=(10,8))
ax = plt.subplot(projection=w)
ax.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
ax.set_xlabel("RA")
ax.set_ylabel("Dec")
plt.show()

ra = 198.509684
dec = 36.683235

# Convert sky → pixel
x_pix, y_pix = w.world_to_pixel_values(ra, dec)
print(x_pix, y_pix)

fig = plt.figure(figsize=(10,8))
ax = plt.subplot(projection=w)

vmin, vmax = np.percentile(data, [5, 99])
ax.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)

# Plot the star
ax.plot(x_pix, y_pix, marker="o", markersize=12,
        markeredgecolor="yellow", markerfacecolor="none", linewidth=2)

ax.set_xlabel("RA")
ax.set_ylabel("Dec")
plt.title("Star overlay check")
plt.show()



# Load image
data = fits.getdata("wcs_solution.fits")

# Estimate background
mean, median, std = sigma_clipped_stats(data, sigma=3.0)

# Detect stars
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
sources = daofind(data - median)

print(sources)



# Stretch the image for visibility
vmin, vmax = np.percentile(data, [5, 99])

fig, ax = plt.subplots(figsize=(10, 8))
ax.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)

# Plot each detected star
for star in sources:
    x = star['xcentroid']
    y = star['ycentroid']
    circ = Circle((x, y), radius=8, edgecolor='yellow', facecolor='none', linewidth=1.5)
    ax.add_patch(circ)

plt.title("Detected Stars Overlay")
plt.xlabel("X pixel")
plt.ylabel("Y pixel")
plt.show()




# query_apass_to_csv(astro_results[0], astro_results[1], astro_results[2], "apass_subset.csv")
