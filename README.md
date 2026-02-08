# 🌌 Stellar Photometry & Calibration Portal  
Precision. Automation. Scientific Accuracy.

This project is a complete web platform for automated stellar photometry, WCS calibration, and magnitude standardization using Sloan **g/r** or Johnson **B/V** filters. It integrates:

- Astrometry.net plate solving  
- APASS DR10 catalog querying  
- Instrumental photometry with background subtraction  
- Color‑term and zero‑point calibration  
- Target‑object analysis  
- Star‑cluster calibration  
- Diagnostic plots for transparency and verification  

Whether you're analyzing a single star or an entire cluster, this tool handles the entire workflow end‑to‑end with scientific rigor and a clean, modern UI.

---

# 🚀 Quick Start Guide

This guide walks you through:

1. Cloning the repository  
2. Creating a virtual environment  
3. Installing all dependencies  
4. Setting up environment variables  
5. Running the website locally  

Everything is designed to be reproducible and beginner‑friendly.

---

# 1. Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

# 2. Clone the Repository
macOS / Linux:
```bash
python3 -m venv venv
source venv/bin/activate
```
Windows (PowerShell)
```bash
python -m venv venv
venv\Scripts\Activate
```

# 3. Install All Required Libraries
Install everything with:
```bash
pip install -r requirements.txt
```
If you ever need to regenerate the file:
```bash
pip freeze > requirements.txt
```

# 4. Set Up Your Environment Variables

This project requires an Astrometry.net API key to perform WCS plate solving.

## 1. Create a `.env` file

Create the file in the root of the project (same folder as `app.py`):

macOS / Linux:
    touch .env

Windows (PowerShell):
    New-Item -Path .env -ItemType File

## 2. Get your Astrometry.net API key

If you don’t have an API key, generate one at:

https://nova.astrometry.net

1. go to https://nova.astrometry.net
2. Press "sign in" on the top right corner 
3. Choose any account you would like to sign in with
4. Go to "Profile" seen on the left side of the screen 
5. You will see "my API key: _____________" (It will look like letters: asdfsdfgghjgho)
6. insert that API key like this: ASTRO_LOGIN=asdfsdfgghjgho

## 3. Add the key to `.env`

Open `.env` and paste:

ASTRO_LOGIN=YOUR_ASTROMETRY_API_KEY


Replace YOUR_ASTROMETRY_API_KEY with your actual key.

## 4. Ensure Flask loads the `.env` file

If needed, add this to the top of `app.py`:
```
from dotenv import load_dotenv
load_dotenv()
```. 
## 5. Verify it works

Run:
```
python app.py
```
Then inside your code, print:
```
print(os.environ.get("ASTRO_LOGIN"))
```
If it prints your key, the environment is configured correctly.





# 5. Run the Website
```bash
python app.py
```
Or:
```bash
flask run
```
You should see something like this:
Running on http://127.0.0.1:5000

Open that link in your browser


# 6. Using the Platform
⭐ Object Calibration
Enter RA/Dec (decimal or HMS/DMS)

Choose photometric system (Sloan g/r or Johnson B/V) 
*NOTE* Johnson B/V is unavailable in this current version

Upload FITS files or provide file paths

The system performs:

WCS solving

APASS DR10 querying

Star detection

Flux extraction

Color‑term calibration

Zero‑point calibration

Final standardized magnitudes

Results + plots appear automatically

⭐ Star Cluster Calibration
Upload two FITS images

Enter cluster center

Choose inner/outer radii

The system:

Splits cluster vs calibration stars

Computes calibrated magnitudes

Generates a CMD plot

Outputs calibration parameters


## The Project Structure:
```
project/
│
├── app.py                     # Main Flask application
├── templates/                 # HTML templates
│   ├── index.html
│   ├── object_calibration.html
│   └── star_cluster_calibration.html
│
├── static/                    # Saved plots and assets
│
├── uploads/                   # User FITS uploads
│
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```



