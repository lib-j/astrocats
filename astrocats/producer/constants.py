"""
"""
from collections import OrderedDict
from palettable import cubehelix
import sys

TRAVIS_LIMIT = 100

RADIO_SIGMA = 3.0

DIR_OUT = "astrocats/supernovae/output/"
DIR_CACHE = "cache/"
DIR_JSON = "json/"
DIR_HTML = "html/"

TOOLS_LIST = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"

GOOGLE_PING_URL = "http://www.google.com/webmasters/tools/ping?sitemap=https%3A%2F%2Fsne.space%2Fsitemap.xml"

SNE_LINK_DIR = "https://sne.space/sne/"

DEF_COLORS = cubehelix.perceptual_rainbow_16.hex_colors[:14]

COLUMN_KEYS = [
    "check", "name", "alias", "discoverdate", "maxdate", "maxappmag",
    "maxabsmag", "host", "ra", "dec", "hostra", "hostdec", "hostoffsetang",
    "hostoffsetdist", "instruments", "redshift", "velocity", "lumdist",
    "claimedtype", "ebv", "photolink", "spectralink", "radiolink", "xraylink",
    "references", "download", "responsive"
]

EVENT_IGNORE_KEY = ["download"]

HEADER = [
    "", "Name", "Aliases", "Disc. Date", "Max Date",
    r"<em>m</em><sub>max</sub>", r"<em>M</em><sub>max</sub>", "Host Name",
    "R.A.", "Dec.", "Host R.A.", "Host Dec.", "Host Offset (\")",
    "Host Offset (kpc)", "Instruments/Bands", r"<em>z</em>",
    r"<em>v</em><sub>&#9737;</sub> (km/s)", r"<em>d</em><sub>L</sub> (Mpc)",
    "Type", "E(B-V)", "Phot.", "Spec.", "Radio", "X-ray", "References", "Data",
    ""
]

EVENT_PAGE_HEADER = [
    "", "Name", "Aliases", "Discovery Date", "Maximum Date [band]",
    r"<em>m</em><sub>max</sub> [band]", r"<em>M</em><sub>max</sub> [band]",
    "Host Name", "R.A.", "Dec.", "Host R.A.", "Host Dec.", "Host Offset (\")",
    "Host Offset (kpc)", "Instruments/Bands", r"<em>z</em>",
    r"<em>v</em><sub>&#9737;</sub> (km/s)", r"<em>d</em><sub>L</sub> (Mpc)",
    "Claimed Type", "E(B-V)", "Photometry", "Spectra", "Radio", "X-ray",
    "References", "Download", ""
]

DEF_TITLES = [
    "", "Name (IAU name preferred)", "Aliases",
    "Discovey Date (year-month-day)", "Date of Maximum (year-month-day)",
    "Maximum apparent AB magnitude", "Maximum absolute AB magnitude",
    "Host Name", "Supernova J2000 Right Ascension (h:m:s)",
    "Supernova J2000 Declination (d:m:s)",
    "Host J2000 Right Ascension (h:m:s)", "Host J2000 Declination (d:m:s)",
    "Host Offset (Arcseconds)", "Host Offset (kpc)",
    "List of Instruments and Bands", "Redshift",
    "Heliocentric velocity (km/s)", "Luminosity distance (Mpc)",
    "Claimed Type", "Milky Way Reddening", "Photometry", "pectra", "Radio",
    "X-rays", "Bibcodes of references with most data on event",
    "Download and edit data", ""
]

SNE_PAGES = [["# name", "aliases", "max apparent mag", "max mag date",
              "claimed type", "redshift", "redshift kind", "ra", "dec",
              "# of photometric obs.", "URL"]]


with open('astrocats/supernovae/html/sitemap-template.xml', 'r') as f:
    SITEMAP_TEMPLATE = f.read()

if len(COLUMN_KEYS) != len(HEADER):
    raise (ValueError('Header not same length as key list.'))
    sys.exit(0)

if len(COLUMN_KEYS) != len(EVENT_PAGE_HEADER):
    raise (ValueError('Event page HEADER not same length as key list.'))
    sys.exit(0)

HEADER = OrderedDict(list(zip(COLUMN_KEYS, HEADER)))
EVENT_PAGE_HEADER = OrderedDict(list(zip(COLUMN_KEYS, EVENT_PAGE_HEADER)))
DEF_TITLES = OrderedDict(list(zip(COLUMN_KEYS, DEF_TITLES)))
