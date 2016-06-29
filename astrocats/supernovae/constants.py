"""Supernovae specific Constant variables.
"""
import os
from astropy import constants as const
from astropy import units as un


class FILENAME:
    PATH_BASE = os.path.abspath(os.path.dirname(__file__))
    PATH_INPUT = os.path.join(PATH_BASE, 'input', '')
    PATH_OUTPUT = os.path.join(PATH_BASE, 'output', '')
    # critical datafiles
    REPOS = os.path.join(PATH_INPUT, 'repos.json')
    TASK_LIST = os.path.join(PATH_INPUT, 'tasks.json')
    # auxiliary datafiles
    TYPE_SYNONYMS = os.path.join(PATH_INPUT, 'type-synonyms.json')
    SOURCE_SYNONYMS = os.path.join(PATH_INPUT, 'source-synonyms.json')
    NON_SNE_TYPES = os.path.join(PATH_INPUT, 'non-sne-types.json')
    NON_SNE_PREFIXES = os.path.join(PATH_INPUT, 'non-sne-prefixes.json')
    BIBERRORS = os.path.join(PATH_INPUT, 'biberrors.json')

    BIBAUTHORS = os.path.join(PATH_OUTPUT, 'cache', 'bibauthors.json')
    EXTINCT = os.path.join(PATH_OUTPUT, 'cache', 'extinctions.json')


CLIGHT = const.c.cgs.value
KM = (1.0 * un.km).cgs.value

PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic',
              'photometric', 'host', 'cluster', '']

REPR_BETTER_QUANTITY = {
    'redshift',
    'ebv',
    'velocity',
    'lumdist',
    'discoverdate',
    'maxdate'
}

MAX_BANDS = [
    ['B', 'b', 'g'],  # B-like bands first
    ['V', 'G'],       # if not, V-like bands
    ['R', 'r']        # if not, R-like bands
]
