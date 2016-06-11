#!/usr/local/bin/python3.5

import csv
import os
import re
import urllib
import requests
import calendar
import sys
import json
import codecs
import resource
import argparse
import warnings
from datetime import timedelta, datetime
from glob import glob
from html import unescape
from cdecimal import Decimal
from astroquery.vizier import Vizier
from copy import deepcopy
from astropy.time import Time as astrotime
from collections import OrderedDict
from math import log10, floor, isnan, ceil
from bs4 import BeautifulSoup, Tag, NavigableString
from string import ascii_letters

from .. utils import get_sig_digits, pretty_num, pbar, pbar_strings, is_number, round_sig, tprint, \
    repo_file_list
from . funcs import add_event, add_photometry, add_quantity, add_source, add_spectrum, \
    archived_task, convert_aq_output, \
    derive_and_sanitize, delete_old_event_files, \
    do_task, event_exists, get_aliases, \
    get_bibauthor_dict, get_extinctions_dict, \
    get_max_light, \
    get_preferred_name, has_task, \
    jd_to_mjd, journal_events, load_event_from_file, \
    load_cached_url, make_date_string, merge_duplicates, \
    set_preferred_names, \
    clean_snname, uniq_cdl, utf8, write_all_events
from scripts import PATH
from . constants import KM, CLIGHT, ACKN_CFA, TRAVIS_QUERY_LIMIT


def import_main():
    """
    """
    args = load_args()
    current_task = ''
    # eventnames = []
    events = OrderedDict()
    warnings.filterwarnings('ignore', r'Warning: converting a masked element to nan.')

    tasks = OrderedDict([
        ('deleteoldevents', {'nicename': 'Deleting old events',          'update': False}),
        ('internal',        {'nicename': '%pre metadata and photometry', 'update': False}),
        ('radio',           {'nicename': '%pre radio data',              'update': False}),
        ('xray',            {'nicename': '%pre X-ray data',              'update': False}),
        ('simbad',          {'nicename': '%pre SIMBAD',                  'update': False}),
        ('vizier',          {'nicename': '%pre VizieR',                  'update': False}),
        ('donations',       {'nicename': '%pre donations',               'update': False}),
        ('pessto-dr1',      {'nicename': '%pre PESSTO DR1',              'update': False}),
        ('scp',             {'nicename': '%pre SCP',                     'update': False}),
        ('ascii',           {'nicename': '%pre ASCII',                   'update': False}),
        ('cccp',            {'nicename': '%pre CCCP',                    'update': False, 'archived': True}),
        ('suspect',         {'nicename': '%pre SUSPECT',                 'update': False}),
        ('cfa',             {'nicename': '%pre CfA archive photometry',  'update': False}),
        ('ucb',             {'nicename': '%pre UCB photometry',          'update': False, 'archived': True}),
        ('sdss',            {'nicename': '%pre SDSS photometry',         'update': False}),
        ('csp',             {'nicename': '%pre CSP photometry',          'update': False}),
        ('itep',            {'nicename': '%pre ITEP',                    'update': False}),
        ('asiago',          {'nicename': '%pre Asiago metadata',         'update': False}),
        ('tns',             {'nicename': '%pre TNS metadata',            'update': True,  'archived': True}),
        # ('rochester',       {'nicename': '%pre Latest Supernovae',       'update': True,  'archived': False}),
        ('lennarz',         {'nicename': '%pre Lennarz',                 'update': False}),
        ('fermi',           {'nicename': '%pre Fermi',                   'update': False}),
        ('gaia',            {'nicename': '%pre GAIA',                    'update': True,  'archived': False}),
        ('ogle',            {'nicename': '%pre OGLE',                    'update': True,  'archived': False}),
        ('snls',            {'nicename': '%pre SNLS',                    'update': False}),
        ('psthreepi',       {'nicename': '%pre Pan-STARRS 3π',           'update': True,  'archived': False}),
        ('psmds',           {'nicename': '%pre Pan-STARRS MDS',          'update': False}),
        ('crts',            {'nicename': '%pre CRTS',                    'update': True,  'archived': False}),
        ('snhunt',          {'nicename': '%pre SNhunt',                  'update': True,  'archived': False}),
        ('nedd',            {'nicename': '%pre NED-D',                   'update': False}),
        ('cpcs',            {'nicename': '%pre CPCS',                    'update': True,  'archived': False}),
        ('ptf',             {'nicename': '%pre PTF',                     'update': False, 'archived': False}),
        ('des',             {'nicename': '%pre DES',                     'update': False, 'archived': False}),
        ('asassn',          {'nicename': '%pre ASASSN',                  'update': True}),
        ('asiagospectra',   {'nicename': '%pre Asiago spectra',          'update': True}),
        ('wiserepspectra',  {'nicename': '%pre WISeREP spectra',         'update': False}),
        ('cfaspectra',      {'nicename': '%pre CfA archive spectra',     'update': False}),
        ('snlsspectra',     {'nicename': '%pre SNLS spectra',            'update': False}),
        ('cspspectra',      {'nicename': '%pre CSP spectra',             'update': False}),
        ('ucbspectra',      {'nicename': '%pre UCB spectra',             'update': True,  'archived': True}),
        ('suspectspectra',  {'nicename': '%pre SUSPECT spectra',         'update': False}),
        ('snfspectra',      {'nicename': '%pre SNH spectra',             'update': False}),
        ('superfitspectra', {'nicename': '%pre Superfit spectra',        'update': False}),
        ('mergeduplicates', {'nicename': 'Merging duplicates',           'update': False}),
        ('setprefnames',    {'nicename': 'Setting preferred names',      'update': False}),
        ('writeevents',     {'nicename': 'Writing events',               'update': True})
    ])

    for task in tasks:
        if do_task(tasks, args, task, 'deleteoldevents'):
            # Delete `current_task` here and wrap deletion in `pbar_strings` ??
            current_task = 'Deleting old events'
            delete_old_event_files()

        # Import data provided directly to OSC, in standard format
        if do_task(tasks, args, task, 'internal'):
            events = do_internal(events, args, tasks)

        if do_task(tasks, args, task, 'radio'):
            events = do_external_radio(events, args, tasks)

        if do_task(tasks, args, task, 'xray'):
            events = do_external_xray(events, args, tasks)

        # if do_task(tasks, args, task, 'simbad'):
        #     events = do_simbad(events, args, tasks)

        # Import primary data sources from Vizier
        if do_task(tasks, args, task, 'vizier'):
            Vizier.ROW_LIMIT = -1

            # 2012ApJS..200...12H
            result = Vizier.get_catalogs('J/ApJS/200/12/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            oldname = ''
            for row in pbar(table, current_task):
                name = row['SN']
                if is_number(name[:4]):
                    name = 'SN' + name
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012ApJS..200...12H')
                add_quantity(events, name, 'alias', name, source)
                if '[' not in row['Gal']:
                    add_quantity(events, name, 'host', row['Gal'].replace('_', ' '), source)
                add_quantity(events, name, 'redshift', str(row['z']), source, kind='heliocentric')
                add_quantity(events, name, 'redshift', str(row['zCMB']), source, kind='cmb')
                add_quantity(events, name, 'ebv', str(row['E_B-V_']), source, error=str(row['e_E_B-V_']) if row['e_E_B-V_'] else '')
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)

            # 2012ApJ...746...85S
            result = Vizier.get_catalogs('J/ApJ/746/85/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            oldname = ''
            for row in pbar(table, current_task):
                name = row['Name'].replace('SCP', 'SCP-')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012ApJ...746...85S')
                add_quantity(events, name, 'alias', name, source)
                if row['f_Name']:
                    add_quantity(events, name, 'claimedtype', 'Ia', source)
                if row['z']:
                    add_quantity(events, name, 'redshift', str(row['z']), source, kind='spectroscopic')
                else:
                    add_quantity(events, name, 'redshift', str(row['zCl']), source, kind='cluster')
                add_quantity(events, name, 'ebv', str(row['E_B-V_']), source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)

            result = Vizier.get_catalogs('J/ApJ/746/85/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            oldname = ''
            for row in pbar(table, current_task):
                name = row['Name'].replace('SCP', 'SCP-')
                flux = Decimal(float(row['Flux']))
                if flux <= 0.0:
                    continue
                err = Decimal(float(row['e_Flux']))
                zp = Decimal(float(row['Zero']))
                sig = get_sig_digits(str(row['Flux']))+1
                magnitude = pretty_num(zp-Decimal(2.5)*(flux.log10()), sig=sig)
                e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig=sig)
                if float(e_magnitude) > 5.0:
                    continue
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012ApJ...746...85S')
                add_quantity(events, name, 'alias', name, source)
                add_photometry(
                    events, name, time=str(row['MJD']), band=row['Filter'], instrument=row['Inst'],
                    magnitude=magnitude, e_magnitude=e_magnitude, source=source)

            # 2004ApJ...602..571B
            result = Vizier.get_catalogs('J/ApJ/602/571/table8')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            oldname = ''
            for row in pbar(table, current_task):
                name = 'SN'+row['SN']
                flux = Decimal(float(row['Flux']))
                if flux <= 0.0:
                    continue
                err = Decimal(float(row['e_Flux']))
                sig = get_sig_digits(str(row['Flux']))+1
                magnitude = pretty_num(Decimal(25.0)-Decimal(2.5)*(flux.log10()), sig=sig)
                e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig=sig)
                if float(e_magnitude) > 5.0:
                    continue
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2004ApJ...602..571B')
                add_quantity(events, name, 'alias', name, source)
                band = row['Filt']
                system = ''
                telescope = ''
                if band in ['R', 'I']:
                    system = 'Cousins'
                if band == 'Z':
                    telescope = 'Subaru'
                add_photometry(
                    events, name, time=str(row['MJD']), band=band, system=system, telescope=telescope,
                    magnitude=magnitude, e_magnitude=e_magnitude, source=source)

            # 2014MNRAS.444.3258M
            result = Vizier.get_catalogs('J/MNRAS/444/3258/SNe')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            oldname = ''
            for row in pbar(table, current_task):
                name = row['SN']
                if name == oldname:
                    continue
                oldname = name
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014MNRAS.444.3258M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'redshift', str(row['z']), source, kind='heliocentric', error=str(row['e_z']))
                add_quantity(events, name, 'ra', str(row['_RA']), source, unit='floatdegrees')
                add_quantity(events, name, 'dec', str(row['_DE']), source, unit='floatdegrees')
            journal_events(tasks, args, events)

            # 2014MNRAS.438.1391P
            result = Vizier.get_catalogs('J/MNRAS/438/1391/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014MNRAS.438.1391P')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'redshift', str(row['zh']), source, kind='heliocentric')
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
            journal_events(tasks, args, events)

            # 2012ApJ...749...18B
            result = Vizier.get_catalogs('J/ApJ/749/18/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = row['Name'].replace(' ', '')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012ApJ...749...18B')
                add_quantity(events, name, 'alias', name, source)
                mjd = str(astrotime(2450000.+row['JD'], format='jd').mjd)
                band = row['Filt']
                magnitude = str(row['mag'])
                e_magnitude = str(row['e_mag'])
                e_magnitude = '' if e_magnitude == '--' else e_magnitude
                upperlimit = True if row['l_mag'] == '>' else False
                add_photometry(
                    events, name, time=mjd, band=band, magnitude=magnitude, e_magnitude=e_magnitude, instrument='UVOT',
                    source=source, upperlimit=upperlimit, telescope='Swift', system='Swift')
            journal_events(tasks, args, events)

            # 2010A&A...523A...7G
            result = Vizier.get_catalogs('J/A+A/523/A7/table9')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = 'SNLS-' + row['SNLS']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2010A&A...523A...7G')
                add_quantity(events, name, 'alias', name, source)
                astrot = astrotime(2450000.+row['Date1'], format='jd').datetime
                add_quantity(events, name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
                add_quantity(events, name, 'ebv', str(row['E_B-V_']), source)
                add_quantity(events, name, 'redshift', str(row['z']), source, kind='heliocentric')
                type_str = row['Type'].replace('*', '?').replace('SN', '').replace('(pec)', ' P').replace('Ia? P?', 'Ia P?')
                add_quantity(
                    events, name, 'claimedtype', type_str, source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
            journal_events(tasks, args, events)

            # 2004A&A...415..863G
            result = Vizier.get_catalogs('J/A+A/415/863/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = 'SN' + row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2004A&A...415..863G')
                add_quantity(events, name, 'alias', name, source)
                datesplit = row['Date'].split('-')
                date_str = make_date_string(datesplit[0], datesplit[1].lstrip('0'), datesplit[2].lstrip('0'))
                add_quantity(events, name, 'discoverdate', date_str, source)
                add_quantity(events, name, 'host', 'Abell ' + str(row['Abell']), source)
                add_quantity(events, name, 'claimedtype', row['Type'], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
                if row['zSN']:
                    add_quantity(events, name, 'redshift', str(row['zSN']), source, kind='spectroscopic')
                else:
                    add_quantity(events, name, 'redshift', str(row['zCl']), source, kind='cluster')
            journal_events(tasks, args, events)

            # 2008AJ....136.2306H
            result = Vizier.get_catalogs('J/AJ/136/2306/sources')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = 'SDSS-II ' + str(row['SNID'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2008AJ....136.2306H')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'claimedtype', row['SpType'].replace('SN.', '').strip(':'), source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)

            # 2010ApJ...708..661D
            result = Vizier.get_catalogs('J/ApJ/708/661/sn')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = row['SN']
                if not name:
                    name = 'SDSS-II ' + str(row['SDSS-II'])
                else:
                    name = 'SN' + name
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2010ApJ...708..661D')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'alias', 'SDSS-II ' + str(row['SDSS-II']), source)
                add_quantity(events, name, 'claimedtype', 'II P', source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)

            result = Vizier.get_catalogs('J/ApJ/708/661/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                if row['f_SN'] == 'a':
                    name = 'SDSS-II ' + str(row['SN'])
                else:
                    name = 'SN' + row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2010ApJ...708..661D')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'redshift', str(row['z']), source, error=str(row['e_z']))
            journal_events(tasks, args, events)

            # 2014ApJ...795...44R
            result = Vizier.get_catalogs('J/ApJ/795/44/ps1_snIa')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014ApJ...795...44R')
                add_quantity(events, name, 'alias', name, source)
                astrot = astrotime(row['tdisc'], format='mjd').datetime
                add_quantity(events, name, 'discoverdate',  make_date_string(astrot.year, astrot.month, astrot.day), source)
                add_quantity(events, name, 'redshift', str(row['z']), source, error=str(row['e_z']), kind='heliocentric')
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
                add_quantity(events, name, 'claimedtype', 'Ia', source)

            result = Vizier.get_catalogs('J/ApJ/795/44/table6')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                name = row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014ApJ...795...44R')
                add_quantity(events, name, 'alias', name, source)
                if row['mag'] != '--':
                    add_photometry(
                        events, name, time=str(row['MJD']), band=row['Filt'], magnitude=str(row['mag']),
                        e_magnitude=str(row['e_mag']), source=source, system='AB', telescope='PS1', instrument='PS1')
            journal_events(tasks, args, events)

            # 1990A&AS...82..145C
            result = Vizier.get_catalogs('II/189/mag')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)

            with open(os.path.join(PATH.REPO_EXTERNAL, 'II_189_refs.csv')) as f:
                tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                ii189bibdict = {}
                ii189refdict = {}
                for r, row in enumerate(tsvin):
                    if row[0] != '0':
                        ii189bibdict[r+1] = row[1]
                    else:
                        ii189refdict[r+1] = row[2]

            for row in pbar(table, current_task):
                if row['band'][0] == '(':
                    continue
                name = 'SN' + row['SN']
                name = add_event(tasks, args, events, name)
                source = ''
                secsource = add_source(events, name, bibcode='1990A&AS...82..145C', secondary=True)
                mjd = str(jd_to_mjd(Decimal(row['JD'])))
                mag = str(row['m'])
                band = row['band'].strip("'")
                if row['r_m'] in ii189bibdict:
                    source = add_source(events, name, bibcode=ii189bibdict[row['r_m']])
                else:
                    source = add_source(events, name, refname=ii189refdict[row['r_m']])
                add_quantity(events, name, 'alias', name, source)

                add_photometry(events, name, time=mjd, band=band, magnitude=mag, source=uniq_cdl([source, secsource]))
            journal_events(tasks, args, events)

            # 2014yCat.7272....0G
            result = Vizier.get_catalogs('VII/272/snrs')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)

            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = ''
                if row['Names']:
                    names = row['Names'].split(',')
                    for nam in names:
                        if nam.strip()[:2] == 'SN':
                            name = nam.strip()
                            if is_number(name[2:]):
                                name = name + 'A'
                    if not name:
                        for nam in names:
                            if nam.strip('()') == nam:
                                name = nam.strip()
                                break
                if not name:
                    name = row['SNR'].strip()

                name = add_event(tasks, args, events, name)
                source = (add_source(events, name, bibcode='2014BASI...42...47G') + ',' +
                          add_source(events, name, refname='Galactic SNRs',
                                     url='https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html'))
                add_quantity(events, name, 'alias', name, source)

                add_quantity(events, name, 'alias', row['SNR'].strip(), source)
                add_quantity(events, name, 'alias', 'MWSNR '+row['SNR'].strip('G '), source)

                if row['Names']:
                    names = row['Names'].split(',')
                    for nam in names:
                        add_quantity(events, name, 'alias', nam.replace('Vela (XYZ)', 'Vela').strip('()').strip(), source)
                        if nam.strip()[:2] == 'SN':
                            add_quantity(events, name, 'discoverdate', nam.strip()[2:], source)

                add_quantity(events, name, 'host', 'Milky Way', source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
            journal_events(tasks, args, events)

            # 2014MNRAS.442..844F
            result = Vizier.get_catalogs('J/MNRAS/442/844/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = 'SN' + row['SN']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014MNRAS.442..844F')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'redshift', str(row['zhost']), source, kind='host')
                add_quantity(events, name, 'ebv', str(row['E_B-V_']), source)
            journal_events(tasks, args, events)

            result = Vizier.get_catalogs('J/MNRAS/442/844/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = 'SN' + str(row['SN'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014MNRAS.442..844F')
                add_quantity(events, name, 'alias', name, source)
                for band in ['B', 'V', 'R', 'I']:
                    bandtag = band + 'mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=row['MJD'], band=band, magnitude=row[bandtag],
                            e_magnitude=row['e_' + bandtag], source=source, telescope='KAIT', instrument='KAIT')
            journal_events(tasks, args, events)

            # 2012MNRAS.425.1789S
            result = Vizier.get_catalogs('J/MNRAS/425/1789/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = ''.join(row['SimbadName'].split(' '))
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012MNRAS.425.1789S')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'alias', 'SN' + row['SN'], source)
                add_quantity(events, name, 'host', row['Gal'], source)
                red_str = str(round_sig(float(row['cz'])*KM/CLIGHT, sig=get_sig_digits(str(row['cz']))))
                if is_number(row['cz']):
                    add_quantity(events, name, 'redshift', red_str, source, kind='heliocentric')
                add_quantity(events, name, 'ebv', str(row['E_B-V_']), source)
            journal_events(tasks, args, events)

            # 2015ApJS..219...13W
            result = Vizier.get_catalogs('J/ApJS/219/13/table3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = u'LSQ' + str(row['LSQ'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015ApJS..219...13W')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source)
                add_quantity(events, name, 'dec', row['DEJ2000'], source)
                add_quantity(events, name, 'redshift', row['z'], source, error=row['e_z'], kind='heliocentric')
                add_quantity(events, name, 'ebv', row['E_B-V_'], source)
                add_quantity(events, name, 'claimedtype', 'Ia', source)
            result = Vizier.get_catalogs('J/ApJS/219/13/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = 'LSQ' + row['LSQ']
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015ApJS..219...13W')
                add_quantity(events, name, 'alias', name, source)
                add_photometry(
                    events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument='QUEST', telescope='ESO Schmidt',
                    observatory='La Silla', band=row['Filt'],
                    magnitude=row['mag'], e_magnitude=row['e_mag'], system='Swope', source=source)
            journal_events(tasks, args, events)

            # 2012Natur.491..228C
            result = Vizier.get_catalogs('J/other/Nat/491.228/tablef1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'SN2213-1745'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2012Natur.491..228C')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'claimedtype', 'SLSN-R', source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                for band in ['g', 'r', 'i']:
                    bandtag = band + '_mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=row['MJD' + band + '_'], band=band + "'", magnitude=row[bandtag],
                            e_magnitude=row['e_' + bandtag], source=source)

            result = Vizier.get_catalogs('J/other/Nat/491.228/tablef2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'SN1000+0216'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2012Natur.491..228C')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'claimedtype', 'SLSN-II?', source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                for band in ['g', 'r', 'i']:
                    bandtag = band + '_mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=row['MJD' + band + '_'], band=band + "'", magnitude=row[bandtag],
                            e_magnitude=row['e_' + bandtag], source=source)
            journal_events(tasks, args, events)

            # 2011Natur.474..484Q
            result = Vizier.get_catalogs('J/other/Nat/474.484/tables1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2011Natur.474..484Q')
                add_quantity(events, name, 'alias', name, source)
                add_photometry(
                    events, name, time=row['MJD'], band=row['Filt'], telescope=row['Tel'],
                    magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
            journal_events(tasks, args, events)

            # 2011ApJ...736..159G
            result = Vizier.get_catalogs('J/ApJ/736/159/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'PTF10vdl'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2011ApJ...736..159G')
            add_quantity(events, name, 'alias', name, source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                add_photometry(
                    events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'],
                    telescope=row['Tel'], magnitude=row['mag'],
                    e_magnitude=row['e_mag'] if is_number(row['e_mag']) else '',
                    upperlimit=(not is_number(row['e_mag'])), source=source)
            journal_events(tasks, args, events)

            # 2012ApJ...760L..33B
            result = Vizier.get_catalogs('J/ApJ/760/L33/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'PTF12gzk'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2012ApJ...760L..33B')
            add_quantity(events, name, 'alias', name, source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                # Fixing a typo in VizieR table
                if str(row['JD']) == '2455151.456':
                    row['JD'] = '2456151.456'
                add_photometry(
                    events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'],
                    telescope=row['Inst'], magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
            journal_events(tasks, args, events)

            # 2013ApJ...769...39S
            result = Vizier.get_catalogs('J/ApJ/769/39/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'PS1-12sk'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2013ApJ...769...39S')
            add_quantity(events, name, 'alias', name, source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                instrument = ''
                telescope = ''
                if row['Inst'] == 'RATCam':
                    instrument = row['Inst']
                else:
                    telescope = row['Inst']
                add_photometry(
                    events, name, time=row['MJD'], band=row['Filt'], telescope=telescope,
                    instrument=instrument, magnitude=row['mag'],
                    e_magnitude=row['e_mag'] if not row['l_mag'] else '',
                    upperlimit=(row['l_mag'] == '>'), source=source)
            journal_events(tasks, args, events)

            # 2009MNRAS.394.2266P
            # Note: Instrument info available via links in VizieR, can't auto-parse just yet.
            name = 'SN2005cs'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2009MNRAS.394.2266P')
            add_quantity(events, name, 'alias', name, source)
            result = Vizier.get_catalogs('J/MNRAS/394/2266/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                for band in ['U', 'B', 'V', 'R', 'I']:
                    bandtag = band + 'mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=band, magnitude=row[bandtag],
                            e_magnitude=(row['e_' + bandtag] if row['l_' + bandtag] != '>' else ''),
                            source=source, upperlimit=(row['l_' + bandtag] == '>'))
                if 'zmag' in row and is_number(row['zmag']) and not isnan(float(row['zmag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band='z', magnitude=row['zmag'],
                                   e_magnitude=row['e_zmag'], source=source)

            result = Vizier.get_catalogs('J/MNRAS/394/2266/table3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                for band in ['B', 'V', 'R']:
                    bandtag = band + 'mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=band, magnitude=row[bandtag],
                            e_magnitude=(row['e_' + bandtag] if row['l_' + bandtag] != '>' else ''),
                            source=source, upperlimit=(row['l_' + bandtag] == '>'))

            result = Vizier.get_catalogs('J/MNRAS/394/2266/table4')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                for band in ['J', 'H', 'K']:
                    bandtag = band + 'mag'
                    if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                        add_photometry(
                            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=band,
                            magnitude=row[bandtag],
                            e_magnitude=row['e_' + bandtag], source=source)
            journal_events(tasks, args, events)

            # 2013AJ....145...99A
            result = Vizier.get_catalogs('J/AJ/145/99/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            name = 'SN2003ie'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2013AJ....145...99A')
            add_quantity(events, name, 'alias', name, source)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=row['MJD'], band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'] if not row['l_Bmag'] else '',
                                   upperlimit=(row['l_Bmag'] == '>'), source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=row['MJD'], band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'] if is_number(row['e_Vmag']) else '',
                                   upperlimit=(not is_number(row['e_Vmag'])), source=source)
                if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                    add_photometry(events, name, time=row['MJD'], band='R', magnitude=row['Rmag'],
                                   e_magnitude=row['e_Rmag'] if not row['l_Rmag'] else '',
                                   upperlimit=(row['l_Rmag'] == '>'), source=source)
                if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                    add_photometry(events, name, time=row['MJD'], band='I', magnitude=row['Imag'],
                                   e_magnitude=row['e_Imag'], source=source)
            journal_events(tasks, args, events)

            # 2011ApJ...729..143C
            name = 'SN2008am'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2011ApJ...729..143C')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/ApJ/729/143/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                add_photometry(
                    events, name, time=row['MJD'], band='ROTSE', telescope='ROTSE', magnitude=row['mag'],
                    e_magnitude=row['e_mag'] if not row['l_mag'] else '', upperlimit=(row['l_mag'] == '<'), source=source)

            result = Vizier.get_catalogs('J/ApJ/729/143/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                if 'Jmag' in row and is_number(row['Jmag']) and not isnan(float(row['Jmag'])):
                    add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='J', magnitude=row['Jmag'],
                                   e_magnitude=row['e_Jmag'], source=source)
                if 'Hmag' in row and is_number(row['Hmag']) and not isnan(float(row['Hmag'])):
                    add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='H', magnitude=row['Hmag'],
                                   e_magnitude=row['e_Hmag'], source=source)
                if 'Ksmag' in row and is_number(row['Ksmag']) and not isnan(float(row['Ksmag'])):
                    add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='Ks', magnitude=row['Ksmag'],
                                   e_magnitude=row['e_Ksmag'], source=source)

            result = Vizier.get_catalogs('J/ApJ/729/143/table4')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                add_photometry(
                    events, name, time=row['MJD'], band=row['Filt'], telescope='P60',
                    magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)

            result = Vizier.get_catalogs('J/ApJ/729/143/table5')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                add_photometry(
                    events, name, time=row['MJD'], band=row['Filt'], instrument='UVOT', telescope='Swift',
                    magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
            journal_events(tasks, args, events)

            # 2011ApJ...728...14P
            name = 'SN2009bb'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2011ApJ...728...14P')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/ApJ/728/14/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(
                        events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                        band='B', magnitude=row['Bmag'], e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(
                        events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                        band='V', magnitude=row['Vmag'], e_magnitude=row['e_Vmag'], source=source)
                if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                    add_photometry(
                        events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                        band='R', magnitude=row['Rmag'], e_magnitude=row['e_Rmag'], source=source)
                if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                    add_photometry(
                        events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                        band='I', magnitude=row['Imag'], e_magnitude=row['e_Imag'], source=source)

            result = Vizier.get_catalogs('J/ApJ/728/14/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                if 'u_mag' in row and is_number(row['u_mag']) and not isnan(float(row['u_mag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='u', magnitude=row['u_mag'],
                                   e_magnitude=row['e_u_mag'], source=source)
                if 'g_mag' in row and is_number(row['g_mag']) and not isnan(float(row['g_mag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='g', magnitude=row['g_mag'],
                                   e_magnitude=row['e_g_mag'], source=source)
                if 'r_mag' in row and is_number(row['r_mag']) and not isnan(float(row['r_mag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='r', magnitude=row['r_mag'],
                                   e_magnitude=row['e_r_mag'], source=source)
                if 'i_mag' in row and is_number(row['i_mag']) and not isnan(float(row['i_mag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='i', magnitude=row['i_mag'],
                                   e_magnitude=row['e_i_mag'], source=source)
                if 'z_mag' in row and is_number(row['z_mag']) and not isnan(float(row['z_mag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='z', magnitude=row['z_mag'],
                                   e_magnitude=row['e_z_mag'], source=source)

            result = Vizier.get_catalogs('J/ApJ/728/14/table3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                if 'Ymag' in row and is_number(row['Ymag']) and not isnan(float(row['Ymag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='Y', magnitude=row['Ymag'],
                                   e_magnitude=row['e_Ymag'], source=source)
                if 'Jmag' in row and is_number(row['Jmag']) and not isnan(float(row['Jmag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='J', magnitude=row['Jmag'],
                                   e_magnitude=row['e_Jmag'], source=source)
                if 'Hmag' in row and is_number(row['Hmag']) and not isnan(float(row['Hmag'])):
                    add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='H', magnitude=row['Hmag'],
                                   e_magnitude=row['e_Hmag'], source=source)
            journal_events(tasks, args, events)

            # 2011PAZh...37..837T
            name = 'SN2009nr'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2011PAZh...37..837T')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/PAZh/37/837/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = str(jd_to_mjd(Decimal(row['JD']) + 2455000))
                if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='U', magnitude=row['Umag'],
                                   e_magnitude=row['e_Umag'], source=source)
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'], source=source)
                if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='R', magnitude=row['Rmag'],
                                   e_magnitude=row['e_Rmag'], source=source)
                if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='I', magnitude=row['Imag'],
                                   e_magnitude=row['e_Imag'], source=source)
            journal_events(tasks, args, events)

            # 2013MNRAS.433.1871B
            name = 'SN2012aw'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2013MNRAS.433.1871B')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/MNRAS/433/1871/table3a')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = str(jd_to_mjd(Decimal(row['JD']) + 2456000))
                if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='U', magnitude=row['Umag'],
                                   e_magnitude=row['e_Umag'], source=source)
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'], source=source)
                if 'Rcmag' in row and is_number(row['Rcmag']) and not isnan(float(row['Rcmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='Rc', magnitude=row['Rcmag'],
                                   e_magnitude=row['e_Rcmag'], source=source)
                if 'Icmag' in row and is_number(row['Icmag']) and not isnan(float(row['Icmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='Ic', magnitude=row['Icmag'],
                                   e_magnitude=row['e_Icmag'], source=source)

            result = Vizier.get_catalogs('J/MNRAS/433/1871/table3b')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = str(jd_to_mjd(Decimal(row['JD']) + 2456000))
                if 'gmag' in row and is_number(row['gmag']) and not isnan(float(row['gmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='g', magnitude=row['gmag'],
                                   e_magnitude=row['e_gmag'], source=source)
                if 'rmag' in row and is_number(row['rmag']) and not isnan(float(row['rmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='r', magnitude=row['rmag'],
                                   e_magnitude=row['e_rmag'], source=source)
                if 'imag' in row and is_number(row['imag']) and not isnan(float(row['imag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='i', magnitude=row['imag'],
                                   e_magnitude=row['e_imag'], source=source)
                if 'zmag' in row and is_number(row['zmag']) and not isnan(float(row['zmag'])):
                    add_photometry(events, name, time=mjd, telescope=row['Tel'], band='z', magnitude=row['zmag'],
                                   e_magnitude=row['e_zmag'], source=source)
            journal_events(tasks, args, events)

            # 2014AJ....148....1Z
            name = 'SN2012fr'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2014AJ....148....1Z')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/AJ/148/1/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = row['MJD']
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'], source=source)
                if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='R', magnitude=row['Rmag'],
                                   e_magnitude=row['e_Rmag'], source=source)
                if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='I', magnitude=row['Imag'],
                                   e_magnitude=row['e_Imag'], source=source)

            result = Vizier.get_catalogs('J/AJ/148/1/table3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = row['MJD']
                if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='U', magnitude=row['Umag'],
                                   e_magnitude=row['e_Umag'], source=source)
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'], source=source)
                if 'UVW1' in row and is_number(row['UVW1']) and not isnan(float(row['UVW1'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='W1', magnitude=row['UVW1'],
                                   e_magnitude=row['e_UVW1'], source=source)
                if 'UVW2' in row and is_number(row['UVW2']) and not isnan(float(row['UVW2'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='W2', magnitude=row['UVW2'],
                                   e_magnitude=row['e_UVW2'], source=source)
                if 'UVM2' in row and is_number(row['UVM2']) and not isnan(float(row['UVM2'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='M2', magnitude=row['UVM2'],
                                   e_magnitude=row['e_UVM2'], source=source)

            result = Vizier.get_catalogs('J/AJ/148/1/table5')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = row['MJD']
                if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', band='B', magnitude=row['Bmag'],
                                   e_magnitude=row['e_Bmag'], source=source)
                if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', band='V', magnitude=row['Vmag'],
                                   e_magnitude=row['e_Vmag'], source=source)
                if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', band='R', magnitude=row['Rmag'],
                                   e_magnitude=row['e_Rmag'], source=source)
                if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                    add_photometry(events, name, time=mjd, telescope='LJT', band='I', magnitude=row['Imag'],
                                   e_magnitude=row['e_Imag'], source=source)
            journal_events(tasks, args, events)

            # 2015ApJ...805...74B
            name = 'SN2014J'
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2014AJ....148....1Z')
            add_quantity(events, name, 'alias', name, source)

            result = Vizier.get_catalogs('J/ApJ/805/74/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                mjd = row['MJD']
                if 'mag' in row and is_number(row['mag']) and not isnan(float(row['mag'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band=row['Filt'], magnitude=row['mag'],
                                   e_magnitude=row['e_mag'], source=source)
                elif 'maglim' in row and is_number(row['maglim']) and not isnan(float(row['maglim'])):
                    add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band=row['Filt'], magnitude=row['maglim'],
                                   upperlimit=True, source=source)
            journal_events(tasks, args, events)

            # 2011ApJ...741...97D
            result = Vizier.get_catalogs('J/ApJ/741/97/table2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['SN'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2011ApJ...741...97D')
                add_quantity(events, name, 'alias', name, source)
                add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'], magnitude=row['mag'],
                               e_magnitude=row['e_mag'] if is_number(row['e_mag']) else '', upperlimit=(not is_number(row['e_mag'])), source=source)
            journal_events(tasks, args, events)

            # 2015MNRAS.448.1206M
            # Note: Photometry from two SN can also be added from this source.
            result = Vizier.get_catalogs('J/MNRAS/448/1206/table3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'redshift', row['zsp'], source, kind='spectroscopic')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', 'Ia', source)
            result = Vizier.get_catalogs('J/MNRAS/448/1206/table4')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'redshift', row['zph'], source, error=row['e_zph'], kind='photometric')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', 'Ia?', source)
            result = Vizier.get_catalogs('J/MNRAS/448/1206/table5')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'redshift', row['zsp'], source, kind='spectroscopic')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', row['Type'], source)
            result = Vizier.get_catalogs('J/MNRAS/448/1206/table6')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', row['Type'], source)
            result = Vizier.get_catalogs('J/MNRAS/448/1206/tablea2')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', row['Typesoft']+'?', source)
                add_quantity(events, name, 'claimedtype', row['Typepsnid']+'?', source)
            result = Vizier.get_catalogs('J/MNRAS/448/1206/tablea3')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = str(row['Name'])
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2015MNRAS.448.1206M')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'discoverdate', '20' + name[4:6], source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
                add_quantity(events, name, 'maxband', 'r', source)
                add_quantity(events, name, 'claimedtype', 'Candidate', source)
            journal_events(tasks, args, events)

            # 2012AJ....143..126B
            result = Vizier.get_catalogs('J/AJ/143/126/table4')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                if not row['Wcl'] or row['Wcl'] == 'N':
                    continue
                row = convert_aq_output(row)
                name = str(row['SN']).replace(' ', '')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012AJ....143..126B')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'claimedtype', 'Ia-' + row['Wcl'], source)
            journal_events(tasks, args, events)

            # 2015ApJS..220....9F
            for viztab in ['1', '2']:
                result = Vizier.get_catalogs('J/ApJS/220/9/table' + viztab)
                table = result[list(result.keys())[0]]
                table.convert_bytestring_to_unicode(python3_only=True)
                for row in pbar(table, current_task):
                    row = convert_aq_output(row)
                    name = add_event(tasks, args, events, row['SN'])
                    source = add_source(events, name, bibcode='2015ApJS..220....9F')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'claimedtype', row['Type'], source)
                    add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                    add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                    if '?' not in row['Host']:
                        add_quantity(events, name, 'host', row['Host'].replace('_', ' '), source)
                    kind = ''
                    if 'Host' in row['n_z']:
                        kind = 'host'
                    elif 'Spectrum' in row['n_z']:
                        kind = 'spectroscopic'
                    add_quantity(events, name, 'redshift', row['z'], source, error=row['e_z'], kind=kind)

            result = Vizier.get_catalogs('J/ApJS/220/9/table8')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = add_event(tasks, args, events, row['SN'])
                source = add_source(events, name, bibcode='2015ApJS..220....9F')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'claimedtype', row['Type'], source)
                add_photometry(
                    events, name, time=row['MJD'], band=row['Band'], magnitude=row['mag'],
                    e_magnitude=row['e_mag'], telescope=row['Tel'], source=source)
            journal_events(tasks, args, events)

            result = Vizier.get_catalogs('J/ApJ/673/999/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = add_event(tasks, args, events, 'SN'+row['SN'])
                source = add_source(events, name, bibcode='2008ApJ...673..999P')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'ra', row['RAJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row['DEJ2000'], source, unit='floatdegrees')
                add_quantity(events, name, 'redshift', row['z'], source, kind='host')
                add_quantity(events, name, 'hostra', row['RAGdeg'], source, unit='floatdegrees')
                add_quantity(events, name, 'hostdec', row['DEGdeg'], source, unit='floatdegrees')
                add_quantity(events, name, 'claimedtype', row['Type'].strip(':'), source)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'donations'):
            # Nicholl 04-01-16 donation
            with open(os.path.join(PATH.REPO_EXTERNAL, 'Nicholl-04-01-16/bibcodes.json'), 'r') as f:
                bcs = json.loads(f.read())

            for datafile in sorted(glob(os.path.join(PATH.REPO_EXTERNAL, 'Nicholl-04-01-16/*.txt')), key=lambda s: s.lower()):
                name = os.path.basename(datafile).split('_')[0]
                name = add_event(tasks, args, events, name)
                bibcode = ''
                for bc in bcs:
                    if name in bcs[bc]:
                        bibcode = bc
                if not bibcode:
                    raise(ValueError('Bibcode not found!'))
                source = add_source(events, name, bibcode=bibcode)
                add_quantity(events, name, 'alias', name, source)
                with open(datafile, 'r') as f:
                    tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                    for r, rrow in enumerate(tsvin):
                        row = list(filter(None, rrow))
                        if not row:
                            continue
                        if row[0][0] == '#' and row[0] != '#MJD':
                            continue
                        if row[0] == '#MJD':
                            bands = [x for x in row[1:] if x and 'err' not in x]
                            continue
                        mjd = row[0]
                        if not is_number(mjd):
                            continue
                        for v, val in enumerate(row[1::2]):
                            upperlimit = ''
                            if '>' in val:
                                upperlimit = True
                            mag = val.strip('>')
                            if not is_number(mag) or isnan(float(mag)) or float(mag) > 90.0:
                                continue
                            err = ''
                            if is_number(row[2*v+2]) and not isnan(float(row[2*v+2])):
                                err = row[2*v+2]
                            add_photometry(
                                events, name, time=mjd, band=bands[v], magnitude=mag,
                                e_magnitude=err, upperlimit=upperlimit, source=source)
            journal_events(tasks, args, events)

            # Maggi 04-11-16 donation (MC SNRs)
            with open(os.path.join(PATH.REPO_EXTERNAL, 'Maggi-04-11-16/LMCSNRs_OpenSNe.csv')) as f:
                tsvin = csv.reader(f, delimiter=',')
                for row in tsvin:
                    name = 'MCSNR ' + row[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2016A&A...585A.162M')
                    add_quantity(events, name, 'alias', name, source)
                    if row[1] != 'noname':
                        add_quantity(events, name, 'alias', row[1], source)
                    add_quantity(events, name, 'ra', row[2], source)
                    add_quantity(events, name, 'dec', row[3], source)
                    add_quantity(events, name, 'host', 'LMC', source)
                    if row[4] == '1':
                        add_quantity(events, name, 'claimedtype', 'Ia', source)
                    elif row[4] == '2':
                        add_quantity(events, name, 'claimedtype', 'CC', source)
            with open(os.path.join(PATH.REPO_EXTERNAL, 'Maggi-04-11-16/SMCSNRs_OpenSNe.csv')) as f:
                tsvin = csv.reader(f, delimiter=',')
                for row in tsvin:
                    name = 'MCSNR ' + row[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, refname='Pierre Maggi')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'alias', row[1], source)
                    add_quantity(events, name, 'alias', row[2], source)
                    add_quantity(events, name, 'ra', row[3], source)
                    add_quantity(events, name, 'dec', row[4], source)
                    add_quantity(events, name, 'host', 'SMC', source)
            journal_events(tasks, args, events)

            # Galbany 04-18-16 donation
            folders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/')))[1]
            bibcode = '2016AJ....151...33G'
            for folder in folders:
                infofiles = glob(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/') + folder + '/*.info')
                photfiles = glob(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/') + folder + '/*.out*')

                zhel = ''
                zcmb = ''
                zerr = ''
                for path in infofiles:
                    with open(path, 'r') as f:
                        lines = f.read().splitlines()
                        for line in lines:
                            splitline = line.split(':')
                            field = splitline[0].strip().lower()
                            value = splitline[1].strip()
                            if field == 'name':
                                name = value[:6].upper() + (value[6].upper() if len(value) == 7 else value[6:])
                                name = add_event(tasks, args, events, name)
                                source = add_source(events, name, bibcode=bibcode)
                                add_quantity(events, name, 'alias', name, source)
                            elif field == 'type':
                                claimedtype = value.replace('SN', '')
                                add_quantity(events, name, 'claimedtype', claimedtype, source)
                            elif field == 'zhel':
                                zhel = value
                            elif field == 'redshift_error':
                                zerr = value
                            elif field == 'zcmb':
                                zcmb = value
                            elif field == 'ra':
                                add_quantity(events, name, 'ra', value, source, unit='floatdegrees')
                            elif field == 'dec':
                                add_quantity(events, name, 'dec', value, source, unit='floatdegrees')
                            elif field == 'host':
                                add_quantity(events, name, 'host', value.replace('- ', '-').replace('G ', 'G'), source)
                            elif field == 'e(b-v)_mw':
                                add_quantity(events, name, 'ebv', value, source)

                add_quantity(events, name, 'redshift', zhel, source, error=zerr, kind='heliocentric')
                add_quantity(events, name, 'redshift', zcmb, source, error=zerr, kind='cmb')

                for path in photfiles:
                    with open(path, 'r') as f:
                        band = ''
                        lines = f.read().splitlines()
                        for li, line in enumerate(lines):
                            if li in [0, 2, 3]:
                                continue
                            if li == 1:
                                band = line.split(':')[-1].strip()
                            else:
                                cols = list(filter(None, line.split()))
                                if not cols:
                                    continue
                                add_photometry(events, name, time=cols[0], magnitude=cols[1], e_magnitude=cols[2],
                                               band=band, system=cols[3], telescope=cols[4], source=source)
            journal_events(tasks, args, events)

            # Brown 05-14-16
            files = glob(os.path.join(PATH.REPO_EXTERNAL, 'brown-05-14-16/*.dat'))
            for fi in pbar(files, current_task):
                name = os.path.basename(fi).split('_')[0]
                name = add_event(tasks, args, events, name)
                source = add_source(
                    events, name, refname='Swift Supernovae', bibcode='2014Ap&SS.354...89B',
                    url='http://people.physics.tamu.edu/pbrown/SwiftSN/swift_sn.html')
                add_quantity(events, name, 'alias', name, source)
                with open(fi, 'r') as f:
                    lines = f.read().splitlines()
                    for line in lines:
                        if not line or line[0] == '#':
                            continue
                        cols = list(filter(None, line.split()))
                        band = cols[0]
                        mjd = cols[1]
                        # Skip lower limit entries for now
                        if cols[2] == 'NULL' and cols[6] == 'NULL':
                            continue
                        isupp = cols[2] == 'NULL' and cols[6] != 'NULL'
                        mag = cols[2] if not isupp else cols[4]
                        e_mag = cols[3] if not isupp else ''
                        upp = '' if not isupp else True
                        add_photometry(events, name, time=mjd, magnitude=mag, e_magnitude=e_mag,
                                       upperlimit=upp, band=band, source=source,
                                       telescope='Swift', instrument='UVOT', system='Vega')
            journal_events(tasks, args, events)

            # Nicholl 05-03-16
            files = glob(os.path.join(PATH.REPO_EXTERNAL, 'nicholl-05-03-16/*.txt'))
            name = add_event(tasks, args, events, 'SN2015bn')
            source = add_source(events, name, bibcode='2016arXiv160304748N')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'alias', 'PS15ae', source)
            for fi in pbar(files, current_task):
                telescope = os.path.basename(fi).split('_')[1]
                with open(fi, 'r') as f:
                    lines = f.read().splitlines()
                    for li, line in enumerate(lines):
                        if not line or (line[0] == '#' and li != 0):
                            continue
                        cols = list(filter(None, line.split()))
                        if not cols:
                            continue
                        if li == 0:
                            bands = cols[1:]
                            continue

                        mjd = cols[0]
                        for ci, col in enumerate(cols[1::2]):
                            if not is_number(col):
                                continue

                            emag = cols[2*ci+2]
                            upp = ''
                            if not is_number(emag):
                                emag = ''
                                upp = True
                            add_photometry(events, name, time=mjd, magnitude=col, e_magnitude=emag,
                                           upperlimit=upp, band=bands[ci], source=source,
                                           telescope=telescope, instrument='UVOT' if telescope == 'Swift' else '')
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'pessto-dr1'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'PESSTO_MPHOT.csv'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',')
                for ri, row in enumerate(tsvin):
                    if ri == 0:
                        bands = [x.split('_')[0] for x in row[3::2]]
                        systems = [x.split('_')[1].capitalize().replace('Ab', 'AB') for x in row[3::2]]
                        continue
                    name = row[1]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2015A&A...579A..40S')
                    add_quantity(events, name, 'alias', name, source)
                    for hi, ci in enumerate(range(3, len(row)-1, 2)):
                        if not row[ci]:
                            continue
                        add_photometry(
                            events, name, time=row[2], magnitude=row[ci], e_magnitude=row[ci+1],
                            band=bands[hi], system=systems[hi], telescope='Swift' if systems[hi] == 'Swift' else '',
                            source=source)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'scp'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'SCP09.csv'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',')
                for ri, row in enumerate(pbar(tsvin, current_task)):
                    if ri == 0:
                        continue
                    name = row[0].replace('SCP', 'SCP-')
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, refname='Supernova Cosmology Project', url='http://supernova.lbl.gov/2009ClusterSurvey/')
                    add_quantity(events, name, 'alias', name, source)
                    if row[1]:
                        add_quantity(events, name, 'alias', row[1], source)
                    if row[2]:
                        add_quantity(events, name, 'redshift', row[2], source, kind='spectroscopic' if row[3] == 'sn' else 'host')
                    if row[4]:
                        add_quantity(events, name, 'redshift', row[2], source, kind='cluster')
                    if row[6]:
                        claimedtype = row[6].replace('SN ', '')
                        kind = ('spectroscopic/light curve' if 'a' in row[7] and 'c' in row[7] else
                                'spectroscopic' if 'a' in row[7] else 'light curve' if 'c' in row[7] else '')
                        if claimedtype != '?':
                            add_quantity(events, name, 'claimedtype', claimedtype, source, kind=kind)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ascii'):
            # 2006ApJ...645..841N
            with open(os.path.join(PATH.REPO_EXTERNAL, '2006ApJ...645..841N-table3.csv'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',')
                for ri, row in enumerate(pbar(tsvin, current_task)):
                    name = 'SNLS-' + row[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2006ApJ...645..841N')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'redshift', row[1], source, kind='spectroscopic')
                    astrot = astrotime(float(row[4]) + 2450000., format='jd').datetime
                    add_quantity(events, name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
            journal_events(tasks, args, events)

            # Anderson 2014
            for datafile in pbar_strings(glob(os.path.join(PATH.REPO_EXTERNAL, 'SNII_anderson2014/*.dat')), desc=current_task):
                basename = os.path.basename(datafile)
                if not is_number(basename[:2]):
                    continue
                if basename == '0210_V.dat':
                    name = 'SN0210'
                else:
                    name = ('SN20' if int(basename[:2]) < 50 else 'SN19') + basename.split('_')[0]
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2014ApJ...786...67A')
                add_quantity(events, name, 'alias', name, source)

                if name in ['SN1999ca', 'SN2003dq', 'SN2008aw']:
                    system = 'Swope'
                else:
                    system = 'Landolt'

                with open(datafile, 'r') as f:
                    tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
                    for row in tsvin:
                        if not row[0]:
                            continue
                        add_photometry(events, name, time=str(jd_to_mjd(Decimal(row[0]))), band='V', magnitude=row[1], e_magnitude=row[2], system=system, source=source)
            journal_events(tasks, args, events)

            # stromlo
            stromlobands = ['B', 'V', 'R', 'I', 'VM', 'RM']
            with open(os.path.join(PATH.REPO_EXTERNAL, 'J_A+A_415_863-1/photometry.csv'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',')
                for row in pbar(tsvin, current_task):
                    name = row[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2004A&A...415..863G')
                    add_quantity(events, name, 'alias', name, source)
                    mjd = str(jd_to_mjd(Decimal(row[1])))
                    for ri, ci in enumerate(range(2, len(row), 3)):
                        if not row[ci]:
                            continue
                        band = stromlobands[ri]
                        upperlimit = True if (not row[ci+1] and row[ci+2]) else False
                        e_upper_magnitude = str(abs(Decimal(row[ci+1]))) if row[ci+1] else ''
                        e_lower_magnitude = str(abs(Decimal(row[ci+2]))) if row[ci+2] else ''
                        add_photometry(
                            events, name, time=mjd, band=band, magnitude=row[ci],
                            e_upper_magnitude=e_upper_magnitude, e_lower_magnitude=e_lower_magnitude,
                            upperlimit=upperlimit, telescope='MSSSO 1.3m' if band in ['VM', 'RM'] else 'CTIO',
                            instrument='MaCHO' if band in ['VM', 'RM'] else '', source=source)
            journal_events(tasks, args, events)

            # 2015MNRAS.449..451W
            with open(os.path.join(PATH.REPO_EXTERNAL, '2015MNRAS.449..451W.dat'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                for r, row in enumerate(pbar(data, current_task)):
                    if r == 0:
                        continue
                    namesplit = row[0].split('/')
                    name = namesplit[-1]
                    if name.startswith('SN'):
                        name = name.replace(' ', '')
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2015MNRAS.449..451W')
                    add_quantity(events, name, 'alias', name, source)
                    if len(namesplit) > 1:
                        add_quantity(events, name, 'alias', namesplit[0], source)
                    add_quantity(events, name, 'claimedtype', row[1], source)
                    add_photometry(events, name, time=row[2], band=row[4], magnitude=row[3], source=source)
            journal_events(tasks, args, events)

            # 2016MNRAS.459.1039T
            with open(os.path.join(PATH.REPO_EXTERNAL, '2016MNRAS.459.1039T.tsv'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                name = add_event(tasks, args, events, 'LSQ13zm')
                source = add_source(events, name, bibcode='2016MNRAS.459.1039T')
                add_quantity(events, name, 'alias', name, source)
                for r, row in enumerate(pbar(data, current_task)):
                    if row[0][0] == '#':
                        bands = [x.replace('(err)', '') for x in row[3:-1]]
                        continue
                    mjd = row[1]
                    mags = [re.sub(r'\([^)]*\)', '', x) for x in row[3:-1]]
                    upps = [True if '>' in x else '' for x in mags]
                    mags = [x.replace('>', '') for x in mags]
                    errs = [x[x.find('(')+1:x.find(')')] if '(' in x else '' for x in row[3:-1]]
                    for mi, mag in enumerate(mags):
                        if not is_number(mag):
                            continue
                        add_photometry(
                            events, name, time=mjd, band=bands[mi], magnitude=mag, e_magnitude=errs[mi],
                            instrument=row[-1], upperlimit=upps[mi], source=source)
            journal_events(tasks, args, events)

            # 2015ApJ...804...28G
            with open(os.path.join(PATH.REPO_EXTERNAL, '2015ApJ...804...28G.tsv'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                name = add_event(tasks, args, events, 'PS1-13arp')
                source = add_source(events, name, bibcode='2015ApJ...804...28G')
                add_quantity(events, name, 'alias', name, source)
                for r, row in enumerate(pbar(data, current_task)):
                    if r == 0:
                        continue
                    mjd = row[1]
                    mag = row[3]
                    upp = True if '<' in mag else ''
                    mag = mag.replace('<', '')
                    err = row[4] if is_number(row[4]) else ''
                    ins = row[5]
                    add_photometry(
                        events, name, time=mjd, band=row[0], magnitude=mag, e_magnitude=err,
                        instrument=ins, upperlimit=upp, source=source)
            journal_events(tasks, args, events)

            # 2016ApJ...819...35A
            with open(os.path.join(PATH.REPO_EXTERNAL, '2016ApJ...819...35A.tsv'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                for r, row in enumerate(pbar(data, current_task)):
                    if row[0][0] == '#':
                        continue
                    name = add_event(tasks, args, events, row[0])
                    source = add_source(events, name, bibcode='2016ApJ...819...35A')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'ra', row[1], source)
                    add_quantity(events, name, 'dec', row[2], source)
                    add_quantity(events, name, 'redshift', row[3], source)
                    add_quantity(
                        events, name, 'discoverdate',
                        datetime.strptime(row[4], '%Y %b %d').isoformat().split('T')[0].replace('-', '/'), source)
            journal_events(tasks, args, events)

            # 2014ApJ...784..105W
            with open(os.path.join(PATH.REPO_EXTERNAL, '2014ApJ...784..105W.tsv'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                for r, row in enumerate(pbar(data, current_task)):
                    if row[0][0] == '#':
                        continue
                    name = add_event(tasks, args, events, row[0])
                    source = add_source(events, name, bibcode='2014ApJ...784..105W')
                    add_quantity(events, name, 'alias', name, source)
                    mjd = row[1]
                    band = row[2]
                    mag = row[3]
                    err = row[4]
                    add_photometry(
                        events, name, time=mjd, band=row[2], magnitude=mag, e_magnitude=err,
                        instrument='WHIRC', telescope='WIYN 3.5 m', observatory='NOAO',
                        system='WHIRC', source=source)
            journal_events(tasks, args, events)

            # 2012MNRAS.425.1007B
            with open(os.path.join(PATH.REPO_EXTERNAL, '2012MNRAS.425.1007B.tsv'), 'r') as f:
                data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
                for r, row in enumerate(pbar(data, current_task)):
                    if row[0][0] == '#':
                        bands = row[2:]
                        continue
                    name = add_event(tasks, args, events, row[0])
                    source = add_source(events, name, bibcode='2012MNRAS.425.1007B')
                    add_quantity(events, name, 'alias', name, source)
                    mjd = row[1]
                    mags = [x.split('±')[0].strip() for x in row[2:]]
                    errs = [x.split('±')[1].strip() if '±' in x else '' for x in row[2:]]
                    if row[0] == 'PTF09dlc':
                        ins = 'HAWK-I'
                        tel = 'VLT 8.1m'
                        obs = 'ESO'
                    else:
                        ins = 'NIRI'
                        tel = 'Gemini North 8.2m'
                        obs = 'Gemini'

                    for mi, mag in enumerate(mags):
                        if not is_number(mag):
                            continue
                        add_photometry(
                            events, name, time=mjd, band=bands[mi], magnitude=mag, e_magnitude=errs[mi],
                            instrument=ins, telescope=tel, observatory=obs,
                            system='Natural', source=source)
            journal_events(tasks, args, events)

        # CCCP
        if do_task(tasks, args, task, 'cccp'):
            cccpbands = ['B', 'V', 'R', 'I']
            for datafile in sorted(glob(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/apj407397*.txt')), key=lambda s: s.lower()):
                with open(datafile, 'r') as f:
                    tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                    for r, row in enumerate(tsvin):
                        if r == 0:
                            continue
                        elif r == 1:
                            name = 'SN' + row[0].split('SN ')[-1]
                            name = add_event(tasks, args, events, name)
                            source = add_source(events, name, bibcode='2012ApJ...744...10K')
                            add_quantity(events, name, 'alias', name, source)
                        elif r >= 5:
                            mjd = str(Decimal(row[0]) + 53000)
                            for b, band in enumerate(cccpbands):
                                if row[2*b + 1]:
                                    add_photometry(
                                        events, name, time=mjd, band=band, magnitude=row[2*b + 1].strip('>'),
                                        e_magnitude=row[2*b + 2], upperlimit=(not row[2*b + 2]), source=source)

            if archived_task(tasks, args, 'cccp'):
                with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'r') as f:
                    html = f.read()
            else:
                session = requests.Session()
                response = session.get('https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
                html = response.text
                with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'w') as f:
                    f.write(html)

            soup = BeautifulSoup(html, 'html5lib')
            links = soup.body.findAll("a")
            for link in pbar(links, current_task):
                if 'sc_sn' in link['href']:
                    name = add_event(tasks, args, events, link.text.replace(' ', ''))
                    source = add_source(events, name, refname='CCCP', url='https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
                    add_quantity(events, name, 'alias', name, source)

                    if archived_task(tasks, args, 'cccp'):
                        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link['href'].split('/')[-1], 'r') as f:
                            html2 = f.read()
                    else:
                        response2 = session.get('https://webhome.weizmann.ac.il/home/iair/' + link['href'])
                        html2 = response2.text
                        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link['href'].split('/')[-1], 'w') as f:
                            f.write(html2)

                    soup2 = BeautifulSoup(html2, 'html5lib')
                    links2 = soup2.body.findAll("a")
                    for link2 in links2:
                        if '.txt' in link2['href'] and '_' in link2['href']:
                            band = link2['href'].split('_')[1].split('.')[0].upper()
                            if archived_task(tasks, args, 'cccp'):
                                fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link2['href'].split('/')[-1]
                                if not os.path.isfile(fname):
                                    continue
                                with open(fname, 'r') as f:
                                    html3 = f.read()
                            else:
                                response3 = session.get('https://webhome.weizmann.ac.il/home/iair/cccp/' + link2['href'])
                                if response3.status_code == 404:
                                    continue
                                html3 = response3.text
                                with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link2['href'].split('/')[-1], 'w') as f:
                                    f.write(html3)
                            table = [[str(Decimal(y.strip())).rstrip('0') for y in x.split(',')] for x in list(filter(None, html3.split('\n')))]
                            for row in table:
                                add_photometry(events, name, time=str(Decimal(row[0]) + 53000), band=band, magnitude=row[1], e_magnitude=row[2], source=source)
            journal_events(tasks, args, events)

        # Suspect catalog
        if do_task(tasks, args, task, 'suspect'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'suspectreferences.csv'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',', skipinitialspace=True)
                suspectrefdict = {}
                for row in tsvin:
                    suspectrefdict[row[0]] = row[1]

            for datafile in pbar_strings(glob(os.path.join(PATH.REPO_EXTERNAL, 'SUSPECT/*.html')), desc=current_task):
                basename = os.path.basename(datafile)
                basesplit = basename.split('-')
                name = basesplit[1]
                name = add_event(tasks, args, events, name)
                if name.startswith('SN') and is_number(name[2:]):
                    name = name + 'A'
                band = basesplit[3].split('.')[0]
                ei = int(basesplit[2])
                bandlink = 'file://' + os.path.abspath(datafile)
                bandresp = urllib.request.urlopen(bandlink)
                bandsoup = BeautifulSoup(bandresp, 'html5lib')
                bandtable = bandsoup.find('table')

                names = bandsoup.body.findAll(text=re.compile('Name'))
                reference = ''
                for link in bandsoup.body.findAll('a'):
                    if 'adsabs' in link['href']:
                        reference = str(link).replace('"', "'")

                bibcode = unescape(suspectrefdict[reference])
                source = add_source(events, name, bibcode=bibcode)

                secondaryreference = 'SUSPECT'
                secondaryrefurl = 'https://www.nhn.ou.edu/~suspect/'
                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)

                if ei == 1:
                    year = re.findall(r'\d+', name)[0]
                    add_quantity(events, name, 'discoverdate', year, secondarysource)
                    add_quantity(events, name, 'host', names[1].split(':')[1].strip(), secondarysource)

                    redshifts = bandsoup.body.findAll(text=re.compile('Redshift'))
                    if redshifts:
                        add_quantity(events, name, 'redshift', redshifts[0].split(':')[1].strip(), secondarysource, kind='heliocentric')
                    # hvels = bandsoup.body.findAll(text=re.compile('Heliocentric Velocity'))
                    # if hvels:
                    #    add_quantity(events, name, 'velocity', hvels[0].split(':')[1].strip().split(' ')[0],
                    #        secondarysource, kind='heliocentric')
                    types = bandsoup.body.findAll(text=re.compile('Type'))

                    add_quantity(events, name, 'claimedtype', types[0].split(':')[1].strip().split(' ')[0], secondarysource)

                for r, row in enumerate(bandtable.findAll('tr')):
                    if r == 0:
                        continue
                    col = row.findAll('td')
                    mjd = str(jd_to_mjd(Decimal(col[0].contents[0])))
                    mag = col[3].contents[0]
                    if mag.isspace():
                        mag = ''
                    else:
                        mag = str(mag)
                    e_magnitude = col[4].contents[0]
                    if e_magnitude.isspace():
                        e_magnitude = ''
                    else:
                        e_magnitude = str(e_magnitude)
                    add_photometry(events, name, time=mjd, band=band, magnitude=mag, e_magnitude=e_magnitude, source=secondarysource + ',' + source)
            journal_events(tasks, args, events)

        # CfA data
        if do_task(tasks, args, task, 'cfa'):
            for fname in pbar_strings(glob(os.path.join(PATH.REPO_EXTERNAL, 'cfa-input/*.dat')), desc=current_task):
                f = open(fname, 'r')
                tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
                csv_data = []
                for r, row in enumerate(tsvin):
                    new = []
                    for item in row:
                        new.extend(item.split('\t'))
                    csv_data.append(new)

                for r, row in enumerate(csv_data):
                    for c, col in enumerate(row):
                        csv_data[r][c] = col.strip()
                    csv_data[r] = [_f for _f in csv_data[r] if _f]

                eventname = os.path.basename(os.path.splitext(fname)[0])

                eventparts = eventname.split('_')

                name = clean_snname(eventparts[0])
                name = add_event(tasks, args, events, name)
                secondaryname = 'CfA Supernova Archive'
                secondaryurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
                secondarysource = add_source(events, name, refname=secondaryname, url=secondaryurl, secondary=True, acknowledgment=ACKN_CFA)
                add_quantity(events, name, 'alias', name, secondarysource)

                year = re.findall(r'\d+', name)[0]
                add_quantity(events, name, 'discoverdate', year, secondarysource)

                eventbands = list(eventparts[1])

                tu = 'MJD'
                jdoffset = Decimal(0.)
                for rc, row in enumerate(csv_data):
                    if len(row) > 0 and row[0][0] == "#":
                        if len(row[0]) > 2 and row[0][:3] == '#JD':
                            tu = 'JD'
                            rowparts = row[0].split('-')
                            jdoffset = Decimal(rowparts[1])
                        elif len(row[0]) > 6 and row[0][:7] == '#Julian':
                            tu = 'JD'
                            jdoffset = Decimal(0.)
                        elif len(row) > 1 and row[1].lower() == 'photometry':
                            for ci, col in enumerate(row[2:]):
                                if col[0] == "(":
                                    refstr = ' '.join(row[2+ci:])
                                    refstr = refstr.replace('(', '').replace(')', '')
                                    bibcode = unescape(refstr)
                                    source = add_source(events, name, bibcode=bibcode)

                        elif len(row) > 1 and row[1] == 'HJD':
                            tu = 'HJD'

                        continue
                    elif len(row) > 0:
                        mjd = row[0]
                        for v, val in enumerate(row):
                            if v == 0:
                                if tu == 'JD':
                                    mjd = str(jd_to_mjd(Decimal(val) + jdoffset))
                                    tuout = 'MJD'
                                elif tu == 'HJD':
                                    mjd = str(jd_to_mjd(Decimal(val)))
                                    tuout = 'MJD'
                                else:
                                    mjd = val
                                    tuout = tu
                            elif v % 2 != 0:
                                if float(row[v]) < 90.0:
                                    add_photometry(events, name, u_time=tuout, time=mjd, band=eventbands[(v-1)//2], magnitude=row[v], e_magnitude=row[v+1], source=secondarysource + ',' + source)
                f.close()

            # Hicken 2012
            f = open(os.path.join(PATH.REPO_EXTERNAL, 'hicken-2012-standard.dat'), 'r')
            tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
            for r, row in enumerate(pbar(tsvin, current_task)):
                if r <= 47:
                    continue

                if row[0][:2] != 'sn':
                    name = 'SN' + row[0].strip()
                else:
                    name = row[0].strip()

                name = add_event(tasks, args, events, name)

                source = add_source(events, name, bibcode='2012ApJS..200...12H')
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'claimedtype', 'Ia', source)
                add_photometry(
                    events, name, u_time='MJD', time=row[2].strip(), band=row[1].strip(),
                    magnitude=row[6].strip(), e_magnitude=row[7].strip(), source=source)

            # Bianco 2014
            tsvin = open(os.path.join(PATH.REPO_EXTERNAL, 'bianco-2014-standard.dat'), 'r')
            tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
            for row in pbar(tsvin, current_task):
                name = 'SN' + row[0]
                name = add_event(tasks, args, events, name)

                source = add_source(events, name, bibcode='2014ApJS..213...19B')
                add_quantity(events, name, 'alias', name, source)
                add_photometry(
                    events, name, u_time='MJD', time=row[2], band=row[1], magnitude=row[3],
                    e_magnitude=row[4], telescope=row[5], system='Standard', source=source)
            f.close()
            journal_events(tasks, args, events)

        # New UCB import
        if do_task(tasks, args, task, 'ucb'):
            secondaryreference = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
            secondaryrefurl = 'http://heracles.astro.berkeley.edu/sndb/info'
            secondaryrefbib = '2012MNRAS.425.1789S'

            jsontxt = load_cached_url(
                args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubphot',
                os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
            if not jsontxt:
                continue

            photom = json.loads(jsontxt)
            photom = sorted(photom, key=lambda k: k['ObjName'])
            for phot in pbar(photom, desc=current_task):
                name = phot['ObjName']
                name = add_event(tasks, args, events, name)

                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondaryrefbib, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                sources = [secondarysource]
                if phot['Reference']:
                    sources += [add_source(events, name, bibcode=phot['Reference'])]
                sources = uniq_cdl(sources)

                if phot['Type'] and phot['Type'].strip() != 'NoMatch':
                    for ct in phot['Type'].strip().split(','):
                        add_quantity(events, name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
                if phot['DiscDate']:
                    add_quantity(events, name, 'discoverdate', phot['DiscDate'].replace('-', '/'), sources)
                if phot['HostName']:
                    add_quantity(events, name, 'host', urllib.parse.unquote(phot['HostName']).replace('*', ''), sources)
                filename = phot['Filename'] if phot['Filename'] else ''

                if not filename:
                    raise(ValueError('Filename not found for SNDB phot!'))
                if not phot['PhotID']:
                    raise(ValueError('ID not found for SNDB phot!'))

                filepath = os.path.join(PATH.REPO_EXTERNAL, 'SNDB/') + filename
                if archived_task(tasks, args, 'ucb') and os.path.isfile(filepath):
                    with open(filepath, 'r') as f:
                        phottxt = f.read()
                else:
                    session = requests.Session()
                    response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=dp:' + str(phot['PhotID']))
                    phottxt = response.text
                    with open(filepath, 'w') as f:
                        f.write(phottxt)

                tsvin = csv.reader(phottxt.splitlines(), delimiter=' ', skipinitialspace=True)

                for r, row in enumerate(tsvin):
                    if len(row) > 0 and row[0] == "#":
                        continue
                    mjd = row[0]
                    magnitude = row[1]
                    if magnitude and float(magnitude) > 99.0:
                        continue
                    e_magnitude = row[2]
                    band = row[4]
                    telescope = row[5]
                    add_photometry(
                        events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                        e_magnitude=e_magnitude, source=sources)

            journal_events(tasks, args, events)

        # Import SDSS
        if do_task(tasks, args, task, 'sdss'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/2010ApJ...708..661D.txt'), 'r') as f:
                bibcodes2010 = f.read().split('\n')
            sdssbands = ['u', 'g', 'r', 'i', 'z']
            for fname in pbar_strings(glob(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/*.sum')), desc=current_task):
                f = open(fname, 'r')
                tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

                basename = os.path.basename(fname)
                if basename in bibcodes2010:
                    bibcode = '2010ApJ...708..661D'
                else:
                    bibcode = '2008AJ....136.2306H'

                for r, row in enumerate(tsvin):
                    if r == 0:
                        if row[5] == 'RA:':
                            name = 'SDSS-II ' + row[3]
                        else:
                            name = 'SN' + row[5]
                        name = add_event(tasks, args, events, name)
                        source = add_source(events, name, bibcode=bibcode)
                        add_quantity(events, name, 'alias', name, source)
                        add_quantity(events, name, 'alias', 'SDSS-II ' + row[3], source)

                        if row[5] != 'RA:':
                            year = re.findall(r'\d+', name)[0]
                            add_quantity(events, name, 'discoverdate', year, source)

                        add_quantity(events, name, 'ra', row[-4], source, unit='floatdegrees')
                        add_quantity(events, name, 'dec', row[-2], source, unit='floatdegrees')
                    if r == 1:
                        error = row[4] if float(row[4]) >= 0.0 else ''
                        add_quantity(events, name, 'redshift', row[2], source, error=error, kind='heliocentric')
                    if r >= 19:
                        # Skip bad measurements
                        if int(row[0]) > 1024:
                            continue

                        mjd = row[1]
                        band = sdssbands[int(row[2])]
                        magnitude = row[3]
                        e_magnitude = row[4]
                        telescope = 'SDSS'
                        add_photometry(
                            events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                            e_magnitude=e_magnitude, source=source, system='SDSS')
                f.close()
            journal_events(tasks, args, events)

        # Import GAIA
        if do_task(tasks, args, task, 'gaia'):
            fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/alerts.csv')
            csvtxt = load_cached_url(args, 'http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv', fname)
            if not csvtxt:
                continue
            tsvin = csv.reader(csvtxt.splitlines(), delimiter=',', skipinitialspace=True)
            reference = 'Gaia Photometric Science Alerts'
            refurl = 'http://gsaweb.ast.cam.ac.uk/alerts/alertsindex'
            for ri, row in enumerate(pbar(tsvin, current_task)):
                if ri == 0 or not row:
                    continue
                name = add_event(tasks, args, events, row[0])
                source = add_source(events, name, refname=reference, url=refurl)
                add_quantity(events, name, 'alias', name, source)
                year = '20' + re.findall(r'\d+', row[0])[0]
                add_quantity(events, name, 'discoverdate', year, source)
                add_quantity(events, name, 'ra', row[2], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row[3], source, unit='floatdegrees')
                if row[7] and row[7] != 'unknown':
                    add_quantity(events, name, 'claimedtype', row[7].replace('SNe', '').replace('SN', '').strip(), source)
                elif (True in [x in row[9].upper() for x in ['SN CANDIATE', 'CANDIDATE SN', 'HOSTLESS SN']]):
                    add_quantity(events, name, 'claimedtype', 'Candidate', source)

                if 'aka' in row[9].replace('gakaxy', 'galaxy').lower() and 'AKARI' not in row[9]:
                    commentsplit = (row[9].replace('_', ' ').replace('MLS ', 'MLS').replace('CSS ', 'CSS').
                                    replace('SN iPTF', 'iPTF').replace('SN ', 'SN').replace('AT ', 'AT').split())
                    for csi, cs in enumerate(commentsplit):
                        if 'aka' in cs.lower() and csi < len(commentsplit) - 1:
                            alias = commentsplit[csi+1].strip('(),:.').replace('PSNJ', 'PSN J')
                            if alias[:6] == 'ASASSN' and alias[6] != '-':
                                alias = 'ASASSN-' + alias[6:]
                            add_quantity(events, name, 'alias', alias, source)
                            break

                fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/') + row[0] + '.csv'
                if not args.fullrefresh and archived_task(tasks, args, 'gaia') and os.path.isfile(fname):
                    with open(fname, 'r') as f:
                        csvtxt = f.read()
                else:
                    response = urllib.request.urlopen('http://gsaweb.ast.cam.ac.uk/alerts/alert/' + row[0] + '/lightcurve.csv')
                    with open(fname, 'w') as f:
                        csvtxt = response.read().decode('utf-8')
                        f.write(csvtxt)

                tsvin2 = csv.reader(csvtxt.splitlines())
                for ri2, row2 in enumerate(tsvin2):
                    if ri2 <= 1 or not row2:
                        continue
                    mjd = str(jd_to_mjd(Decimal(row2[1].strip())))
                    magnitude = row2[2].strip()
                    if magnitude == 'null':
                        continue
                    e_magnitude = 0.
                    telescope = 'GAIA'
                    band = 'G'
                    add_photometry(events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude, e_magnitude=e_magnitude, source=source)
                if args.update:
                    journal_events(tasks, args, events)
            journal_events(tasks, args, events)

        # Import CSP
        # VizieR catalogs exist for this: J/AJ/139/519, J/AJ/142/156. Should replace eventually.
        if do_task(tasks, args, task, 'csp'):
            cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
            for fname in pbar_strings(glob(os.path.join(PATH.REPO_EXTERNAL, 'CSP/*.dat')), desc=current_task):
                f = open(fname, 'r')
                tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)

                eventname = os.path.basename(os.path.splitext(fname)[0])

                eventparts = eventname.split('opt+')

                name = clean_snname(eventparts[0])
                name = add_event(tasks, args, events, name)

                reference = 'Carnegie Supernova Project'
                refbib = '2010AJ....139..519C'
                refurl = 'http://csp.obs.carnegiescience.edu/data'
                source = add_source(events, name, bibcode=refbib, refname=reference, url=refurl)
                add_quantity(events, name, 'alias', name, source)

                year = re.findall(r'\d+', name)[0]
                add_quantity(events, name, 'discoverdate', year, source)

                for r, row in enumerate(tsvin):
                    if len(row) > 0 and row[0][0] == "#":
                        if r == 2:
                            add_quantity(events, name, 'redshift', row[0].split(' ')[-1], source, kind='cmb')
                            add_quantity(events, name, 'ra', row[1].split(' ')[-1], source)
                            add_quantity(events, name, 'dec', row[2].split(' ')[-1], source)
                        continue
                    for v, val in enumerate(row):
                        if v == 0:
                            mjd = val
                        elif v % 2 != 0:
                            if float(row[v]) < 90.0:
                                add_photometry(
                                    events, name, time=mjd, observatory='LCO', band=cspbands[(v-1)//2],
                                    system='CSP', magnitude=row[v], e_magnitude=row[v+1], source=source)
                f.close()
            journal_events(tasks, args, events)

        # Import ITEP
        if do_task(tasks, args, task, 'itep'):
            itepbadsources = ['2004ApJ...602..571B']

            needsbib = []
            with open(os.path.join(PATH.REPO_EXTERNAL, 'itep-refs.txt'), 'r') as f:
                refrep = f.read().splitlines()
            refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
            f = open(os.path.join(PATH.REPO_EXTERNAL, 'itep-lc-cat-28dec2015.txt'), 'r')
            tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
            curname = ''
            for r, row in enumerate(pbar(tsvin, current_task)):
                if r <= 1 or len(row) < 7:
                    continue
                name = 'SN' + row[0].strip()
                mjd = str(jd_to_mjd(Decimal(row[1].strip())))
                band = row[2].strip()
                magnitude = row[3].strip()
                e_magnitude = row[4].strip()
                reference = row[6].strip().strip(',')

                if curname != name:
                    curname = name
                    name = add_event(tasks, args, events, name)

                    secondaryreference = 'Sternberg Astronomical Institute Supernova Light Curve Catalogue'
                    secondaryrefurl = 'http://dau.itep.ru/sn/node/72'
                    secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, secondary=True)
                    add_quantity(events, name, 'alias', name, secondarysource)

                    year = re.findall(r'\d+', name)[0]
                    add_quantity(events, name, 'discoverdate', year, secondarysource)
                if reference in refrepf:
                    bibcode = unescape(refrepf[reference])
                    source = add_source(events, name, bibcode=bibcode)
                else:
                    needsbib.append(reference)
                    source = add_source(events, name, refname=reference) if reference else ''

                if bibcode not in itepbadsources:
                    add_photometry(events, name, time=mjd, band=band, magnitude=magnitude,
                                   e_magnitude=e_magnitude, source=secondarysource + ',' + source)
            f.close()

            # Write out references that could use a bibcode
            needsbib = list(OrderedDict.fromkeys(needsbib))
            with open('../itep-needsbib.txt', 'w') as f:
                f.writelines(['%s\n' % i for i in needsbib])
            journal_events(tasks, args, events)

        # Now import the Asiago catalog
        if do_task(tasks, args, task, 'asiago'):
            # response = urllib.request.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
            path = os.path.abspath(os.path.join(PATH.REPO_EXTERNAL, 'asiago-cat.php'))
            response = urllib.request.urlopen('file://' + path)
            html = response.read().decode('utf-8')
            html = html.replace('\r', "")

            soup = BeautifulSoup(html, 'html5lib')
            table = soup.find('table')

            records = []
            for r, row in enumerate(table.findAll('tr')):
                if r == 0:
                    continue

                col = row.findAll('td')
                records.append([utf8(x.renderContents()) for x in col])

            for record in pbar(records, current_task):
                if len(record) > 1 and record[1] != '':
                    name = clean_snname('SN' + record[1]).strip('?')
                    name = add_event(tasks, args, events, name)

                    reference = 'Asiago Supernova Catalogue'
                    refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                    refbib = '1989A&AS...81..421B'
                    source = add_source(
                        events, name, refname=reference, url=refurl, bibcode=refbib, secondary=True)
                    add_quantity(events, name, 'alias', name, source)

                    year = re.findall(r'\d+', name)[0]
                    add_quantity(events, name, 'discoverdate', year, source)

                    hostname = record[2]
                    hostra = record[3]
                    hostdec = record[4]
                    ra = record[5].strip(':')
                    dec = record[6].strip(':')
                    redvel = record[11].strip(':')
                    discoverer = record[19]

                    datestring = year

                    monthday = record[18]
                    if "*" in monthday:
                        datekey = 'discover'
                    else:
                        datekey = 'max'

                    if monthday.strip() != '':
                        monthstr = ''.join(re.findall('[a-zA-Z]+', monthday))
                        monthstr = str(list(calendar.month_abbr).index(monthstr))
                        datestring = datestring + '/' + monthstr

                        dayarr = re.findall(r'\d+', monthday)
                        if dayarr:
                            daystr = dayarr[0]
                            datestring = datestring + '/' + daystr

                    add_quantity(events, name, datekey + 'date', datestring, source)

                    velocity = ''
                    redshift = ''
                    if redvel != '':
                        if round(float(redvel)) == float(redvel):
                            velocity = int(redvel)
                        else:
                            redshift = float(redvel)
                        redshift = str(redshift)
                        velocity = str(velocity)

                    claimedtype = record[17].replace(':', '').replace('*', '').strip()

                    if (hostname != ''):
                        add_quantity(events, name, 'host', hostname, source)
                    if (claimedtype != ''):
                        add_quantity(events, name, 'claimedtype', claimedtype, source)
                    if (redshift != ''):
                        add_quantity(events, name, 'redshift', redshift, source, kind='host')
                    if (velocity != ''):
                        add_quantity(events, name, 'velocity', velocity, source, kind='host')
                    if (hostra != ''):
                        add_quantity(events, name, 'hostra', hostra, source, unit='nospace')
                    if (hostdec != ''):
                        add_quantity(events, name, 'hostdec', hostdec, source, unit='nospace')
                    if (ra != ''):
                        add_quantity(events, name, 'ra', ra, source, unit='nospace')
                    if (dec != ''):
                        add_quantity(events, name, 'dec', dec, source, unit='nospace')
                    if (discoverer != ''):
                        add_quantity(events, name, 'discoverer', discoverer, source)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'lennarz'):
            Vizier.ROW_LIMIT = -1
            result = Vizier.get_catalogs('J/A+A/538/A120/usc')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)

            bibcode = '2012A&A...538A.120L'
            for row in pbar(table, current_task):
                row = convert_aq_output(row)
                name = 'SN' + row['SN']
                name = add_event(tasks, args, events, name)

                source = add_source(events, name, bibcode=bibcode)
                add_quantity(events, name, 'alias', name, source)

                if row['RAJ2000']:
                    add_quantity(events, name, 'ra', row['RAJ2000'], source)
                if row['DEJ2000']:
                    add_quantity(events, name, 'dec', row['DEJ2000'], source)
                if row['RAG']:
                    add_quantity(events, name, 'hostra', row['RAG'], source)
                if row['DEG']:
                    add_quantity(events, name, 'hostdec', row['DEG'], source)
                if row['Gal']:
                    add_quantity(events, name, 'host', row['Gal'], source)
                if row['Type']:
                    claimedtypes = row['Type'].split('|')
                    for claimedtype in claimedtypes:
                        add_quantity(events, name, 'claimedtype', claimedtype.strip(' -'), source)
                if row['z']:
                    if name not in ['SN1985D', 'SN2004cq']:
                        add_quantity(events, name, 'redshift', row['z'], source, kind='host')
                if row['Dist']:
                    if row['e_Dist']:
                        add_quantity(events, name, 'lumdist', row['Dist'], source, error=row['e_Dist'], kind='host')
                    else:
                        add_quantity(events, name, 'lumdist', row['Dist'], source, kind='host')

                if row['Ddate']:
                    datestring = row['Ddate'].replace('-', '/')

                    add_quantity(events, name, 'discoverdate', datestring, source)

                    if 'photometry' not in events[name]:
                        if 'Dmag' in row and is_number(row['Dmag']) and not isnan(float(row['Dmag'])):
                            datesplit = row['Ddate'].strip().split('-')
                            if len(datesplit) == 3:
                                datestr = row['Ddate'].strip()
                            elif len(datesplit) == 2:
                                datestr = row['Ddate'].strip() + '-01'
                            elif len(datesplit) == 1:
                                datestr = row['Ddate'].strip() + '-01-01'
                            mjd = str(astrotime(datestr).mjd)
                            add_photometry(events, name, time=mjd, band=row['Dband'], magnitude=row['Dmag'], source=source)
                if row['Mdate']:
                    datestring = row['Mdate'].replace('-', '/')

                    add_quantity(events, name, 'maxdate', datestring, source)

                    if 'photometry' not in events[name]:
                        if 'MMag' in row and is_number(row['MMag']) and not isnan(float(row['MMag'])):
                            datesplit = row['Mdate'].strip().split('-')
                            if len(datesplit) == 3:
                                datestr = row['Mdate'].strip()
                            elif len(datesplit) == 2:
                                datestr = row['Mdate'].strip() + '-01'
                            elif len(datesplit) == 1:
                                datestr = row['Mdate'].strip() + '-01-01'
                            mjd = str(astrotime(datestr).mjd)
                            add_photometry(events, name, time=mjd, band=row['Mband'], magnitude=row['Mmag'], source=source)
            f.close()
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'fermi'):
            with open(os.path.join(PATH.REPO_EXTERNAL, '1SC_catalog_v01.asc'), 'r') as f:
                tsvin = csv.reader(f, delimiter=',')
                for ri, row in enumerate(pbar(tsvin, current_task)):
                    if row[0].startswith('#'):
                        if len(row) > 1 and 'UPPER_LIMITS' in row[1]:
                            break
                        continue
                    if 'Classified' not in row[1]:
                        continue
                    name = row[0].replace('SNR', 'G')
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2016ApJS..224....8A')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'alias', row[0].replace('SNR', 'MWSNR'), source)
                    add_quantity(events, name, 'ra', row[2], source, unit='floatdegrees')
                    add_quantity(events, name, 'dec', row[3], source, unit='floatdegrees')
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'tns'):
            session = requests.Session()
            csvtxt = load_cached_url(
                args, 'https://wis-tns.weizmann.ac.il/search?&num_page=1&format=html&sort=desc&order=id&format=csv&page=0',
                os.path.join(PATH.REPO_EXTERNAL, 'TNS/index.csv'))
            if not csvtxt:
                continue
            maxid = csvtxt.splitlines()[1].split(',')[0].strip('"')
            maxpages = ceil(int(maxid)/1000.)

            for page in pbar(range(maxpages), current_task):
                fname = os.path.join(PATH.REPO_EXTERNAL, 'TNS/page-') + str(page).zfill(2) + '.csv'
                if archived_task(tasks, args, 'tns') and os.path.isfile(fname) and page < 7:
                    with open(fname, 'r') as f:
                        csvtxt = f.read()
                else:
                    with open(fname, 'w') as f:
                        session = requests.Session()
                        response = session.get('https://wis-tns.weizmann.ac.il/search?&num_page=1000&format=html&edit[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1&display[hostname]=1&display[host_redshift]=1&display[source_group_name]=1&display[programs_name]=1&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1&display[end_pop_period]=0&display[spectra_count]=1&display[discoverymag]=1&display[discmagfilter]=1&display[discoverydate]=1&display[discoverer]=1&display[sources]=1&display[bibcode]=1&format=csv&page=' + str(page))
                        csvtxt = response.text
                        f.write(csvtxt)

                tsvin = csv.reader(csvtxt.splitlines(), delimiter=',')
                for ri, row in enumerate(pbar(tsvin, current_task, leave=False)):
                    if ri == 0:
                        continue
                    if row[4] and 'SN' not in row[4]:
                        continue
                    name = row[1].replace(' ', '')
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, refname='Transient Name Server', url='https://wis-tns.weizmann.ac.il')
                    add_quantity(events, name, 'alias', name, source)
                    if row[2] and row[2] != '00:00:00.00':
                        add_quantity(events, name, 'ra', row[2], source)
                    if row[3] and row[3] != '+00:00:00.00':
                        add_quantity(events, name, 'dec', row[3], source)
                    if row[4]:
                        add_quantity(events, name, 'claimedtype', row[4].replace('SN', '').strip(), source)
                    if row[5]:
                        add_quantity(events, name, 'redshift', row[5], source, kind='spectroscopic')
                    if row[6]:
                        add_quantity(events, name, 'host', row[6], source)
                    if row[7]:
                        add_quantity(events, name, 'redshift', row[7], source, kind='host')
                    if row[8]:
                        add_quantity(events, name, 'discoverer', row[8], source)
                    # Currently, all events listing all possible observers. TNS bug?
                    # if row[9]:
                    #    observers = row[9].split(',')
                    #    for observer in observers:
                    #        add_quantity(events, name, 'observer', observer.strip(), source)
                    if row[10]:
                        add_quantity(events, name, 'alias', row[10], source)
                    if row[8] and row[14] and row[15] and row[16]:
                        survey = row[8]
                        magnitude = row[14]
                        band = row[15].split('-')[0]
                        mjd = astrotime(row[16]).mjd
                        add_photometry(events, name, time=mjd, magnitude=magnitude, band=band, survey=survey, source=source)
                    if row[16]:
                        date = row[16].split()[0].replace('-', '/')
                        if date != '0000/00/00':
                            date = date.replace('/00', '')
                            time = row[16].split()[1]
                            if time != '00:00:00':
                                ts = time.split(':')
                                date += pretty_num(timedelta(hours=int(ts[0]), minutes=int(ts[1]), seconds=int(ts[2])).total_seconds()/(24*60*60), sig=6).lstrip('0')
                            add_quantity(events, name, 'discoverdate', date, source)
                    if args.update:
                        journal_events(tasks, args, events)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'rochester'):
            rochesterpaths = ['http://www.rochesterastronomy.org/snimages/snredshiftall.html', 'http://www.rochesterastronomy.org/sn2016/snredshift.html']
            rochesterupdate = [False, True]

            for p, path in enumerate(pbar(rochesterpaths, current_task)):
                if args.update and not rochesterupdate[p]:
                    continue

                filepath = os.path.join(PATH.REPO_EXTERNAL, 'rochester/') + os.path.basename(path)
                html = load_cached_url(args, path, filepath)
                if not html:
                    continue

                soup = BeautifulSoup(html, 'html5lib')
                rows = soup.findAll('tr')
                secondaryreference = 'Latest Supernovae'
                secondaryrefurl = 'http://www.rochesterastronomy.org/snimages/snredshiftall.html'
                for r, row in enumerate(pbar(rows, current_task)):
                    if r == 0:
                        continue
                    cols = row.findAll('td')
                    if not len(cols):
                        continue

                    name = ''
                    if cols[14].contents:
                        aka = str(cols[14].contents[0]).strip()
                        if is_number(aka.strip('?')):
                            aka = 'SN' + aka.strip('?') + 'A'
                            name = add_event(tasks, args, events, aka)
                        elif len(aka) >= 4 and is_number(aka[:4]):
                            aka = 'SN' + aka
                            name = add_event(tasks, args, events, aka)

                    ra = str(cols[3].contents[0]).strip()
                    dec = str(cols[4].contents[0]).strip()

                    sn = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
                    if is_number(sn.strip('?')):
                        sn = 'SN' + sn.strip('?') + 'A'
                    elif len(sn) >= 4 and is_number(sn[:4]):
                        sn = 'SN' + sn
                    if not name:
                        if not sn:
                            continue
                        if sn[:8] == 'MASTER J':
                            sn = sn.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                        if 'POSSIBLE' in sn.upper() and ra and dec:
                            sn = 'PSN J' + ra.replace(':', '').replace('.', '') + dec.replace(':', '').replace('.', '')
                        name = add_event(tasks, args, events, sn)

                    reference = cols[12].findAll('a')[0].contents[0].strip()
                    refurl = cols[12].findAll('a')[0]['href'].strip()
                    source = add_source(events, name, refname=reference, url=refurl)
                    secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, secondary=True)
                    sources = uniq_cdl(list(filter(None, [source, secondarysource])))
                    add_quantity(events, name, 'alias', name, sources)
                    add_quantity(events, name, 'alias', sn, sources)

                    if cols[14].contents:
                        if aka == 'SNR G1.9+0.3':
                            aka = 'G001.9+00.3'
                        if aka[:4] == 'PS1 ':
                            aka = 'PS1-' + aka[4:]
                        if aka[:8] == 'MASTER J':
                            aka = aka.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                        if 'POSSIBLE' in aka.upper() and ra and dec:
                            aka = 'PSN J' + ra.replace(':', '').replace('.', '') + dec.replace(':', '').replace('.', '')
                        add_quantity(events, name, 'alias', aka, sources)

                    if str(cols[1].contents[0]).strip() != 'unk':
                        add_quantity(events, name, 'claimedtype', str(cols[1].contents[0]).strip(' :,'), sources)
                    if str(cols[2].contents[0]).strip() != 'anonymous':
                        add_quantity(events, name, 'host', str(cols[2].contents[0]).strip(), sources)
                    add_quantity(events, name, 'ra', ra, sources)
                    add_quantity(events, name, 'dec', dec, sources)
                    if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
                        astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd').datetime
                        add_quantity(events, name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), sources)
                    if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
                        astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
                        if ((float(str(cols[8].contents[0]).strip()) <= 90.0 and
                             not any('GRB' in x for x in get_aliases(name)))):
                            add_photometry(events, name, time=str(astrot.mjd), magnitude=str(cols[8].contents[0]).strip(), source=sources)
                    if cols[11].contents[0] != 'n/a':
                        add_quantity(events, name, 'redshift', str(cols[11].contents[0]).strip(), sources)
                    add_quantity(events, name, 'discoverer', str(cols[13].contents[0]).strip(), sources)
                    if args.update:
                        journal_events(tasks, args, events)

            if not args.update:
                vsnetfiles = ['latestsne.dat']
                for vsnetfile in vsnetfiles:
                    f = open(os.path.join(PATH.REPO_EXTERNAL, "" + vsnetfile, 'r', encoding='latin1'))
                    tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
                    for r, row in enumerate(tsvin):
                        if not row or row[0][:4] in ['http', 'www.'] or len(row) < 3:
                            continue
                        name = row[0].strip()
                        if name[:4].isdigit():
                            name = 'SN' + name
                        if name.startswith('PSNJ'):
                            name = 'PSN J' + name[4:]
                        if name.startswith('MASTEROTJ'):
                            name = name.replace('MASTEROTJ', 'MASTER OT J')
                        name = add_event(tasks, args, events, name)
                        secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, secondary=True)
                        add_quantity(events, name, 'alias', name, secondarysource)

                        if not is_number(row[1]):
                            continue
                        year = row[1][:4]
                        month = row[1][4:6]
                        day = row[1][6:]
                        if '.' not in day:
                            day = day[:2] + '.' + day[2:]
                        mjd = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
                        magnitude = row[2].rstrip(ascii_letters)
                        if not is_number(magnitude):
                            continue
                        if magnitude.isdigit():
                            if int(magnitude) > 100:
                                magnitude = magnitude[:2] + '.' + magnitude[2:]

                        if float(str(cols[8].contents[0]).strip()) >= 90.0:
                            continue

                        if len(row) >= 4:
                            if is_number(row[3]):
                                e_magnitude = row[3]
                                refind = 4
                            else:
                                e_magnitude = ''
                                refind = 3

                            if refind >= len(row):
                                sources = secondarysource
                            else:
                                reference = ' '.join(row[refind:])
                                source = add_source(events, name, refname=reference)
                                add_quantity(events, name, 'alias', name, secondarysource)
                                sources = uniq_cdl([source, secondarysource])
                        else:
                            sources = secondarysource

                        band = row[2].lstrip('1234567890.')

                        add_photometry(events, name, time=mjd, band=band, magnitude=magnitude, e_magnitude=e_magnitude, source=sources)
                    f.close()
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ogle'):
            basenames = ['transients', 'transients/2014b', 'transients/2014', 'transients/2013', 'transients/2012']
            oglenames = []
            ogleupdate = [True, False, False, False, False]
            for b, bn in enumerate(pbar(basenames, current_task)):
                if args.update and not ogleupdate[b]:
                    continue

                filepath = os.path.join(PATH.REPO_EXTERNAL, 'OGLE-') + bn.replace('/', '-') + '-transients.html'
                htmltxt = load_cached_url(args, 'http://ogle.astrouw.edu.pl/ogle4/' + bn + '/transients.html', filepath)
                if not htmltxt:
                    continue

                soup = BeautifulSoup(htmltxt, 'html5lib')
                links = soup.findAll('a')
                breaks = soup.findAll('br')
                datalinks = []
                datafnames = []
                for a in links:
                    if a.has_attr('href'):
                        if '.dat' in a['href']:
                            datalinks.append('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/' + a['href'])
                            datafnames.append(bn.replace('/', '-') + '-' + a['href'].replace('/', '-'))

                ec = -1
                reference = 'OGLE-IV Transient Detection System'
                refurl = 'http://ogle.astrouw.edu.pl/ogle4/transients/transients.html'
                for br in pbar(breaks, current_task):
                    sibling = br.nextSibling
                    if 'Ra,Dec=' in sibling:
                        line = sibling.replace('\n', '').split('Ra,Dec=')
                        name = line[0].strip()
                        ec += 1

                        if 'NOVA' in name or 'dupl' in name:
                            continue

                        if name in oglenames:
                            continue
                        oglenames.append(name)

                        name = add_event(tasks, args, events, name)

                        mySibling = sibling.nextSibling
                        atelref = ''
                        claimedtype = ''
                        while 'Ra,Dec=' not in mySibling:
                            if isinstance(mySibling, NavigableString):
                                if 'Phot.class=' in str(mySibling):
                                    claimedtype = re.sub(r'\([^)]*\)', '', str(mySibling).split('=')[-1]).replace('SN', '').strip()
                            if isinstance(mySibling, Tag):
                                atela = mySibling
                                if atela and atela.has_attr('href') and 'astronomerstelegram' in atela['href']:
                                    atelref = atela.contents[0].strip()
                                    atelurl = atela['href']
                            mySibling = mySibling.nextSibling
                            if mySibling is None:
                                break

                        nextSibling = sibling.nextSibling
                        if isinstance(nextSibling, Tag) and nextSibling.has_attr('alt') and nextSibling.contents[0].strip() != 'NED':
                            radec = nextSibling.contents[0].strip().split()
                        else:
                            radec = line[-1].split()
                        ra = radec[0]
                        dec = radec[1]

                        fname = os.path.join(PATH.REPO_EXTERNAL, 'OGLE/') + datafnames[ec]
                        if not args.fullrefresh and archived_task(tasks, args, 'ogle') and os.path.isfile(fname):
                            with open(fname, 'r') as f:
                                csvtxt = f.read()
                        else:
                            response = urllib.request.urlopen(datalinks[ec])
                            with open(fname, 'w') as f:
                                csvtxt = response.read().decode('utf-8')
                                f.write(csvtxt)

                        lcdat = csvtxt.splitlines()
                        sources = [add_source(events, name, refname=reference, url=refurl)]
                        add_quantity(events, name, 'alias', name, sources[0])
                        if atelref and atelref != 'ATel#----':
                            sources.append(add_source(events, name, refname=atelref, url=atelurl))
                        sources = uniq_cdl(sources)

                        if name.startswith('OGLE'):
                            if name[4] == '-':
                                if is_number(name[5:9]):
                                    add_quantity(events, name, 'discoverdate', name[5:9], sources)
                            else:
                                if is_number(name[4:6]):
                                    add_quantity(events, name, 'discoverdate', '20' + name[4:6], sources)

                        # RA and Dec from OGLE pages currently not reliable
                        # add_quantity(events, name, 'ra', ra, sources)
                        # add_quantity(events, name, 'dec', dec, sources)
                        if claimedtype and claimedtype != '-':
                            add_quantity(events, name, 'claimedtype', claimedtype, sources)
                        elif 'SN' not in name and 'claimedtype' not in events[name]:
                            add_quantity(events, name, 'claimedtype', 'Candidate', sources)
                        for row in lcdat:
                            row = row.split()
                            mjd = str(jd_to_mjd(Decimal(row[0])))
                            magnitude = row[1]
                            if float(magnitude) > 90.0:
                                continue
                            e_magnitude = row[2]
                            upperlimit = False
                            if e_magnitude == '-1' or float(e_magnitude) > 10.0:
                                e_magnitude = ''
                                upperlimit = True
                            add_photometry(
                                events, name, time=mjd, band='I', magnitude=magnitude, e_magnitude=e_magnitude,
                                system='Vega', source=sources, upperlimit=upperlimit)
                        if args.update:
                            journal_events(tasks, args, events)
                journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snls'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'SNLS-ugriz.dat'), 'r') as f:
                data = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace=True)
                for row in data:
                    flux = row[3]
                    err = row[4]
                    # Being extra strict here with the flux constraint, see note below.
                    if float(flux) < 3.0*float(err):
                        continue
                    name = 'SNLS-' + row[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2010A&A...523A...7G')
                    add_quantity(events, name, 'alias', name, source)
                    band = row[1]
                    mjd = row[2]
                    sig = get_sig_digits(flux.split('E')[0])+1
                    # Conversion comes from SNLS-Readme
                    # NOTE: Datafiles available for download suggest different zeropoints than 30, need to inquire.
                    magnitude = pretty_num(30.0-2.5*log10(float(flux)), sig=sig)
                    e_magnitude = pretty_num(2.5*log10(1.0 + float(err)/float(flux)), sig=sig)
                    # e_magnitude=pretty_num(2.5*(log10(float(flux) + float(err)) - log10(float(flux))), sig=sig)
                    add_photometry(
                        events, name, time=mjd, band=band, magnitude=magnitude, e_magnitude=e_magnitude, counts=flux,
                        e_counts=err, source=source)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'psthreepi'):
            fname = os.path.join(PATH.REPO_EXTERNAL, '3pi/page00.html')
            html = load_cached_url(args, 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=1&sort=followup_flag_date', fname, write=False)
            if not html:
                continue

            bs = BeautifulSoup(html, 'html5lib')
            div = bs.find('div', {'class': 'pagination'})
            offline = False
            if not div:
                offline = True
            else:
                links = div.findAll('a')
                if not links:
                    offline = True

            if offline:
                if args.update:
                    continue
                warnings.warn('Pan-STARRS 3pi offline, using local files only.')
                with open(fname, 'r') as f:
                    html = f.read()
                bs = BeautifulSoup(html, 'html5lib')
                div = bs.find('div', {'class': 'pagination'})
                links = div.findAll('a')
            else:
                with open(fname, 'w') as f:
                    f.write(html)

            numpages = int(links[-2].contents[0])
            oldnumpages = len(glob(os.path.join(PATH.REPO_EXTERNAL, '3pi/page*')))
            for page in pbar(range(1, numpages), current_task):
                fname = os.path.join(PATH.REPO_EXTERNAL, '3pi/page') + str(page).zfill(2) + '.html'
                if not args.fullrefresh and archived_task(tasks, args, 'psthreepi') and os.path.isfile(fname) and page < oldnumpages:
                    with open(fname, 'r') as f:
                        html = f.read()
                elif not offline:
                    response = urllib.request.urlopen('http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=' + str(page) + '&sort=followup_flag_date')
                    with open(fname, 'w') as f:
                        html = response.read().decode('utf-8')
                        f.write(html)
                else:
                    continue

                bs = BeautifulSoup(html, 'html5lib')
                trs = bs.findAll('tr')
                for tr in pbar(trs, current_task):
                    tds = tr.findAll('td')
                    if not tds:
                        continue
                    refs = []
                    aliases = []
                    ttype = ''
                    ctype = ''
                    for tdi, td in enumerate(tds):
                        if tdi == 0:
                            psname = td.contents[0]
                            pslink = psname['href']
                            psname = psname.text
                        elif tdi == 1:
                            ra = td.contents[0]
                        elif tdi == 2:
                            dec = td.contents[0]
                        elif tdi == 3:
                            ttype = td.contents[0]
                            if ttype != 'sn' and ttype != 'orphan':
                                break
                        elif tdi == 5:
                            if not td.contents:
                                continue
                            ctype = td.contents[0]
                            if ctype == 'Observed':
                                ctype = ''
                        elif tdi == 16:
                            if td.contents:
                                crossrefs = td.findAll('a')
                                for cref in crossrefs:
                                    if 'atel' in cref.contents[0].lower():
                                        refs.append([cref.contents[0], cref['href']])
                                    elif is_number(cref.contents[0][:4]):
                                        continue
                                    else:
                                        aliases.append(cref.contents[0])

                    if ttype != 'sn' and ttype != 'orphan':
                        continue

                    name = ''
                    for alias in aliases:
                        if alias[:2] == 'SN':
                            name = alias
                    if not name:
                        name = psname
                    name = add_event(tasks, args, events, name)
                    sources = [add_source(events, name, refname='Pan-STARRS 3Pi', url='http://psweb.mp.qub.ac.uk/ps1threepi/psdb/')]
                    add_quantity(events, name, 'alias', name, sources[0])
                    for ref in refs:
                        sources.append(add_source(events, name, refname=ref[0], url=ref[1]))
                    source = uniq_cdl(sources)
                    for alias in aliases:
                        newalias = alias
                        if alias[:3] in ['CSS', 'SSS', 'MLS']:
                            newalias = alias.replace('-', ':', 1)
                        newalias = newalias.replace('PSNJ', 'PSN J')
                        add_quantity(events, name, 'alias', newalias, source)
                    add_quantity(events, name, 'ra', ra, source)
                    add_quantity(events, name, 'dec', dec, source)
                    add_quantity(events, name, 'claimedtype', ctype, source)

                    fname2 = os.path.join(PATH.REPO_EXTERNAL, '3pi/candidate-') + pslink.rstrip('/').split('/')[-1] + '.html'
                    if archived_task(tasks, args, 'psthreepi') and os.path.isfile(fname2):
                        with open(fname2, 'r') as f:
                            html2 = f.read()
                    elif not offline:
                        pslink = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/' + pslink
                        with open(fname2, 'w') as f:
                            response2 = urllib.request.urlopen(pslink)
                            html2 = response2.read().decode('utf-8')
                            f.write(html2)
                    else:
                        continue

                    bs2 = BeautifulSoup(html2, 'html5lib')
                    scripts = bs2.findAll('script')
                    nslines = []
                    nslabels = []
                    for script in scripts:
                        if 'jslcdata.push' not in script.text:
                            continue
                        slines = script.text.splitlines()
                        for line in slines:
                            if 'jslcdata.push' in line:
                                nslines.append(json.loads(line.strip().replace('jslcdata.push(', '').replace(');', '')))
                            if 'jslabels.push' in line and 'blanks' not in line and 'non det' not in line:
                                nslabels.append(json.loads(line.strip().replace('jslabels.push(', '').replace(');', ''))['label'])
                    for li, line in enumerate(nslines[:len(nslabels)]):
                        if not line:
                            continue
                        for obs in line:
                            add_photometry(
                                events, name, time=str(obs[0]), band=nslabels[li], magnitude=str(obs[1]), e_magnitude=str(obs[2]), source=source,
                                telescope='Pan-STARRS1')
                    for li, line in enumerate(nslines[2*len(nslabels):]):
                        if not line:
                            continue
                        for obs in line:
                            add_photometry(
                                events, name, time=str(obs[0]), band=nslabels[li], magnitude=str(obs[1]), upperlimit=True, source=source,
                                telescope='Pan-STARRS1')
                    assoctab = bs2.find('table', {'class': 'generictable'})
                    hostname = ''
                    redshift = ''
                    if assoctab:
                        trs = assoctab.findAll('tr')
                        headertds = [x.contents[0] for x in trs[1].findAll('td')]
                        tds = trs[1].findAll('td')
                        for tdi, td in enumerate(tds):
                            if tdi == 1:
                                hostname = td.contents[0].strip()
                            elif tdi == 4:
                                if 'z' in headertds:
                                    redshift = td.contents[0].strip()
                    # Skip galaxies with just SDSS id
                    if is_number(hostname):
                        continue
                    add_quantity(events, name, 'host', hostname, source)
                    if redshift:
                        add_quantity(events, name, 'redshift', redshift, source, kind='host')
                    if args.update:
                        journal_events(tasks, args, events)
                journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'psmds'):
            with open(os.path.join(PATH.REPO_EXTERNAL, 'MDS/apj506838t1_mrt.txt')) as f:
                for ri, row in enumerate(pbar(f.read().splitlines(), current_task)):
                    if ri < 35:
                        continue
                    cols = [x.strip() for x in row.split(',')]
                    name = add_event(tasks, args, events, cols[0])
                    source = add_source(events, name, bibcode='2015ApJ...799..208S')
                    add_quantity(events, name, 'alias', name, source)
                    add_quantity(events, name, 'ra', cols[2], source)
                    add_quantity(events, name, 'dec', cols[3], source)
                    astrot = astrotime(float(cols[4]), format='mjd').datetime
                    add_quantity(events, name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
                    add_quantity(events, name, 'redshift', cols[5], source, kind='spectroscopic')
                    add_quantity(events, name, 'claimedtype', 'II P', source)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'crts'):
            crtsnameerrors = ['2011ax']

            folders = ['catalina', 'MLS', 'SSS']
            for fold in pbar(folders, current_task):
                html = load_cached_url(args, 'http://nesssi.cacr.caltech.edu/' + fold + '/AllSN.html', os.path.join(PATH.REPO_EXTERNAL, 'CRTS/') + fold + '.html')
                if not html:
                    continue
                bs = BeautifulSoup(html, 'html5lib')
                trs = bs.findAll('tr')
                for tr in pbar(trs, current_task):
                    tds = tr.findAll('td')
                    if not tds:
                        continue
                    refs = []
                    aliases = []
                    ttype = ''
                    ctype = ''
                    for tdi, td in enumerate(tds):
                        if tdi == 0:
                            crtsname = td.contents[0].text.strip()
                        elif tdi == 1:
                            ra = td.contents[0]
                        elif tdi == 2:
                            dec = td.contents[0]
                        elif tdi == 11:
                            lclink = td.find('a')['onclick']
                            lclink = lclink.split("'")[1]
                        elif tdi == 13:
                            aliases = re.sub('[()]', '', re.sub('<[^<]+?>', '', td.contents[-1].strip()))
                            aliases = [x.strip('; ') for x in list(filter(None, aliases.split(' ')))]

                    name = ''
                    hostmag = ''
                    hostupper = False
                    validaliases = []
                    for ai, alias in enumerate(aliases):
                        if alias in ['SN', 'SDSS']:
                            continue
                        if alias in crtsnameerrors:
                            continue
                        if alias == 'mag':
                            if ai < len(aliases) - 1:
                                ind = ai+1
                                if aliases[ai+1] in ['SDSS']:
                                    ind = ai+2
                                elif aliases[ai+1] in ['gal', 'obj', 'object', 'source']:
                                    ind = ai-1
                                if '>' in aliases[ind]:
                                    hostupper = True
                                hostmag = aliases[ind].strip('>~').replace(',', '.')
                            continue
                        if is_number(alias[:4]) and alias[:2] == '20' and len(alias) > 4:
                            name = 'SN' + alias
                        # lalias = alias.lower()
                        if ((('asassn' in alias and len(alias) > 6) or ('ptf' in alias and len(alias) > 3) or
                             ('ps1' in alias and len(alias) > 3) or 'snhunt' in alias or
                             ('mls' in alias and len(alias) > 3) or 'gaia' in alias or ('lsq' in alias and len(alias) > 3))):
                            alias = alias.replace('SNHunt', 'SNhunt')
                            validaliases.append(alias)
                    if not name:
                        name = crtsname
                    name = add_event(tasks, args, events, name)
                    source = add_source(
                        events, name, refname='Catalina Sky Survey', bibcode='2009ApJ...696..870D',
                        url='http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
                    add_quantity(events, name, 'alias', name, source)
                    for alias in validaliases:
                        add_quantity(events, name, 'alias', alias, source)
                    add_quantity(events, name, 'ra', ra, source, unit='floatdegrees')
                    add_quantity(events, name, 'dec', dec, source, unit='floatdegrees')

                    if hostmag:
                        # 1.0 magnitude error based on Drake 2009 assertion that SN are only considered real if they are 2 mags brighter than host.
                        add_photometry(
                            events, name, band='C', magnitude=hostmag, e_magnitude=1.0, source=source, host=True,
                            telescope='Catalina Schmidt', upperlimit=hostupper)

                    fname2 = PATH.REPO_EXTERNAL + '/' + fold + '/' + lclink.split('.')[-2].rstrip('p').split('/')[-1] + '.html'
                    if not args.fullrefresh and archived_task(tasks, args, 'crts') and os.path.isfile(fname2):
                        with open(fname2, 'r') as f:
                            html2 = f.read()
                    else:
                        with open(fname2, 'w') as f:
                            response2 = urllib.request.urlopen(lclink)
                            html2 = response2.read().decode('utf-8')
                            f.write(html2)

                    lines = html2.splitlines()
                    for line in lines:
                        if 'javascript:showx' in line:
                            mjdstr = re.search("showx\('(.*?)'\)", line).group(1).split('(')[0].strip()
                            if not is_number(mjdstr):
                                continue
                            mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                        else:
                            continue
                        if 'javascript:showy' in line:
                            mag = re.search("showy\('(.*?)'\)", line).group(1)
                        if 'javascript:showz' in line:
                            err = re.search("showz\('(.*?)'\)", line).group(1)
                        add_photometry(
                            events, name, time=mjd, band='C', magnitude=mag, source=source, includeshost=True,
                            telescope='Catalina Schmidt', e_magnitude=err if float(err) > 0.0 else '', upperlimit=(float(err) == 0.0))
                    if args.update:
                        journal_events(tasks, args, events)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snhunt'):
            html = load_cached_url(args, 'http://nesssi.cacr.caltech.edu/catalina/current.html', os.path.join(PATH.REPO_EXTERNAL, 'SNhunt/current.html'))
            if not html:
                continue
            text = html.splitlines()
            findtable = False
            for ri, row in enumerate(text):
                if 'Supernova Discoveries' in row:
                    findtable = True
                if findtable and '<table' in row:
                    tstart = ri+1
                if findtable and '</table>' in row:
                    tend = ri-1
            tablestr = '<html><body><table>'
            for row in text[tstart:tend]:
                if row[:3] == 'tr>':
                    tablestr = tablestr + '<tr>' + row[3:]
                else:
                    tablestr = tablestr + row
            tablestr = tablestr + '</table></body></html>'
            bs = BeautifulSoup(tablestr, 'html5lib')
            trs = bs.find('table').findAll('tr')
            for tr in pbar(trs, current_task):
                cols = [str(x.text) for x in tr.findAll('td')]
                if not cols:
                    continue
                name = re.sub('<[^<]+?>', '', cols[4]).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, refname='Supernova Hunt', url='http://nesssi.cacr.caltech.edu/catalina/current.html')
                add_quantity(events, name, 'alias', name, source)
                host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
                add_quantity(events, name, 'host', host, source)
                add_quantity(events, name, 'ra', cols[2], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', cols[3], source, unit='floatdegrees')
                dd = cols[0]
                discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
                add_quantity(events, name, 'discoverdate', discoverdate, source)
                discoverers = cols[5].split('/')
                for discoverer in discoverers:
                    add_quantity(events, name, 'discoverer', 'CRTS', source)
                    add_quantity(events, name, 'discoverer', discoverer, source)
                if args.update:
                    journal_events(tasks, args, events)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'nedd'):
            f = open(os.path.join(PATH.REPO_EXTERNAL, 'NED25.12.1-D-10.4.0-20151123.csv'), 'r')
            data = csv.reader(f, delimiter=',', quotechar='"')
            reference = 'NED-D'
            refurl = 'http://ned.ipac.caltech.edu/Library/Distances/'
            nedd_dict = OrderedDict()
            oldhostname = ''
            for r, row in enumerate(data):
                if r <= 12:
                    continue
                hostname = row[3]
                if args.update and oldhostname != hostname:
                    journal_events(tasks, args, events)
                # distmod = row[4]
                # moderr = row[5]
                dist = row[6]
                bibcode = unescape(row[8])
                name = ''
                if hostname.startswith('SN '):
                    if is_number(hostname[3:7]):
                        name = 'SN' + hostname[3:]
                    else:
                        name = hostname[3:]
                elif hostname.startswith('SNLS '):
                    name = 'SNLS-' + hostname[5:].split()[0]
                else:
                    cleanhost = hostname.replace('MESSIER 0', 'M').replace('MESSIER ', 'M').strip()
                    if True in [x in cleanhost for x in ['UGC', 'PGC', 'IC']]:
                        cleanhost = ' '.join([x.lstrip('0') for x in cleanhost.split()])
                    if 'ESO' in cleanhost:
                        cleanhost = cleanhost.replace(' ', '').replace('ESO', 'ESO ')
                    nedd_dict.setdefault(cleanhost, []).append(Decimal(dist))

                if name:
                    name = add_event(tasks, args, events, name)
                    secondarysource = add_source(events, name, refname=reference, url=refurl, secondary=True)
                    add_quantity(events, name, 'alias', name, secondarysource)
                    if bibcode:
                        source = add_source(events, name, bibcode=bibcode)
                        sources = uniq_cdl([source, secondarysource])
                    else:
                        sources = secondarysource
                    add_quantity(events, name, 'comovingdist', dist, sources)
                oldhostname = hostname
            journal_events(tasks, args, events)

        # Import CPCS
        if do_task(tasks, args, task, 'cpcs'):
            jsontxt = load_cached_url(
                args, 'http://gsaweb.ast.cam.ac.uk/followup/list_of_alerts?format=json&num=100000&published=1&observed_only=1&hashtag=JG_530ad9462a0b8785bfb385614bf178c6',
                os.path.join(PATH.REPO_EXTERNAL, 'CPCS/index.json'))
            if not jsontxt:
                continue
            alertindex = json.loads(jsontxt, object_pairs_hook=OrderedDict)
            ids = [x['id'] for x in alertindex]

            for i, ai in enumerate(pbar(ids, current_task)):
                name = alertindex[i]['ivorn'].split('/')[-1].strip()
                # Skip a few weird entries
                if name == 'ASASSNli':
                    continue
                # Just use a whitelist for now since naming seems inconsistent
                if True in [x in name.upper() for x in ['GAIA', 'OGLE', 'ASASSN', 'MASTER', 'OTJ', 'PS1', 'IPTF']]:
                    name = name.replace('Verif', '').replace('_', ' ')
                    if 'ASASSN' in name and name[6] != '-':
                        name = 'ASASSN-' + name[6:]
                    if 'MASTEROTJ' in name:
                        name = name.replace('MASTEROTJ', 'MASTER OT J')
                    if 'OTJ' in name:
                        name = name.replace('OTJ', 'MASTER OT J')
                    if name.upper().startswith('IPTF'):
                        name = 'iPTF' + name[4:]
                    # Only add events that are classified as SN.
                    if event_exists(events, name):
                        continue
                    name = add_event(tasks, args, events, name)
                else:
                    continue

                secondarysource = add_source(events, name, refname='Cambridge Photometric Calibration Server', url='http://gsaweb.ast.cam.ac.uk/followup/', secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                add_quantity(events, name, 'ra', str(alertindex[i]['ra']), secondarysource, unit='floatdegrees')
                add_quantity(events, name, 'dec', str(alertindex[i]['dec']), secondarysource, unit='floatdegrees')

                alerturl = 'http://gsaweb.ast.cam.ac.uk/followup/get_alert_lc_data?alert_id=' + str(ai)
                source = add_source(events, name, refname='CPCS Alert ' + str(ai), url=alerturl)
                fname = os.path.join(PATH.REPO_EXTERNAL, 'CPCS/alert-') + str(ai).zfill(2) + '.json'
                if archived_task(tasks, args, 'cpcs') and os.path.isfile(fname):
                    with open(fname, 'r') as f:
                        jsonstr = f.read()
                else:
                    session = requests.Session()
                    response = session.get(alerturl + '&hashtag=JG_530ad9462a0b8785bfb385614bf178c6')
                    with open(fname, 'w') as f:
                        jsonstr = response.text
                        f.write(jsonstr)

                try:
                    cpcsalert = json.loads(jsonstr)
                except:
                    continue

                mjds = [round_sig(x, sig=9) for x in cpcsalert['mjd']]
                mags = [round_sig(x, sig=6) for x in cpcsalert['mag']]
                errs = [round_sig(x, sig=6) if (is_number(x) and float(x) > 0.0) else '' for x in cpcsalert['magerr']]
                bnds = cpcsalert['filter']
                obs  = cpcsalert['observatory']
                for mi, mjd in enumerate(mjds):
                    add_photometry(
                        events, name, time=mjd, magnitude=mags[mi], e_magnitude=errs[mi],
                        band=bnds[mi], observatory=obs[mi], source=uniq_cdl([source, secondarysource]))
                if args.update:
                    journal_events(tasks, args, events)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ptf'):
            # response = urllib.request.urlopen('http://wiserep.weizmann.ac.il/objects/list')
            # bs = BeautifulSoup(response, 'html5lib')
            # select = bs.find('select', {'name': 'objid'})
            # options = select.findAll('option')
            # for option in options:
            #    print(option.text)
            #    name = option.text
            #    if ((name.startswith('PTF') and is_number(name[3:5])) or
            #        name.startswith('PTFS') or name.startswith('iPTF')):
            #        name = add_event(tasks, args, events, name)

            if archived_task(tasks, args, 'ptf'):
                with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'r') as f:
                    html = f.read()
            else:
                session = requests.Session()
                response = session.get('http://wiserep.weizmann.ac.il/spectra/update')
                html = response.text
                with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'w') as f:
                    f.write(html)

            bs = BeautifulSoup(html, 'html5lib')
            select = bs.find('select', {'name': 'objid'})
            options = select.findAll('option')
            for option in options:
                name = option.text
                if (((name.startswith('PTF') and is_number(name[3:5])) or
                     name.startswith('PTFS') or name.startswith('iPTF'))):
                    if '(' in name:
                        alias = name.split('(')[0].strip(' ')
                        name = name.split('(')[-1].strip(') ').replace('sn', 'SN')
                        name = add_event(tasks, args, events, name)
                        source = add_source(events, name, bibcode='2012PASP..124..668Y')
                        add_quantity(events, name, 'alias', alias, source)
                    else:
                        name = add_event(tasks, args, events, name)

            with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/old-ptf-events.csv')) as f:
                for suffix in f.read().splitlines():
                    name = add_event(tasks, args, events, 'PTF' + suffix)
            with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/perly-2016.csv')) as f:
                for row in f.read().splitlines():
                    cols = [x.strip() for x in row.split(',')]
                    alias = ''
                    if cols[8]:
                        name = cols[8]
                        alias = 'PTF' + cols[0]
                    else:
                        name = 'PTF' + cols[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2016arXiv160408207P')
                    add_quantity(events, name, 'alias', name, source)
                    if alias:
                        add_quantity(events, name, 'alias', alias, source)
                    add_quantity(events, name, 'ra', cols[1], source)
                    add_quantity(events, name, 'dec', cols[2], source)
                    add_quantity(events, name, 'claimedtype', 'SLSN-' + cols[3], source)
                    add_quantity(events, name, 'redshift', cols[4], source, kind='spectroscopic')
                    maxdate = cols[6].replace('-', '/')
                    add_quantity(events, name, 'maxdate', maxdate.lstrip('<'), source, upperlimit=maxdate.startswith('<'))
                    add_quantity(events, name, 'ebv', cols[7], source, kind='spectroscopic')
                    name = add_event(tasks, args, events, 'PTF' + suffix)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'des'):
            html = load_cached_url(args, 'https://portal.nersc.gov/des-sn/transients/', os.path.join(PATH.REPO_EXTERNAL, 'DES/transients.html'))
            if not html:
                continue
            bs = BeautifulSoup(html, 'html5lib')
            trs = bs.find('tbody').findAll('tr')
            for tri, tr in enumerate(pbar(trs, current_task)):
                name = ''
                source = ''
                if tri == 0:
                    continue
                tds = tr.findAll('td')
                for tdi, td in enumerate(tds):
                    if tdi == 0:
                        name = add_event(tasks, args, events, td.text.strip())
                    if tdi == 1:
                        (ra, dec) = [x.strip() for x in td.text.split('\xa0')]
                    if tdi == 6:
                        atellink = td.find('a')
                        if atellink:
                            atellink = atellink['href']
                        else:
                            atellink = ''

                sources = [add_source(events, name, url='https://portal.nersc.gov/des-sn/', refname='DES Bright Transients',
                                      acknowledgment='http://www.noao.edu/noao/library/NOAO_Publications_Acknowledgments.html#DESdatause')]
                if atellink:
                    sources.append(add_source(events, name, refname='ATel ' + atellink.split('=')[-1], url=atellink))
                sources += [add_source(events, name, bibcode='2012ApJ...753..152B'),
                            add_source(events, name, bibcode='2015AJ....150..150F'),
                            add_source(events, name, bibcode='2015AJ....150...82G'),
                            add_source(events, name, bibcode='2015AJ....150..172K')]
                sources = ','.join(sources)
                add_quantity(events, name, 'alias', name, sources)
                add_quantity(events, name, 'ra', ra, sources)
                add_quantity(events, name, 'dec', dec, sources)

                html2 = load_cached_url(args, 'https://portal.nersc.gov/des-sn/transients/' + name, os.path.join(PATH.REPO_EXTERNAL, 'DES/') + name + '.html')
                if not html2:
                    continue
                lines = html2.splitlines()
                for line in lines:
                    if 'var data = ' in line:
                        jsontxt = json.loads(line.split('=')[-1].rstrip(';'))
                        for i, band in enumerate(jsontxt['band']):
                            add_photometry(
                                events, name, time=jsontxt['mjd'][i], magnitude=jsontxt['mag'][i], e_magnitude=jsontxt['mag_error'][i],
                                band=band, observatory='CTIO', telescope='Blanco 4m', instrument='DECam',
                                upperlimit=True if float(jsontxt['snr'][i]) <= 3.0 else '', source=sources)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'asassn'):
            html = load_cached_url(args, 'http://www.astronomy.ohio-state.edu/~assassin/sn_list.html', os.path.join(PATH.REPO_EXTERNAL, 'ASASSN/sn_list.html'))
            if not html:
                continue
            bs = BeautifulSoup(html, 'html5lib')
            trs = bs.find('table').findAll('tr')
            for tri, tr in enumerate(pbar(trs, current_task)):
                name = ''
                source = ''
                ra = ''
                dec = ''
                redshift = ''
                hostoff = ''
                claimedtype = ''
                host = ''
                atellink = ''
                typelink = ''
                if tri == 0:
                    continue
                tds = tr.findAll('td')
                for tdi, td in enumerate(tds):
                    if tdi == 1:
                        name = add_event(tasks, args, events, td.text.strip())
                        atellink = td.find('a')
                        if atellink:
                            atellink = atellink['href']
                        else:
                            atellink = ''
                    if tdi == 2:
                        discdate = td.text.replace('-', '/')
                    if tdi == 3:
                        ra = td.text
                    if tdi == 4:
                        dec = td.text
                    if tdi == 5:
                        redshift = td.text
                    if tdi == 8:
                        hostoff = td.text
                    if tdi == 9:
                        claimedtype = td.text
                        typelink = td.find('a')
                        if typelink:
                            typelink = typelink['href']
                        else:
                            typelink = ''
                    if tdi == 12:
                        host = td.text

                sources = [add_source(events, name, url='http://www.astronomy.ohio-state.edu/~assassin/sn_list.html', refname='ASAS-SN Supernovae')]
                typesources = sources[:]
                if atellink:
                    sources.append(add_source(events, name, refname='ATel ' + atellink.split('=')[-1], url=atellink))
                if typelink:
                    typesources.append(add_source(events, name, refname='ATel ' + typelink.split('=')[-1], url=typelink))
                sources = ','.join(sources)
                typesources = ','.join(typesources)
                add_quantity(events, name, 'alias', name, sources)
                add_quantity(events, name, 'discoverdate', discdate, sources)
                add_quantity(events, name, 'ra', ra, sources, unit='floatdegrees')
                add_quantity(events, name, 'dec', dec, sources, unit='floatdegrees')
                add_quantity(events, name, 'redshift', redshift, sources)
                add_quantity(events, name, 'hostoffset', hostoff, sources, unit='arcseconds')
                for ct in claimedtype.split('/'):
                    if ct != 'Unk':
                        add_quantity(events, name, 'claimedtype', ct, typesources)
                if host != 'Uncatalogued':
                    add_quantity(events, name, 'host', host, sources)
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'asiagospectra'):
            html = load_cached_url(args, 'http://sngroup.oapd.inaf.it./cgi-bin/output_class.cgi?sn=1990',
                                   os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Asiago/spectra.html'))
            if not html:
                continue
            bs = BeautifulSoup(html, 'html5lib')
            trs = bs.findAll('tr')
            for tr in pbar(trs, current_task):
                tds = tr.findAll('td')
                name = ''
                host = ''
                # fitsurl = ''
                source = ''
                reference = ''
                for tdi, td in enumerate(tds):
                    if tdi == 0:
                        butt = td.find('button')
                        if not butt:
                            break
                        alias = butt.text.strip()
                        alias = alias.replace('PSNJ', 'PSN J').replace('GAIA', 'Gaia')
                    elif tdi == 1:
                        name = td.text.strip().replace('PSNJ', 'PSN J').replace('GAIA', 'Gaia')
                        if name.startswith('SN '):
                            name = 'SN' + name[3:]
                        if not name:
                            name = alias
                        if is_number(name[:4]):
                            name = 'SN' + name
                        name = add_event(tasks, args, events, name)
                        reference = 'Asiago Supernova Catalogue'
                        refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                        secondarysource = add_source(events, name, refname=reference, url=refurl, secondary=True)
                        add_quantity(events, name, 'alias', name, secondarysource)
                        if alias != name:
                            add_quantity(events, name, 'alias', alias, secondarysource)
                    elif tdi == 2:
                        host = td.text.strip()
                        if host == 'anonymous':
                            host = ''
                    elif tdi == 3:
                        discoverer = td.text.strip()
                    elif tdi == 5:
                        ra = td.text.strip()
                    elif tdi == 6:
                        dec = td.text.strip()
                    elif tdi == 7:
                        claimedtype = td.text.strip()
                    elif tdi == 8:
                        redshift = td.text.strip()
                    elif tdi == 9:
                        epochstr = td.text.strip()
                        if epochstr:
                            mjd = (astrotime(epochstr[:4] + '-' + epochstr[4:6] + '-' + str(floor(float(epochstr[6:]))).zfill(2)).mjd +
                                   float(epochstr[6:]) - floor(float(epochstr[6:])))
                        else:
                            mjd = ''
                    elif tdi == 10:
                        refs = td.findAll('a')
                        source = ''
                        reference = ''
                        refurl = ''
                        for ref in refs:
                            if ref.text != 'REF':
                                reference = ref.text
                                refurl = ref['href']
                        if reference:
                            source = add_source(events, name, refname=reference, url=refurl)
                        add_quantity(events, name, 'alias', name, secondarysource)
                        sources = uniq_cdl(list(filter(None, [source, secondarysource])))
                    elif tdi == 12:
                        pass
                        # fitslink = td.find('a')
                        # if fitslink:
                        #     fitsurl = fitslink['href']
                if name:
                    add_quantity(events, name, 'claimedtype', claimedtype, sources)
                    add_quantity(events, name, 'ra', ra, sources)
                    add_quantity(events, name, 'dec', dec, sources)
                    add_quantity(events, name, 'redshift', redshift, sources)
                    add_quantity(events, name, 'discoverer', discoverer, sources)
                    add_quantity(events, name, 'host', host, sources)

                    # if fitsurl:
                    #    response = urllib.request.urlopen('http://sngroup.oapd.inaf.it./' + fitsurl)
                    #    compressed = io.BytesIO(response.read())
                    #    decompressed = gzip.GzipFile(fileobj=compressed)
                    #    hdulist = fits.open(decompressed)
                    #    scidata = hdulist[0].data
                    #    print(hdulist[0].header)
                    #
                    #    print(scidata[3])
                    #    sys.exit()
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'wiserepspectra'):
            secondaryreference = 'WISeREP'
            secondaryrefurl = 'http://wiserep.weizmann.ac.il/'
            secondarybibcode = '2012PASP..124..668Y'
            wiserepcnt = 0

            # These are known to be in error on the WISeREP page, either fix or ignore them.
            wiserepbibcorrectdict = {'2000AJ....120..367G]': '2000AJ....120..367G',
                                     'Harutyunyan+et+al.+2008': '2008A&A...488..383H',
                                     '0609268': '2007AJ....133...58K',
                                     '2006ApJ...636...400Q': '2006ApJ...636..400Q',
                                     '2011ApJ...741...76': '2011ApJ...741...76C',
                                     '2016PASP...128...961': '2016PASP..128...961',
                                     '2002AJ....1124..417H': '2002AJ....1124.417H',
                                     '2013ApJ…774…58D': '2013ApJ...774...58D',
                                     '2011Sci.333..856S': '2011Sci...333..856S',
                                     '2014MNRAS.438,368': '2014MNRAS.438..368T',
                                     '2012MNRAS.420.1135': '2012MNRAS.420.1135S',
                                     '2012Sci..337..942D': '2012Sci...337..942D',
                                     'stt1839': ''}

            oldname = ''
            file_names = next(os.walk(PATH.REPO_EXTERNAL_WISEREP))[1]
            for folder in pbar_strings(file_names, current_task):
                files = glob(PATH.REPO_EXTERNAL_WISEREP + '/' + folder + '/*')
                for fname in pbar(files, current_task):
                    if '.html' in fname:
                        lfiles = deepcopy(files)
                        with open(fname, 'r') as f:
                            path = os.path.abspath(fname)
                            response = urllib.request.urlopen('file://' + path)
                            bs = BeautifulSoup(response, 'html5lib')
                            trs = bs.findAll('tr', {'valign': 'top'})
                            for tri, tr in enumerate(trs):
                                if 'Click to show/update object' in str(tr.contents):
                                    claimedtype = ''
                                    instrument = ''
                                    epoch = ''
                                    observer = ''
                                    reducer = ''
                                    specfile = ''
                                    produceoutput = True
                                    specpath = ''
                                    tds = tr.findAll('td')
                                    for tdi, td in enumerate(tds):
                                        if td.contents:
                                            if tdi == 3:
                                                name = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            elif tdi == 5:
                                                claimedtype = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                                if claimedtype == 'SN':
                                                    claimedtype = ''
                                                    continue
                                                if claimedtype[:3] == 'SN ':
                                                    claimedtype = claimedtype[3:].strip()
                                                claimedtype = claimedtype.replace('-like', '').strip()
                                            elif tdi == 9:
                                                instrument = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            elif tdi == 11:
                                                epoch = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            elif tdi == 13:
                                                observer = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                                if observer == 'Unknown' or observer == 'Other':
                                                    observer = ''
                                            elif tdi == 17:
                                                reducer = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                                if reducer == 'Unknown' or reducer == 'Other':
                                                    reducer = ''
                                            elif tdi == 25:
                                                speclinks = td.findAll('a')
                                                try:
                                                    for link in speclinks:
                                                        if 'Ascii' in link['href']:
                                                            specfile = link.contents[0].strip()
                                                            tfiles = deepcopy(lfiles)
                                                            for fi, fname in enumerate(lfiles):
                                                                if specfile in fname:
                                                                    specpath = fname
                                                                    del(tfiles[fi])
                                                                    lfiles = deepcopy(tfiles)
                                                                    raise(StopIteration)
                                                except StopIteration:
                                                    pass
                                                if not specpath:
                                                    warnings.warn('Spectrum file not found, '' + specfile + ''')
                                            else:
                                                continue
                                if 'Spec Type:</span>' in str(tr.contents) and produceoutput:
                                    produceoutput = False

                                    trstr = str(tr)
                                    result = re.search('redshift=(.*?)&amp;', trstr)
                                    redshift = ''
                                    if result:
                                        redshift = result.group(1)
                                        if not is_number(redshift) or float(redshift) > 100.:
                                            redshift = ''

                                    result = re.search('publish=(.*?)&amp;', trstr)
                                    bibcode = ''
                                    if result:
                                        bibcode = unescape(urllib.parse.unquote(urllib.parse.unquote(result.group(1))).split('/')[-1])

                                    if not bibcode:
                                        biblink = tr.find('a', {'title': 'Link to NASA ADS'})
                                        if biblink:
                                            bibcode = biblink.contents[0]

                                    if name.startswith('sn'):
                                        name = 'SN' + name[2:]
                                    if name.startswith(('CSS', 'SSS', 'MLS')) and ':' not in name:
                                        name = name.replace('-', ':', 1)
                                    if name.startswith('MASTERJ'):
                                        name = name.replace('MASTERJ', 'MASTER OT J')
                                    if name.startswith('PSNJ'):
                                        name = name.replace('PSNJ', 'PSN J')
                                    name = get_preferred_name(events, name)
                                    if oldname and name != oldname:
                                        journal_events(tasks, args, events)
                                    oldname = name
                                    name = add_event(tasks, args, events, name)

                                    # print(name + ' ' + claimedtype + ' ' + epoch + ' ' + observer + ' ' + reducer + ' ' + specfile + ' ' + bibcode + ' ' + redshift)

                                    secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondarybibcode, secondary=True)
                                    add_quantity(events, name, 'alias', name, secondarysource)
                                    if bibcode:
                                        newbibcode = bibcode
                                        if bibcode in wiserepbibcorrectdict:
                                            newbibcode = wiserepbibcorrectdict[bibcode]
                                        if newbibcode:
                                            source = add_source(events, name, bibcode=unescape(newbibcode))
                                        else:
                                            source = add_source(events, name, refname=unescape(bibcode))
                                        sources = uniq_cdl([source, secondarysource])
                                    else:
                                        sources = secondarysource

                                    if claimedtype not in ['Other']:
                                        add_quantity(events, name, 'claimedtype', claimedtype, secondarysource)
                                    add_quantity(events, name, 'redshift', redshift, secondarysource)

                                    if not specpath:
                                        continue

                                    with open(specpath, 'r') as f:
                                        data = [x.split() for x in f]

                                        skipspec = False
                                        newdata = []
                                        oldval = ''
                                        for row in data:
                                            if row and '#' not in row[0]:
                                                if len(row) >= 2 and is_number(row[0]) and is_number(row[1]) and row[1] != oldval:
                                                    newdata.append(row)
                                                    oldval = row[1]

                                        if skipspec or not newdata:
                                            warnings.warn('Skipped adding spectrum file ' + specfile)
                                            continue

                                        data = [list(i) for i in zip(*newdata)]
                                        wavelengths = data[0]
                                        fluxes = data[1]
                                        errors = ''
                                        if len(data) == 3:
                                            errors = data[1]
                                        time = str(astrotime(epoch).mjd)

                                        if max([float(x) for x in fluxes]) < 1.0e-5:
                                            fluxunit = 'erg/s/cm^2/Angstrom'
                                        else:
                                            fluxunit = 'Uncalibrated'

                                        add_spectrum(events, name=name, waveunit='Angstrom', fluxunit=fluxunit, errors=errors, errorunit=fluxunit, wavelengths=wavelengths,
                                                     fluxes=fluxes, u_time='MJD', time=time, instrument=instrument, source=sources, observer=observer, reducer=reducer,
                                                     filename=specfile)
                                        wiserepcnt = wiserepcnt + 1

                                        if args.travis and wiserepcnt % TRAVIS_QUERY_LIMIT == 0:
                                            break

                        tprint('Unadded files: ' + str(len(lfiles) - 1) + "/" + str(len(files)-1))
                        tprint('WISeREP spectrum count: ' + str(wiserepcnt))
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'cfaspectra'):
            # Ia spectra
            oldname = ''
            file_names = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIa')))[1]
            for name in pbar_strings(file_names, current_task):
                fullpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIa/') + name
                origname = name
                if name.startswith('sn') and is_number(name[2:6]):
                    name = 'SN' + name[2:]
                if name.startswith('snf') and is_number(name[3:7]):
                    name = 'SNF' + name[3:]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                reference = 'CfA Supernova Archive'
                refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
                source = add_source(events, name, refname=reference, url=refurl, secondary=True, acknowledgment=ACKN_CFA)
                add_quantity(events, name, 'alias', name, source)
                for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
                    filename = os.path.basename(fname)
                    fileparts = filename.split('-')
                    if origname.startswith('sn') and is_number(origname[2:6]):
                        year = fileparts[1][:4]
                        month = fileparts[1][4:6]
                        day = fileparts[1][6:]
                        instrument = fileparts[2].split('.')[0]
                    else:
                        year = fileparts[2][:4]
                        month = fileparts[2][4:6]
                        day = fileparts[2][6:]
                        instrument = fileparts[3].split('.')[0]
                    time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                    f = open(fname, 'r')
                    data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                    data = [list(i) for i in zip(*data)]
                    wavelengths = data[0]
                    fluxes = data[1]
                    errors = data[2]
                    sources = uniq_cdl([source, add_source(events, name, bibcode='2012AJ....143..126B'), add_source(events, name, bibcode='2008AJ....135.1598M')])
                    add_spectrum(
                        name=name, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', filename=filename,
                        wavelengths=wavelengths, fluxes=fluxes, u_time='MJD' if time else '', time=time, instrument=instrument,
                        errorunit='ergs/s/cm^2/Angstrom', errors=errors, source=sources, dereddened=False, deredshifted=False)
                    if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                        break
            journal_events(tasks, args, events)

            # Ibc spectra
            oldname = ''
            file_names = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIbc')))[1]
            for name in pbar(file_names, current_task):
                fullpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIbc/') + name
                if name.startswith('sn') and is_number(name[2:6]):
                    name = 'SN' + name[2:]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                reference = 'CfA Supernova Archive'
                refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
                source = add_source(events, name, refname=reference, url=refurl, secondary=True, acknowledgment=ACKN_CFA)
                add_quantity(events, name, 'alias', name, source)
                for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
                    filename = os.path.basename(fname)
                    fileparts = filename.split('-')
                    instrument = ''
                    year = fileparts[1][:4]
                    month = fileparts[1][4:6]
                    day = fileparts[1][6:].split('.')[0]
                    if len(fileparts) > 2:
                        instrument = fileparts[-1].split('.')[0]
                    time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                    f = open(fname, 'r')
                    data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                    data = [list(i) for i in zip(*data)]
                    wavelengths = data[0]
                    fluxes = data[1]
                    sources = uniq_cdl([source, add_source(events, name, bibcode='2014AJ....147...99M')])
                    add_spectrum(
                        name=name, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths, filename=filename,
                        fluxes=fluxes, u_time='MJD' if time else '', time=time, instrument=instrument, source=sources,
                        dereddened=False, deredshifted=False)
                    if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                        break
            journal_events(tasks, args, events)

            # Other spectra
            oldname = ''
            file_names = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_Extra')))[1]
            for name in pbar_strings(file_names, current_task):
                fullpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_Extra/') + name
                if name.startswith('sn') and is_number(name[2:6]):
                    name = 'SN' + name[2:]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                reference = 'CfA Supernova Archive'
                refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
                source = add_source(events, name, refname=reference, url=refurl, secondary=True, acknowledgment=ACKN_CFA)
                add_quantity(events, name, 'alias', name, source)
                for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
                    if not os.path.isfile(fname):
                        continue
                    filename = os.path.basename(fname)
                    if ((not filename.startswith('sn') or not filename.endswith('flm') or
                         any(x in filename for x in ['-interp', '-z', '-dered', '-obj', '-gal']))):
                        continue
                    fileparts = filename.split('.')[0].split('-')
                    instrument = ''
                    time = ''
                    if len(fileparts) > 1:
                        year = fileparts[1][:4]
                        month = fileparts[1][4:6]
                        day = fileparts[1][6:]
                        if is_number(year) and is_number(month) and is_number(day):
                            if len(fileparts) > 2:
                                instrument = fileparts[-1]
                            time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                    f = open(fname, 'r')
                    data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                    data = [list(i) for i in zip(*data)]
                    wavelengths = data[0]
                    fluxes = [str(Decimal(x)*Decimal(1.0e-15)) for x in data[1]]
                    add_spectrum(
                        name=name, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths, filename=filename,
                        fluxes=fluxes, u_time='MJD' if time else '', time=time, instrument=instrument, source=source,
                        dereddened=False, deredshifted=False)
                    if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                        break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snlsspectra'):
            result = Vizier.get_catalogs('J/A+A/507/85/table1')
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            datedict = {}
            for row in table:
                datedict['SNLS-' + row['SN']] = str(astrotime(row['Date']).mjd)

            oldname = ''
            file_names = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNLS/*'))
            for fi, fname in enumerate(pbar_strings(file_names, current_task=current_task)):
                filename = os.path.basename(fname)
                fileparts = filename.split('_')
                name = 'SNLS-' + fileparts[1]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2009A&A...507...85B')
                add_quantity(events, name, 'alias', name, source)

                add_quantity(events, name, 'discoverdate', '20' + fileparts[1][:2], source)

                f = open(fname, 'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                specdata = []
                for r, row in enumerate(data):
                    if row[0] == '@TELESCOPE':
                        telescope = row[1].strip()
                    elif row[0] == '@REDSHIFT':
                        add_quantity(events, name, 'redshift', row[1].strip(), source)
                    if r < 14:
                        continue
                    specdata.append(list(filter(None, [x.strip(' \t') for x in row])))
                specdata = [list(i) for i in zip(*specdata)]
                wavelengths = specdata[1]

                fluxes = [pretty_num(float(x)*1.e-16, sig=get_sig_digits(x)) for x in specdata[2]]
                errors = [pretty_num(float(x)*1.e-16, sig=get_sig_digits(x)) for x in specdata[3]]

                add_spectrum(
                    name=name, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths,
                    fluxes=fluxes, u_time='MJD' if name in datedict else '',
                    time=datedict[name] if name in datedict else '', telescope=telescope, source=source,
                    filename=filename)
                if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                    break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'cspspectra'):
            oldname = ''
            file_names = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CSP/*')
            for fi, fname in enumerate(pbar_strings(file_names, current_task=current_task)):
                filename = os.path.basename(fname)
                sfile = filename.split('.')
                if sfile[1] == 'txt':
                    continue
                sfile = sfile[0]
                fileparts = sfile.split('_')
                name = 'SN20' + fileparts[0][2:]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                telescope = fileparts[-2]
                instrument = fileparts[-1]
                source = add_source(events, name, bibcode='2013ApJ...773...53F')
                add_quantity(events, name, 'alias', name, source)

                f = open(fname, 'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                specdata = []
                for r, row in enumerate(data):
                    if row[0] == '#JDate_of_observation:':
                        jd = row[1].strip()
                        time = str(jd_to_mjd(Decimal(jd)))
                    elif row[0] == '#Redshift:':
                        add_quantity(events, name, 'redshift', row[1].strip(), source)
                    if r < 7:
                        continue
                    specdata.append(list(filter(None, [x.strip(' ') for x in row])))
                specdata = [list(i) for i in zip(*specdata)]
                wavelengths = specdata[0]
                fluxes = specdata[1]

                add_spectrum(
                    name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths,
                    fluxes=fluxes, telescope=telescope, instrument=instrument, source=source, deredshifted=True, filename=filename)
                if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                    break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ucbspectra'):
            secondaryreference = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
            secondaryrefurl = 'http://heracles.astro.berkeley.edu/sndb/info'
            secondaryrefbib = '2012MNRAS.425.1789S'
            ucbspectracnt = 0

            jsontxt = load_cached_url(
                args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubspec',
                os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
            if not jsontxt:
                continue

            spectra = json.loads(jsontxt)
            spectra = sorted(spectra, key=lambda k: k['ObjName'])
            oldname = ''
            for spectrum in pbar(spectra, desc=current_task):
                name = spectrum['ObjName']
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)

                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondaryrefbib, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                sources = [secondarysource]
                if spectrum['Reference']:
                    sources += [add_source(events, name, bibcode=spectrum['Reference'])]
                sources = uniq_cdl(sources)

                if spectrum['Type'] and spectrum['Type'].strip() != 'NoMatch':
                    for ct in spectrum['Type'].strip().split(','):
                        add_quantity(events, name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
                if spectrum['DiscDate']:
                    add_quantity(events, name, 'discoverdate', spectrum['DiscDate'].replace('-', '/'), sources)
                if spectrum['HostName']:
                    add_quantity(events, name, 'host', urllib.parse.unquote(spectrum['HostName']).replace('*', ''), sources)
                if spectrum['UT_Date']:
                    epoch = str(spectrum['UT_Date'])
                    year = epoch[:4]
                    month = epoch[4:6]
                    day = epoch[6:]
                    sig = get_sig_digits(day) + 5
                    mjd = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig=sig)
                filename = spectrum['Filename'] if spectrum['Filename'] else ''
                instrument = spectrum['Instrument'] if spectrum['Instrument'] else ''
                reducer = spectrum['Reducer'] if spectrum['Reducer'] else ''
                observer = spectrum['Observer'] if spectrum['Observer'] else ''
                snr = str(spectrum['SNR']) if spectrum['SNR'] else ''

                if not filename:
                    raise(ValueError('Filename not found for SNDB spectrum!'))
                if not spectrum['SpecID']:
                    raise(ValueError('ID not found for SNDB spectrum!'))

                filepath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/') + filename
                if archived_task('ucbspectra', tasks) and os.path.isfile(filepath):
                    with open(filepath, 'r') as f:
                        spectxt = f.read()
                else:
                    session = requests.Session()
                    response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=ds:' + str(spectrum['SpecID']))
                    spectxt = response.text
                    with open(filepath, 'w') as f:
                        f.write(spectxt)

                specdata = list(csv.reader(spectxt.splitlines(), delimiter=' ', skipinitialspace=True))
                startrow = 0
                for row in specdata:
                    if row[0][0] == '#':
                        startrow += 1
                    else:
                        break
                specdata = specdata[startrow:]

                haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]

                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]

                if not list(filter(None, errors)):
                    errors = ''

                add_spectrum(
                    name=name, u_time='MJD', time=mjd, waveunit='Angstrom', fluxunit='Uncalibrated',
                    wavelengths=wavelengths, filename=filename, fluxes=fluxes, errors=errors, errorunit='Uncalibrated',
                    instrument=instrument, source=sources, snr=snr, observer=observer, reducer=reducer,
                    deredshifted=('-noz' in filename))
                ucbspectracnt = ucbspectracnt + 1
                if args.travis and ucbspectracnt >= TRAVIS_QUERY_LIMIT:
                    break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'suspectspectra'):
            with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/sources.json'), 'r') as f:
                sourcedict = json.loads(f.read())

            with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/filename-changes.txt'), 'r') as f:
                rows = f.readlines()
                changedict = {}
                for row in rows:
                    if not row.strip() or row[0] == "#":
                        continue
                    items = row.strip().split(' ')
                    changedict[items[1]] = items[0]

            suspectcnt = 0
            folders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect')))[1]
            for folder in pbar(folders, current_task):
                eventfolders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/')+folder))[1]
                oldname = ''
                for eventfolder in pbar(eventfolders, current_task):
                    name = eventfolder
                    if is_number(name[:4]):
                        name = 'SN' + name
                    name = get_preferred_name(events, name)
                    if oldname and name != oldname:
                        journal_events(tasks, args, events)
                    oldname = name
                    name = add_event(tasks, args, events, name)
                    secondaryreference = 'SUSPECT'
                    secondaryrefurl = 'https://www.nhn.ou.edu/~suspect/'
                    secondarybibcode = '2001AAS...199.8408R'
                    secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondarybibcode, secondary=True)
                    add_quantity(events, name, 'alias', name, secondarysource)
                    eventspectra = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/')+folder+'/'+eventfolder))[2]
                    for spectrum in eventspectra:
                        sources = [secondarysource]
                        bibcode = ''
                        if spectrum in changedict:
                            specalias = changedict[spectrum]
                        else:
                            specalias = spectrum
                        if specalias in sourcedict:
                            bibcode = sourcedict[specalias]
                        elif name in sourcedict:
                            bibcode = sourcedict[name]
                        if bibcode:
                            source = add_source(events, name, bibcode=unescape(bibcode))
                            sources += source
                        sources = uniq_cdl(sources)

                        date = spectrum.split('_')[1]
                        year = date[:4]
                        month = date[4:6]
                        day = date[6:]
                        sig = get_sig_digits(day) + 5
                        time = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig=sig)

                        with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/'+folder+'/'+eventfolder+'/'+spectrum)) as f:
                            specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                            specdata = list(filter(None, specdata))
                            newspec = []
                            oldval = ''
                            for row in specdata:
                                if row[1] == oldval:
                                    continue
                                newspec.append(row)
                                oldval = row[1]
                            specdata = newspec
                        haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                        specdata = [list(i) for i in zip(*specdata)]

                        wavelengths = specdata[0]
                        fluxes = specdata[1]
                        errors = ''
                        if haserrors:
                            errors = specdata[2]

                        add_spectrum(
                            name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='Uncalibrated', wavelengths=wavelengths,
                            fluxes=fluxes, errors=errors, errorunit='Uncalibrated', source=sources, filename=spectrum)
                        suspectcnt = suspectcnt + 1
                        if args.travis and suspectcnt % TRAVIS_QUERY_LIMIT == 0:
                            break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snfspectra'):
            eventfolders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory')))[1]
            bibcodes = {'SN2005gj': '2006ApJ...650..510A', 'SN2006D': '2007ApJ...654L..53T', 'SN2007if': '2010ApJ...713.1073S', 'SN2011fe': '2013A&A...554A..27P'}
            oldname = ''
            snfcnt = 0
            for eventfolder in eventfolders:
                name = eventfolder
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                secondaryreference = 'Nearby Supernova Factory'
                secondaryrefurl = 'http://snfactory.lbl.gov/'
                secondarybibcode = '2002SPIE.4836...61A'
                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondarybibcode, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                bibcode = bibcodes[name]
                source = add_source(events, name, bibcode=bibcode)
                sources = uniq_cdl([source, secondarysource])
                eventspectra = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory/')+eventfolder+'/*.dat')
                for spectrum in eventspectra:
                    filename = os.path.basename(spectrum)
                    with open(spectrum) as f:
                        specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                    specdata = list(filter(None, specdata))
                    newspec = []
                    time = ''
                    telescope = ''
                    instrument = ''
                    observer = ''
                    observatory = ''
                    if 'Keck_20060202_R' in spectrum:
                        time = '53768.23469'
                    elif 'Spectrum05_276' in spectrum:
                        time = pretty_num(astrotime('2005-10-03').mjd, sig=5)
                    elif 'Spectrum05_329' in spectrum:
                        time = pretty_num(astrotime('2005-11-25').mjd, sig=5)
                    elif 'Spectrum05_336' in spectrum:
                        time = pretty_num(astrotime('2005-12-02').mjd, sig=5)
                    for row in specdata:
                        if row[0][0] == '#':
                            joinrow = (' '.join(row)).split('=')
                            if len(joinrow) < 2:
                                continue
                            field = joinrow[0].strip('# ')
                            value = joinrow[1].split('/')[0].strip('\' ')
                            if not time:
                                if field == 'JD':
                                    time = str(jd_to_mjd(Decimal(value)))
                                elif field == 'MJD':
                                    time = value
                                elif field == 'MJD-OBS':
                                    time = value
                            if field == 'OBSERVER':
                                observer = value.capitalize()
                            if field == 'OBSERVAT':
                                observatory = value.capitalize()
                            if field == 'TELESCOP':
                                telescope = value.capitalize()
                            if field == 'INSTRUME':
                                instrument = value.capitalize()
                        else:
                            newspec.append(row)
                    if not time:
                        raise(ValueError('Time missing from spectrum.'))
                    specdata = newspec
                    haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                    specdata = [list(i) for i in zip(*specdata)]

                    wavelengths = specdata[0]
                    fluxes = specdata[1]
                    errors = ''
                    if haserrors:
                        errors = specdata[2]

                    add_spectrum(
                        name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom',
                        wavelengths=wavelengths, fluxes=fluxes, errors=errors, observer=observer, observatory=observatory,
                        telescope=telescope, instrument=instrument,
                        errorunit=('Variance' if name == 'SN2011fe' else 'erg/s/cm^2/Angstrom'), source=sources, filename=filename)
                    snfcnt = snfcnt + 1
                    if args.travis and snfcnt % TRAVIS_QUERY_LIMIT == 0:
                        break
            journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'superfitspectra'):
            sfdirs = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'superfit/*'))
            for sfdir in pbar(sfdirs, desc=current_task):
                sffiles = sorted(glob(sfdir + '/*.dat'))
                lastname = ''
                oldname = ''
                for sffile in pbar(sffiles, desc=current_task):
                    basename = os.path.basename(sffile)
                    name = basename.split('.')[0]
                    if name.startswith('sn'):
                        name = 'SN' + name[2:]
                        if len(name) == 7:
                            name = name[:6] + name[6].upper()
                    elif name.startswith('ptf'):
                        name = 'PTF' + name[3:]

                    if 'theory' in name:
                        continue
                    if event_exists(events, name):
                        prefname = get_preferred_name(events, name)
                        if 'spectra' in events[prefname] and lastname != prefname:
                            continue
                    if oldname and name != oldname:
                        journal_events(tasks, args, events)
                    oldname = name
                    name = add_event(tasks, args, events, name)
                    epoch = basename.split('.')[1]
                    (mldt, mlmag, mlband, mlsource) = get_max_light(events, name)
                    if mldt:
                        epoff = Decimal(0.0) if epoch == 'max' else (Decimal(epoch[1:]) if epoch[0] == 'p' else -Decimal(epoch[1:]))
                    else:
                        epoff = ''

                    source = add_source(events, name, refname='Superfit', url='http://www.dahowell.com/superfit.html', secondary=True)
                    add_quantity(events, name, 'alias', name, source)

                    with open(sffile) as f:
                        rows = f.read().splitlines()
                    specdata = []
                    for row in rows:
                        if row.strip():
                            specdata.append(list(filter(None, re.split('\t+|\s+', row, maxsplit=0))))
                    specdata = [[x.replace('D', 'E') for x in list(i)] for i in zip(*specdata)]
                    wavelengths = specdata[0]
                    fluxes = specdata[1]

                    mlmjd = str(Decimal(astrotime('-'.join([str(mldt.year), str(mldt.month), str(mldt.day)])).mjd) + epoff) if (epoff != '') else ''
                    add_spectrum(
                        name, u_time='MJD' if mlmjd else '', time=mlmjd, waveunit='Angstrom', fluxunit='Uncalibrated',
                        wavelengths=wavelengths, fluxes=fluxes, source=source)

                    lastname = name
                journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'mergeduplicates'):
            if args.update and not len(events):
                tprint('No sources changed, event files unchanged in update.')
                sys.exit(1)
            merge_duplicates(tasks, args, events)

        if do_task(tasks, args, task, 'setprefnames'):
            set_preferred_names(tasks, args, events)

    files = repo_file_list()

    bibauthor_dict = get_bibauthor_dict()
    extinctions_dict = get_extinctions_dict()

    for fi in pbar(files, 'Sanitizing and deriving quantities for events'):
        events = OrderedDict()
        name = os.path.basename(os.path.splitext(fi)[0]).replace('.json', '')
        name = add_event(tasks, args, events, args, name, loadifempty=False)
        events, extinctions_dict, bibauthor_dict = derive_and_sanitize(
            tasks, args, events, extinctions_dict, bibauthor_dict, nedd_dict)
        if has_task(tasks, args, 'writeevents'):
            write_all_events(events, args, empty=True, gz=True, bury=True)

    jsonstring = json.dumps(bibauthor_dict, indent='\t', separators=(',', ':'), ensure_ascii=False)
    with codecs.open('../bibauthors.json', 'w', encoding='utf8') as f:
        f.write(jsonstring)
    jsonstring = json.dumps(extinctions_dict, indent='\t', separators=(',', ':'), ensure_ascii=False)
    with codecs.open('../extinctions.json', 'w', encoding='utf8') as f:
        f.write(jsonstring)

    print('Memory used (MBs on Mac, GBs on Linux): ' + '{:,}'.format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))

    sys.exit(0)


def do_internal(events, args, tasks):
    """Load events from files in the 'internal' repository, and save them.
    """
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    for datafile in pbar_strings(files, desc=current_task):
        if args.update:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False, append=True):
                raise IOError('Failed to find specified file.')
        else:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False):
                raise IOError('Failed to find specified file.')

    events = journal_events(tasks, args, events)
    return events


def do_external_radio(events, args, tasks):
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as f:
            for li, line in enumerate([x.strip() for x in f.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    radiosourcedict[line.split()[0]] = add_source(events, name, bibcode=line.split()[-1])
                elif li in [x + len(radiosourcedict) for x in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    add_photometry(
                        events, name, time=cols[0], frequency=cols[2], u_frequency='GHz', fluxdensity=cols[3],
                        e_fluxdensity=cols[4], u_fluxdensity='µJy', instrument=cols[5], source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


def do_external_xray(events, args, tasks):
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        with open(datafile, 'r') as f:
            for li, line in enumerate(f.read().splitlines()):
                if li == 0:
                    source = add_source(events, name, bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    add_photometry(
                        events, name, time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4], flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/s/cm^2',
                        photonindex=cols[15], instrument=cols[17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


def do_simbad(events, args, tasks):
    Simbad.list_votable_fields()
    customSimbad = Simbad()
    customSimbad.add_votable_fields('otype', 'id(opt)')
    result = customSimbad.query_object('SN 20[0-9][0-9]*', wildcard=True)
    for r, row in enumerate(result):
        if row['OTYPE'].decode() != 'SN':
            continue
        name = row['MAIN_ID'].decode()
        aliases = Simbad.query_objectids(name)
        print(aliases)
        if name[:3] == 'SN ':
            name = 'SN' + name[3:]
        if name[:2] == 'SN' and is_number(name[2:]):
            name = name + 'A'
        name = add_event(tasks, args, events, name)
    events = journal_events(tasks, args, events)
    return events


def do_vizier(events, args, tasks):



def load_args():
    parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
    parser.add_argument('--update', '-u',       dest='update',      help='Only update catalog using live sources.',    default=False, action='store_true')
    parser.add_argument('--verbose', '-v',      dest='verbose',     help='Print more messages to the screen.',         default=False, action='store_true')
    parser.add_argument('--refresh', '-r',      dest='refresh',     help='Ignore most task caches.',                   default=False, action='store_true')
    parser.add_argument('--full-refresh', '-f', dest='fullrefresh', help='Ignore all task caches.',                    default=False, action='store_true')
    parser.add_argument('--archived', '-a',     dest='archived',    help='Always use task caches.',                    default=False, action='store_true')
    parser.add_argument('--travis', '-tr',      dest='travis',      help='Run import script in test mode for Travis.', default=False, action='store_true')
    parser.add_argument('--refreshlist', '-rl', dest='refreshlist', help='Comma-delimited list of caches to clear.',   default='')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    import_main()
