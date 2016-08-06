"""
"""
import csv
import filecmp
import gzip
import json
import operator
import os
import re
import shutil
import urllib.parse
import urllib.request
import warnings
from collections import OrderedDict
from copy import deepcopy
from glob import glob
from math import ceil, isnan, pi
from statistics import mean

import inflect
import numpy
from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.time import Time as astrotime
from bokeh.embed import file_html
from bokeh.layouts import row as bokehrow
from bokeh.layouts import column, layout
from bokeh.models import (ColumnDataSource, CustomJS, DatetimeAxis, HoverTool,
                          LinearAxis, Range1d, Slider)
from bokeh.models.widgets import Select
from bokeh.plotting import Figure, reset_output
from bokeh.resources import CDN
from bs4 import BeautifulSoup

from astrocats.catalog.utils import (bandaliasf, bandcolorf,
                                     bandgroupf, bandshortaliasf, bandwavef,
                                     get_sig_digits,
                                     is_number, pretty_num, radiocolorf,
                                     round_sig, tprint, tq, xraycolorf)
from astrocats.supernovae.scripts.events import (get_event_filename,
                                                 get_event_text)
from cdecimal import Decimal

from .utils import touch, label_format, get_first_kind, \
    get_first_value, md5file

from . import webplot_photometry, webplot_spectra, webplot_radio, webplot_xray

from .constants import TRAVIS_LIMIT, RADIO_SIGMA, GOOGLE_PING_URL, SNE_LINK_DIR, DEF_COLORS, \
    COLUMN_KEYS, EVENT_IGNORE_KEY, HEADER, EVENT_PAGE_HEADER, DEF_TITLES, SNE_PAGES, \
    SITEMAP_TEMPLATE, TOOLS_LIST, SAFE_FILES


def main(astro_catalog):
    args = astro_catalog.args
    paths = astro_catalog.paths
    log = astro_catalog.log

    infl = inflect.engine()
    infl.defnoun("spectrum", "spectra")
    testsuffix = '.test' if args.test else ''

    catalog = OrderedDict()
    catalogcopy = OrderedDict()

    sourcedict = {}
    lcspye = []
    lcspno = []
    lconly = []
    sponly = []
    hasalc = []
    hasasp = []
    totalphoto = 0
    totalspectra = 0

    # Host Images File
    host_img_dict = load_dict_file(paths.host_imgs_file, log)

    # MD5 Hash File
    md5_dict = load_dict_file(paths.md5_file, log)

    files = astro_catalog.paths.get_repo_output_file_list(
        normal=(not args.boneyard), bones=args.boneyard)
    files = sorted(files, key=lambda ss: ss.lower())
    for fcnt, eventfile in enumerate(tq(files)):
        event_file_name = os.path.splitext(os.path.basename(eventfile))[0]
        event_file_name = event_file_name.replace('.json', '')
        log.debug("Event file: '{}'".format(eventfile))
        if args.eventlist and event_file_name not in args.eventlist:
            continue

        if args.travis and fcnt >= TRAVIS_LIMIT:
            log.warning("Hit Travis Limit (`TRAVIS_LIMIT` = '{}').".format(TRAVIS_LIMIT))
            break

        entry_changed = False
        checksum = md5file(eventfile)
        if eventfile not in md5_dict or md5_dict[eventfile] != checksum:
            entry_changed = True
            md5_dict[eventfile] = checksum

        file_text = get_event_text(eventfile)

        catalog.update(json.loads(file_text, object_pairs_hook=OrderedDict))
        entry = next(reversed(catalog))

        event_name = entry

        # FIX: why would this be different than the above similar check?
        if args.eventlist and event_name not in args.eventlist:
            continue

        # tprint(eventfile + ' [' + checksum + ']')
        log.info("'{}': [{}]".format(event_name, checksum))

        # FIX: 'internal' is a list, for now just select the first one...
        internal_repo = astro_catalog.paths.repos_dict['internal'][0]
        internal_file = os.path.join(internal_repo, event_file_name + ".json")
        if os.path.isfile(internal_file):
            catalog[entry]['download'] = 'e'
        else:
            catalog[entry]['download'] = ''

        if 'discoverdate' in catalog[entry]:
            for d, date in enumerate(catalog[entry]['discoverdate']):
                catalog[entry]['discoverdate'][d]['value'] = \
                    catalog[entry]['discoverdate'][d]['value'].split('.')[0]
        if 'maxdate' in catalog[entry]:
            for d, date in enumerate(catalog[entry]['maxdate']):
                catalog[entry]['maxdate'][d]['value'] = \
                    catalog[entry]['maxdate'][d]['value'].split('.')[0]

        hostmag = ''
        hosterr = ''
        if 'photometry' in catalog[entry]:
            for photo in catalog[entry]['photometry']:
                if 'host' in photo and ('upperlimit' not in photo or
                                        not photo['upperlimit']):
                    hostmag = float(photo['magnitude'])
                    hosterr = float(photo[
                        'e_magnitude']) if 'e_magnitude' in photo else 0.0

            # Delete the host magnitudes so they are not plotted as points
            catalog[entry]['photometry'][:] = [
                x for x in catalog[entry]['photometry'] if 'host' not in x
            ]

        # FIX: if 'photometry' is not in entry, add empty list... will make everything easier!

        photoavail = 'photometry' in catalog[entry] and any(
            ['magnitude' in x for x in catalog[entry]['photometry']])
        radioavail = 'photometry' in catalog[entry] and any(
            ['fluxdensity' in x for x in catalog[entry]['photometry']])
        xrayavail = 'photometry' in catalog[entry] and any(
            ['counts' in x and 'magnitude' not in x
             for x in catalog[entry]['photometry']])
        spectraavail = 'spectra' in catalog[entry]

        # Must be two sigma above host magnitude, if host magnitude known, to add
        # to phot count.
        numphoto = len(
            [x for x in catalog[entry]['photometry']
             if 'upperlimit' not in x and 'magnitude' in x and
             (not hostmag or 'includeshost' not in x or float(x['magnitude']) <= (
                 hostmag - 2.0 * hosterr))]) if photoavail else 0
        numradio = len([x for x in catalog[entry]['photometry']
                        if 'upperlimit' not in x and 'fluxdensity' in x and
                        (not x['e_fluxdensity'] or float(x['fluxdensity']) >
                         RADIO_SIGMA * float(x['e_fluxdensity'])
                         ) and (not hostmag or 'includeshost' not in x or float(x[
                             'magnitude']) <= (hostmag - 2.0 * hosterr))
                        ]) if photoavail else 0
        numxray = len(
            [x for x in catalog[entry]['photometry']
             if 'upperlimit' not in x and 'counts' in x and
             (not hostmag or 'includeshost' not in x or float(x['magnitude']) <= (
                 hostmag - 2.0 * hosterr))]) if photoavail else 0
        numspectra = len(catalog[entry]['spectra']) if spectraavail else 0

        redshiftfactor = (1.0 / (
            1.0 + float(catalog[entry]['redshift'][0]['value']))) if (
                'redshift' in catalog[entry]) else 1.0
        dayframe = 'Rest frame days' if 'redshift' in catalog[
            entry] else 'Observer frame days'

        mjdmax = ''
        if 'maxdate' in catalog[entry]:
            datestr = catalog[entry]['maxdate'][0]['value']
            datesplit = datestr.split('/')
            if len(datesplit) < 2:
                datestr += "/01"
            if len(datesplit) < 3:
                datestr += "/01"
            try:
                mjdmax = astrotime(datestr.replace('/', '-')).mjd
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                pass

        minphotoep = ''
        maxphotoep = ''
        if mjdmax:
            photoeps = [(Decimal(x['time']) - Decimal(mjdmax + 0.5)) *
                        Decimal(redshiftfactor)
                        for x in catalog[entry]['photometry']
                        if 'upperlimit' not in x and 'includeshost' not in x and
                        'magnitude' in x and 'time' in x] if photoavail else []
            if photoeps:
                minphotoep = pretty_num(float(min(photoeps)), sig=3)
                maxphotoep = pretty_num(float(max(photoeps)), sig=3)

        minspectraep = ''
        maxspectraep = ''
        if mjdmax:
            spectraeps = ([(Decimal(x['time']) - Decimal(mjdmax + 0.5)) *
                           Decimal(redshiftfactor)
                           for x in catalog[entry]['spectra'] if 'time' in x]
                          if spectraavail else [])
            if spectraeps:
                minspectraep = pretty_num(float(min(spectraeps)), sig=3)
                maxspectraep = pretty_num(float(max(spectraeps)), sig=3)

        catalog[entry]['numphoto'] = numphoto
        catalog[entry]['numradio'] = numradio
        catalog[entry]['numxray'] = numxray
        catalog[entry]['numspectra'] = numspectra

        distancemod = 0.0
        if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
            distancemod = float(get_first_value(catalog, entry, 'maxappmag')) - \
                float(get_first_value(catalog, entry, 'maxabsmag'))

        plotlink = "sne/" + event_file_name + "/"
        if photoavail:
            catalog[entry]['photolink'] = (str(numphoto) + (
                (',' + minphotoep + ',' + maxphotoep) if
                (minphotoep and maxphotoep and minphotoep != maxphotoep) else ''))
        if radioavail:
            catalog[entry]['radiolink'] = str(numradio)
        if xrayavail:
            catalog[entry]['xraylink'] = str(numxray)
        if spectraavail:
            catalog[entry]['spectralink'] = (str(numspectra) + (
                (',' + minspectraep + ',' + maxspectraep)
                if (minspectraep and maxspectraep and minspectraep != maxspectraep
                    ) else ''))

        prange = list(range(len(catalog[entry][
            'photometry']))) if 'photometry' in catalog[entry] else []

        instrulist = {catalog[entry]['photometry'][x]['instrument']
                      if 'instrument' in catalog[entry]['photometry'][x]
                      else None
                      for x in prange}
        instrulist = sorted([_f for _f in list(instrulist) if _f])
        if len(instrulist) > 0:
            instruments = ''
            for i, instru in enumerate(instrulist):
                instruments += instru
                bandlist = sorted(
                    [_f
                     for _f in list({bandshortaliasf(catalog[entry]['photometry'][
                         x]['band'] if 'band' in catalog[entry]['photometry'][
                             x] else '') if 'instrument' in catalog[entry][
                                 'photometry'][x] and catalog[entry]['photometry'][
                                     x]['instrument'] == instru else ""
                                     for x in prange}) if _f],
                    key=lambda y: (bandwavef(y), y))
                if bandlist:
                    instruments += ' (' + ", ".join(bandlist) + ')'
                if i < len(instrulist) - 1:
                    instruments += ', '

            # Now add bands without attached instrument
            obandlist = sorted(
                [_f
                 for _f in list({bandshortaliasf(catalog[entry]['photometry'][x][
                     'band'] if 'band' in catalog[entry]['photometry'][x] else '')
                                 if 'instrument' not in catalog[entry][
                                     'photometry'][x] else ""
                                 for x in prange}) if _f],
                key=lambda y: (bandwavef(y), y))
            if obandlist:
                instruments += ", " + ", ".join(obandlist)
            catalog[entry]['instruments'] = instruments
        else:
            bandlist = sorted(
                [_f
                 for _f in list({bandshortaliasf(catalog[entry]['photometry'][x][
                     'band'] if 'band' in catalog[entry]['photometry'][x] else '')
                                 for x in prange}) if _f],
                key=lambda y: (bandwavef(y), y))
            if len(bandlist) > 0:
                catalog[entry]['instruments'] = ", ".join(bandlist)

        # Check file modification times before constructing .html files, which is
        # expensive
        dohtml = True
        if not args.forcehtml:
            if os.path.isfile(paths.output_html + event_file_name + ".html"):
                if not entry_changed:
                    dohtml = False

        # Copy JSON files up a directory if they've changed
        if dohtml:
            shutil.copy2(eventfile, paths.output_json + os.path.basename(eventfile))

        if (photoavail or radioavail or xrayavail) and dohtml and args.writehtml:
            phototime = [
                (mean([float(y) for y in x['time']])
                 if isinstance(x['time'], list) else float(x['time']))
                for x in catalog[entry]['photometry']
                if any([y in x for y in ['fluxdensity', 'magnitude', 'flux']])
            ]
            phototimelowererrs = [
                float(x['e_lower_time'])
                if ('e_lower_time' in x and 'e_upper_time' in x) else
                (float(x['e_time']) if 'e_time' in x else 0.)
                for x in catalog[entry]['photometry']
                if any([y in x for y in ['fluxdensity', 'magnitude', 'flux']])
            ]
            phototimeuppererrs = [
                float(x['e_upper_time'])
                if ('e_lower_time' in x and 'e_upper_time' in x) in x else
                (float(x['e_time']) if 'e_time' in x else 0.)
                for x in catalog[entry]['photometry']
                if any([y in x for y in ['fluxdensity', 'magnitude', 'flux']])
            ]

            x_buffer = 0.1 * (
                max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0

            min_x_range = -0.5 * x_buffer + \
                min([x - y for x, y in list(zip(phototime, phototimeuppererrs))])
            max_x_range = 2.0 * x_buffer + \
                max([x + y for x, y in list(zip(phototime, phototimelowererrs))])

        if photoavail and dohtml and args.writehtml:
            webplot_photometry.plot_photo(catalog, entry)

        if spectraavail and dohtml and args.writehtml:
            webplot_spectra.plot_spectra(catalog, entry)

        if radioavail and dohtml and args.writehtml:
            webplot_radio.plot_radio(catalog, entry)

        if xrayavail and dohtml and args.writehtml:
            phototime = [(mean([float(y) for y in x['time']])
                          if isinstance(x['time'], list) else float(x['time']))
                         for x in catalog[entry]['photometry'] if 'flux' in x]
            phototimelowererrs = [float(x['e_lower_time'])
                                  if ('e_lower_time' in x and 'e_upper_time' in x)
                                  else (float(x['e_time'])
                                        if 'e_time' in x else 0.)
                                  for x in catalog[entry]['photometry']
                                  if 'flux' in x]
            phototimeuppererrs = [
                float(x['e_upper_time'])
                if ('e_lower_time' in x and 'e_upper_time' in x) in x else
                (float(x['e_time']) if 'e_time' in x else 0.)
                for x in catalog[entry]['photometry'] if 'flux' in x
            ]
            photofl = [float(x['flux'])
                       if ('e_flux' not in x or
                           float(x['flux']) > RADIO_SIGMA * float(x['e_flux'])) else
                       round_sig(
                           RADIO_SIGMA * float(x['e_flux']),
                           sig=get_sig_digits(x['e_flux']))
                       for x in catalog[entry]['photometry'] if 'flux' in x]
            photoflerrs = [(float(x['e_flux']) if 'e_flux' in x else 0.)
                           for x in catalog[entry]['photometry'] if 'flux' in x]
            photoufl = [(x['u_flux'] if 'flux' in x else '')
                        for x in catalog[entry]['photometry'] if 'flux' in x]
            photoener = [((' - '.join([y.rstrip('.') for y in x['energy']])
                           if isinstance(x['energy'], list) else x['energy'])
                          if 'flux' in x else '')
                         for x in catalog[entry]['photometry'] if 'flux' in x]
            photouener = [(x['u_energy'] if 'flux' in x else '')
                          for x in catalog[entry]['photometry'] if 'flux' in x]
            photoinstru = [(x['instrument'] if 'instrument' in x else '')
                           for x in catalog[entry]['photometry'] if 'flux' in x]
            photosource = [', '.join(
                str(j)
                for j in sorted(
                    int(i)
                    for i in catalog[entry]['photometry'][x]['source'].split(',')))
                           for x, y in enumerate(catalog[entry]['photometry'])
                           if 'flux' in y]
            phototype = [
                (True if 'upperlimit' in x or
                 RADIO_SIGMA * float(x['e_flux']) >= float(x['flux']) else False)
                for x in catalog[entry]['photometry'] if 'flux' in x
            ]

            photoutime = catalog[entry]['photometry'][0][
                'u_time'] if 'u_time' in catalog[entry]['photometry'][0] else 'MJD'
            if distancemod:
                dist = (10.0**(1.0 + 0.2 * distancemod) * un.pc).cgs.value
                areacorr = 4.0 * pi * dist**2.0

            x_buffer = 0.1 * (
                max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0

            hastimeerrs = (len(list(filter(None, phototimelowererrs))) and
                           len(list(filter(None, phototimeuppererrs))))
            hasfl = len(list(filter(None, photofl)))
            hasflerrs = len(list(filter(None, photoflerrs)))
            yaxis = 'Flux'
            if not hasfl:
                yaxis = 'Counts'
                photofl = [float(x['counts'])
                           if ('e_counts' not in x or
                               float(x['counts']) > RADIO_SIGMA * float(x['e_counts'])) else
                           round_sig(
                               RADIO_SIGMA * float(x['e_counts']),
                               sig=get_sig_digits(x['e_counts']))
                           for x in catalog[entry]['photometry'] if 'counts' in x]
                photoflerrs = [(float(x['e_counts']) if 'e_counts' in x else 0.)
                               for x in catalog[entry]['photometry'] if 'counts' in x]
                photoufl = ['' for x in photofl]
                hasfl = len(list(filter(None, photofl)))
                hasflerrs = len(list(filter(None, photoflerrs)))
            tt = [
                ("Source ID(s)", "@src"),
                ("Epoch (" + photoutime + ")",
                 "@x{1.11}" + ("<sub>-@xle{1}</sub><sup>+@xue{1}</sup>"
                               if hastimeerrs else ""))
            ]
            if hasfl:
                tt += [(yaxis + " (" + photoufl[0].replace("ergs/s/cm^2",
                                                       "ergs s⁻¹ cm⁻²") + ")",
                        "@y" + ("&nbsp;±&nbsp;@err" if hasflerrs else ""))]
                if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
                    tt += [("Iso. Lum. (ergs s⁻¹)",
                            "@yabs" + ("&nbsp;±&nbsp;@abserr"
                                       if hasflerrs else ""))]
            if len(list(filter(None, photoener))):
                tt += [("Frequency (" + photouener[0] + ")", "@desc")]
            if len(list(filter(None, photoinstru))):
                tt += [("Instrument", "@instr")]
            hover = HoverTool(tooltips=tt)

            if photoavail:
                x_range = p1.x_range
            else:
                x_range = (min_x_range, max_x_range)
            min_y_range = min([x - y for x, y in list(zip(photofl, photoflerrs))])
            max_y_range = max([x + y for x, y in list(zip(photofl, photoflerrs))])
            [min_y_range, max_y_range] = [
                min_y_range - 0.1 * (max_y_range - min_y_range),
                max_y_range + 0.1 * (max_y_range - min_y_range)
            ]

            p4 = Figure(
                title='X-ray Observations of ' + event_name,
                active_drag='box_zoom',
                # sizing_mode = "scale_width",
                y_axis_label='Flux (ergs s⁻¹ cm⁻²)',
                tools=TOOLS_LIST,
                plot_width=485,
                plot_height=485,
                x_range=x_range,
                y_range=(min_y_range, max_y_range),
                toolbar_location='above',
                toolbar_sticky=False)
            p4.xaxis.axis_label_text_font = 'futura'
            p4.yaxis.axis_label_text_font = 'futura'
            p4.xaxis.major_label_text_font = 'futura'
            p4.yaxis.major_label_text_font = 'futura'
            p4.xaxis.axis_label_text_font_size = '11pt'
            p4.yaxis.axis_label_text_font_size = '11pt'
            p4.xaxis.major_label_text_font_size = '8pt'
            p4.yaxis.major_label_text_font_size = '8pt'
            p4.yaxis[0].formatter.precision = 1
            p4.title.align = 'center'
            p4.title.text_font_size = '16pt'
            p4.title.text_font = 'futura'

            min_x_date = astrotime(min_x_range, format='mjd').datetime
            max_x_date = astrotime(max_x_range, format='mjd').datetime

            p4.extra_x_ranges = {"gregorian date": Range1d(
                start=min_x_date, end=max_x_date)}
            p4.add_layout(
                DatetimeAxis(
                    major_label_text_font_size='8pt',
                    axis_label='Time (' + photoutime + '/Gregorian)',
                    major_label_text_font='futura',
                    axis_label_text_font='futura',
                    major_tick_in=0,
                    x_range_name="gregorian date",
                    axis_label_text_font_size='11pt'),
                'below')

            if mjdmax:
                min_xm_range = (min_x_range - mjdmax) * redshiftfactor
                max_xm_range = (max_x_range - mjdmax) * redshiftfactor
                p4.extra_x_ranges["time since max"] = Range1d(
                    start=min_xm_range, end=max_xm_range)
                p4.add_layout(
                    LinearAxis(
                        axis_label="Time since max (" + dayframe + ")",
                        major_label_text_font_size='8pt',
                        major_label_text_font='futura',
                        axis_label_text_font='futura',
                        x_range_name="time since max",
                        axis_label_text_font_size='11pt'),
                    'above')

            if distancemod:
                min_y_absmag = min([(x - y) * areacorr
                                    for x, y in list(zip(photofl, photoflerrs))])
                max_y_absmag = max([(x + y) * areacorr
                                    for x, y in list(zip(photofl, photoflerrs))])
                [min_y_absmag, max_y_absmag] = [
                    min_y_absmag - 0.1 * (max_y_absmag - min_y_absmag),
                    max_y_absmag + 0.1 * (max_y_absmag - min_y_absmag)
                ]
                p4.extra_y_ranges = {"abs mag": Range1d(
                    start=min_y_absmag, end=max_y_absmag)}
                p4.add_layout(
                    LinearAxis(
                        axis_label="Luminosity in band (ergs s⁻¹)",
                        major_label_text_font_size='8pt',
                        major_label_text_font='futura',
                        axis_label_text_font='futura',
                        y_range_name="abs mag",
                        axis_label_text_font_size='11pt'),
                    'right')
                p4.yaxis[1].formatter.precision = 1
            p4.add_tools(hover)

            xs = []
            ys = []
            err_xs = []
            err_ys = []

            for x, y, xlowerr, xupperr, yerr in list(
                    zip(phototime, photofl, phototimelowererrs, phototimeuppererrs,
                        photoflerrs)):
                xs.append(x)
                ys.append(y)
                err_xs.append((x - xlowerr, x + xupperr))
                err_ys.append((y - yerr, y + yerr))

            enerset = set(photoener)
            enerunit = photouener[0] if photouener else ''

            for ener in enerset:
                indb = [i for i, j in enumerate(photoener) if j == ener]
                indt = [i for i, j in enumerate(phototype) if not j]
                # Should always have upper error if have lower error.
                indnex = [i for i, j in enumerate(phototimelowererrs) if j == 0.]
                indyex = [i for i, j in enumerate(phototimelowererrs) if j > 0.]
                indney = [i for i, j in enumerate(photoflerrs) if j == 0.]
                indyey = [i for i, j in enumerate(photoflerrs) if j > 0.]
                indne = set(indb).intersection(indt).intersection(
                    indney).intersection(indnex)
                indye = set(indb).intersection(indt).intersection(
                    set(indyey).union(indyex))

                enerlabel = str(ener) + " " + enerunit

                noerrorlegend = enerlabel if len(indye) == 0 and len(
                    indne) > 0 else ''

                data = dict(
                    x=[phototime[i] for i in indne],
                    y=[photofl[i] for i in indne],
                    err=[photoflerrs[i] for i in indne],
                    desc=[photoener[i] for i in indne],
                    instr=[photoinstru[i] for i in indne],
                    src=[photosource[i] for i in indne])
                if distancemod:
                    data['yabs'] = [
                        str(round_sig(
                            photofl[i] * areacorr, sig=3)) for i in indne
                    ]
                    data['abserr'] = [
                        str(round_sig(
                            photoflerrs[i] * areacorr, sig=3)) for i in indne
                    ]
                if hastimeerrs:
                    data['xle'] = [phototimelowererrs[i] for i in indne]
                    data['xue'] = [phototimeuppererrs[i] for i in indne]

                source = ColumnDataSource(data)
                p4.circle(
                    'x',
                    'y',
                    source=source,
                    color=xraycolorf(ener),
                    fill_color="white",
                    legend=noerrorlegend,
                    size=4)

                yeserrorlegend = enerlabel if len(indye) > 0 else ''

                data = dict(
                    x=[phototime[i] for i in indye],
                    y=[photofl[i] for i in indye],
                    err=[photoflerrs[i] for i in indye],
                    desc=[photoener[i] for i in indye],
                    instr=[photoinstru[i] for i in indye],
                    src=[photosource[i] for i in indye])
                if distancemod:
                    data['yabs'] = [
                        str(round_sig(
                            photofl[i] * areacorr, sig=3)) for i in indye
                    ]
                    data['abserr'] = [
                        str(round_sig(
                            photoflerrs[i] * areacorr, sig=3)) for i in indye
                    ]
                if hastimeerrs:
                    data['xle'] = [phototimelowererrs[i] for i in indye]
                    data['xue'] = [phototimeuppererrs[i] for i in indye]

                source = ColumnDataSource(data)
                p4.multi_line(
                    [err_xs[x] for x in indye], [[ys[x], ys[x]] for x in indye],
                    color=xraycolorf(ener))
                p4.multi_line(
                    [[xs[x], xs[x]] for x in indye], [err_ys[x] for x in indye],
                    color=xraycolorf(ener))
                p4.circle(
                    'x',
                    'y',
                    source=source,
                    color=xraycolorf(ener),
                    legend=yeserrorlegend,
                    size=4)

                upplimlegend = enerlabel if len(indye) == 0 and len(
                    indne) == 0 else ''

                indt = [i for i, j in enumerate(phototype) if j]
                ind = set(indb).intersection(indt)
                data = dict(
                    x=[phototime[i] for i in ind],
                    y=[photofl[i] for i in ind],
                    err=[photoflerrs[i] for i in ind],
                    desc=[photoener[i] for i in ind],
                    instr=[photoinstru[i] for i in ind],
                    src=[photosource[i] for i in ind])
                if distancemod:
                    data['yabs'] = [
                        str(round_sig(
                            photofl[i] * areacorr, sig=3)) for i in ind
                    ]
                    data['abserr'] = [
                        str(round_sig(
                            photoflerrs[i] * areacorr, sig=3)) for i in ind
                    ]
                if hastimeerrs:
                    data['xle'] = [phototimelowererrs[i] for i in ind]
                    data['xue'] = [phototimeuppererrs[i] for i in ind]

                source = ColumnDataSource(data)
                # Currently Bokeh doesn't support tooltips for inverted_triangle,
                # so hide an invisible circle behind for the tooltip
                p4.circle('x', 'y', source=source, alpha=0.0, size=7)
                p4.inverted_triangle(
                    'x',
                    'y',
                    source=source,
                    color=xraycolorf(ener),
                    legend=upplimlegend,
                    size=7)

            p4.legend.label_text_font = 'futura'
            p4.legend.label_text_font_size = '8pt'
            p4.legend.label_width = 20
            p4.legend.label_height = 14
            p4.legend.glyph_height = 14

        hasimage = False
        skyhtml = ''
        if 'ra' in catalog[entry] and 'dec' in catalog[
                entry] and args.collecthosts:
            snra = catalog[entry]['ra'][0]['value']
            sndec = catalog[entry]['dec'][0]['value']
            try:
                c = coord(ra=snra, dec=sndec, unit=(un.hourangle, un.deg))
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                warnings.warn('Malformed angle for event ' + entry + '.')
            else:
                # if 'lumdist' in catalog[entry] and float(catalog[entry]['lumdist'][0]['value']) > 0.:
                #    if 'host' in catalog[entry] and catalog[entry]['host'][0]['value'] == 'Milky Way':
                #        sdssimagescale = max(0.05,0.4125/float(catalog[entry]['lumdist'][0]['value']))
                #    else:
                #    sdssimagescale = max(0.5,20.6265/float(catalog[entry]['lumdist'][0]['value']))
                # else:
                #    if 'host' in catalog[entry] and catalog[entry]['host'][0]['value'] == 'Milky Way':
                #        sdssimagescale = 0.006
                #    else:
                #    sdssimagescale = 0.5
                # dssimagescale = 0.13889*sdssimagescale
                # At the moment, no way to check if host is in SDSS footprint
                # without comparing to empty image, which is only possible at fixed
                # angular resolution.
                sdssimagescale = 0.3
                dssimagescale = 0.13889 * sdssimagescale

                imgsrc = ''
                hasimage = True
                if event_name in host_img_dict:
                    imgsrc = host_img_dict[event_name]
                else:
                    try:
                        response = urllib.request.urlopen(
                            'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?ra='
                            + str(c.ra.deg) + '&dec=' + str(c.dec.deg) + '&scale='
                            + str(sdssimagescale) + '&width=500&height=500&opt=G',
                            timeout=60)
                        resptxt = response.read()
                    except (KeyboardInterrupt, SystemExit):
                        raise
                    except:
                        hasimage = False
                    else:
                        with open(paths.output_html + event_file_name + '-host.jpg',
                                  'wb') as f:
                            f.write(resptxt)
                        imgsrc = 'SDSS'

                    if hasimage and filecmp.cmp(
                            paths.output_html + event_file_name + '-host.jpg',
                            'astrocats/supernovae/input/missing.jpg'):
                        hasimage = False

                    if not hasimage:
                        hasimage = True
                        url = (
                            "http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Position="
                            + str(urllib.parse.quote_plus(snra + " " + sndec)) +
                            "&coordinates=J2000&coordinates=&projection=Tan&pixels=500&size="
                            + str(dssimagescale) +
                            "&float=on&scaling=Log&resolver=SIMBAD-NED" +
                            "&Sampler=_skip_&Deedger=_skip_&rotation=&Smooth=&lut=colortables%2Fb-w-linear.bin&PlotColor=&grid=_skip_&gridlabels=1"
                            +
                            "&catalogurl=&CatalogIDs=on&RGB=1&survey=DSS2+IR&survey=DSS2+Red&survey=DSS2+Blue&IOSmooth=&contour=&contourSmooth=&ebins=null"
                        )

                        try:
                            response = urllib.request.urlopen(url, timeout=60)
                            bandsoup = BeautifulSoup(response, "html5lib")
                        except (KeyboardInterrupt, SystemExit):
                            raise
                        except:
                            hasimage = False
                        else:
                            images = bandsoup.findAll('img')
                            imgname = ''
                            for image in images:
                                if "Quicklook RGB image" in image.get('alt', ''):
                                    imgname = image.get('src', '').split('/')[-1]

                            if imgname:
                                try:
                                    response = urllib.request.urlopen(
                                        'http://skyview.gsfc.nasa.gov/tempspace/fits/'
                                        + imgname)
                                except (KeyboardInterrupt, SystemExit):
                                    raise
                                except:
                                    hasimage = False
                                else:
                                    with open(paths.output_html + event_file_name +
                                              '-host.jpg', 'wb') as f:
                                        f.write(response.read())
                                    imgsrc = 'DSS'
                            else:
                                hasimage = False

            if hasimage:
                if imgsrc == 'SDSS':
                    host_img_dict[event_name] = 'SDSS'
                    skyhtml = (
                        '<a href="http://skyserver.sdss.org/DR12/en/tools/chart/navi.aspx?opt=G&ra='
                        + str(c.ra.deg) + '&dec=' + str(c.dec.deg) +
                        '&scale=0.15"><img src="' + event_file_name +
                        '-host.jpg" width=250></a>')
                elif imgsrc == 'DSS':
                    host_img_dict[event_name] = 'DSS'
                    url = (
                        "http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Position="
                        + str(urllib.parse.quote_plus(snra + " " + sndec)) +
                        "&coordinates=J2000&coordinates=&projection=Tan&pixels=500&size="
                        + str(dssimagescale) +
                        "float=on&scaling=Log&resolver=SIMBAD-NED" +
                        "&Sampler=_skip_&Deedger=_skip_&rotation=&Smooth=&lut=colortables%2Fb-w-linear.bin&PlotColor=&grid=_skip_&gridlabels=1"
                        +
                        "&catalogurl=&CatalogIDs=on&RGB=1&survey=DSS2+IR&survey=DSS2+Red&survey=DSS2+Blue&IOSmooth=&contour=&contourSmooth=&ebins=null"
                    )
                    skyhtml = ('<a href="' + url + '"><img src="' + event_file_name +
                               '-host.jpg" width=250></a>')
            else:
                host_img_dict[event_name] = 'None'

        if dohtml and args.writehtml:
            # if (photoavail and spectraavail) and dohtml and args.writehtml:
            plots = []
            if photoavail:
                if photochecks:
                    p1box = column(p1, photochecks)
                else:
                    p1box = p1
                plots += [p1box]
            if spectraavail:
                plots += [column(p2, bokehrow(binslider, spacingslider))]
            if radioavail:
                plots += [p3]
            if xrayavail:
                plots += [p4]

            p = layout(
                [plots[i:i + 2] for i in range(0, len(plots), 2)],
                ncols=2,
                toolbar_location=None)

            html = '<html><head><title>' + event_name + '</title>'
            if photoavail or spectraavail or radioavail or xrayavail:
                html = file_html(p, CDN, event_name)
                # html = html + '''<link href="https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.css" rel="stylesheet" type="text/css">
                #    <script src="https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.js"></script>''' + script + '</head><body>'
            else:
                html = '<html><title></title><body></body></html>'

            # if photoavail and spectraavail:
            #    html = html + div['p1'] + div['p2']# + div['binslider'] + div['spacingslider']
            # elif photoavail:
            #    html = html + div['p1']
            # elif spectraavail:
            #    html = html + div['p2'] + div['binslider'] + div['spacingslider']

            # html = html + '</body></html>'

            html = html.replace(
                '<body>',
                '''<body class='event-body'><div style="padding-bottom:8px;"><strong>Disclaimer:</strong> All data collected by the OSC was originally generated by others, if you intend to use this data in a publication, we ask that you please cite the linked sources and/or contact the sources of the data directly. Data sources are revealed by hovering over the data with your cursor.</div>''')
            html = re.sub(r'(\<\/title\>)', r'''\1\n
                <base target="_parent" />\n
                <link rel="stylesheet" href="https://sne.space/astrocats/astrocats/supernovae/html/event.css" type="text/css">\n
                <script type="text/javascript" src="https://sne.space/astrocats/astrocats/supernovae/scripts/marks.js" type="text/css"></script>\n
                <script type="text/javascript">\n
                    if(top==self)\n
                    this.location="''' + event_name + '''"\n
                </script>''', html)

            html = re.sub(
                r'(\<\/body\>)', '<div class="event-download">' + r'<a href="' +
                SNE_LINK_DIR + event_file_name + r'.json" download>' +
                r'Download all data for ' + event_name + r'</a></div>\n\1', html)
            issueargs = '?title=' + ('[' + event_name + '] <Descriptive issue title>').encode('ascii', 'xmlcharrefreplace').decode("utf-8") + '&body=' + \
                ('Please describe the issue with ' + event_name + '\'s data here, be as descriptive as possible! ' +
                 'If you believe the issue appears in other events as well, please identify which other events the issue possibly extends to.').encode('ascii', 'xmlcharrefreplace').decode("utf-8")
            html = re.sub(
                r'(\<\/body\>)', '<div class="event-issue">' +
                r'<a href="https://github.com/astrocatalogs/supernovae/issues/new'
                + issueargs + r'" target="_blank">' + r'Report an issue with ' +
                event_name + r'</a></div>\n\1', html)

            newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Event metadata</h3><table class="event-table"><tr><th width=100px class="event-cell">Quantity</th><th class="event-cell">Value<sup>Sources</sup> [Kind]</th></tr>\n'
            for key in COLUMN_KEYS:
                if key in catalog[entry] and key not in EVENT_IGNORE_KEY and len(
                        catalog[entry][key]) > 0:
                    keyhtml = ''
                    if isinstance(catalog[entry][key], str):
                        if key in ['photolink', 'spectralink', 'radiolink',
                                   'xraylink']:
                            keysplit = catalog[entry][key].split(',')
                            if keysplit:
                                num = int(keysplit[0])
                                keyhtml = keyhtml + keysplit[0] + ' ' + (
                                    infl.plural('spectrum', num)
                                    if key == 'spectralink' else infl.plural(
                                        'detection', num))
                                if len(keysplit) == 3:
                                    keyhtml = keyhtml + \
                                        '<br>[' + keysplit[1] + ' – ' + \
                                        keysplit[2] + ' days from max]'
                        else:
                            subentry = re.sub('<[^<]+?>', '', catalog[entry][key])
                            keyhtml = keyhtml + subentry
                    else:
                        for r, row in enumerate(catalog[entry][key]):
                            if 'value' in row and 'source' in row:
                                sources = [
                                    str(x)
                                    for x in sorted(
                                        [x.strip()
                                         for x in row['source'].split(',')],
                                        key=lambda x: float(x) if is_number(x) else float("inf"))
                                ]
                                sourcehtml = ''
                                for s, source in enumerate(sources):
                                    sourcehtml = sourcehtml + \
                                        (', ' if s > 0 else '') + r'<a href="#source' + \
                                        source + r'">' + source + r'</a>'
                                keyhtml = keyhtml + (r'<br>' if r > 0 else '')
                                keyhtml = keyhtml + "<div class='singletooltip'>"
                                if 'derived' in row and row['derived']:
                                    keyhtml = keyhtml + '<span class="derived">'
                                keyhtml = keyhtml + row['value']
                                if ((key == 'maxdate' or key == 'maxabsmag' or
                                     key == 'maxappmag') and
                                        'maxband' in catalog[entry] and
                                        catalog[entry]['maxband']):
                                    keyhtml = keyhtml + \
                                        r' [' + catalog[entry]['maxband'][0]['value'] + ']'
                                if 'e_value' in row:
                                    keyhtml = keyhtml + r' ± ' + row['e_value']
                                if 'derived' in row and row['derived']:
                                    keyhtml = keyhtml + '</span>'

                                # Mark erroneous button
                                sourceids = []
                                idtypes = []
                                for alias in row['source'].split(','):
                                    for source in catalog[entry]['sources']:
                                        if source['alias'] == alias:
                                            if 'bibcode' in source:
                                                sourceids.append(source['bibcode'])
                                                idtypes.append('bibcode')
                                            else:
                                                sourceids.append(source['name'])
                                                idtypes.append('name')
                                if not sourceids or not idtypes:
                                    raise ValueError(
                                        'Unable to find associated source by alias!')
                                edit = "true" if os.path.isfile(
                                    'astrocats/supernovae/input/sne-internal/' +
                                    get_event_filename(
                                        entry) + '.json') else "false"
                                keyhtml = (
                                    keyhtml +
                                    "<span class='singletooltiptext'><button class='singlemarkerror' type='button' onclick='markError(\""
                                    + entry + "\", \"" + key + "\", \"" +
                                    ','.join(idtypes) + "\", \"" +
                                    ','.join(sourceids) + "\", \"" + edit +
                                    "\")'>Flag as erroneous</button></span>")
                                keyhtml = keyhtml + r'</div><sup>' + sourcehtml + r'</sup>'
                            elif isinstance(row, str):
                                keyhtml = keyhtml + \
                                    (r'<br>' if r > 0 else '') + row.strip()

                    if keyhtml:
                        newhtml = (
                            newhtml + r'<tr><td class="event-cell">' +
                            EVENT_PAGE_HEADER[key] +
                            r'</td><td width=250px class="event-cell">' + keyhtml)

                    newhtml = newhtml + r'</td></tr>\n'

            newhtml = newhtml + r'</table><em>Values that are colored <span class="derived">purple</span> were computed by the OSC using values provided by the specified sources.</em></div>\n\1'
            html = re.sub(r'(\<\/body\>)', newhtml, html)

            if 'sources' in catalog[entry] and len(catalog[entry]['sources']):
                newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Sources of data</h3><table class="event-table"><tr><th width=30px class="event-cell">ID</th><th class="event-cell">Source Info</th></tr><tr><th colspan="2" class="event-cell">Primary Sources</th></tr>\n'
                first_secondary = False
                for source in catalog[entry]['sources']:
                    biburl = ''
                    if 'bibcode' in source:
                        biburl = 'http://adsabs.harvard.edu/abs/' + \
                            source['bibcode']

                    refurl = ''
                    if 'url' in source:
                        refurl = source['url']

                    sourcename = source['name'] if 'name' in source else source[
                        'bibcode']
                    if not first_secondary and source.get('secondary', False):
                        first_secondary = True
                        newhtml += r'<tr><th colspan="2" class="event-cell">Secondary Sources</th></tr>\n'
                    newhtml = (
                        newhtml + r'<tr><td class="event-cell" id="source' +
                        source['alias'] + '">' + source['alias'] +
                        r'</td><td width=250px class="event-cell">' + (
                            (((r'<a href="' + refurl + '">')
                              if refurl else '') + sourcename.encode(
                                  'ascii', 'xmlcharrefreplace').decode("utf-8") + (
                                      (r'</a>\n') if refurl else '') + r'<br>')
                            if 'bibcode' not in source or sourcename !=
                            source['bibcode'] else '') + (
                                (source['reference'] + r'<br>')
                                if 'reference' in source else '') +
                        ((r'\n[' + (
                            ('<a href="' + biburl + '">')
                            if 'reference' in source else '') + source['bibcode'] +
                          (r'</a>' if 'reference' in source else '') + ']')
                         if 'bibcode' in source else '') + r'</td></tr>\n')
                newhtml = newhtml + r'</table><em>Sources are presented in order of importation, not in order of importance.</em></div>\n'

                if hasimage:
                    newhtml = newhtml + \
                        '<div class="event-host-div"><h3 class="event-host-title">Host Image</h3>' + skyhtml
                    newhtml = newhtml + \
                        r'</table><em>Host images are taken from SDSS if available; if not, DSS is used.</em></div>\n'

            newhtml = newhtml + r'\n\1'

            html = re.sub(r'(\<\/body\>)', newhtml, html)

            html_filename = paths.output_html + event_file_name + ".html"
            with gzip.open(html_filename + ".gz", 'wt') as fff:
                touch(html_filename)
                fff.write(html)

        # Necessary to clear Bokeh state
        reset_output()

        # if spectraavail and dohtml:
        #    sys.exit()

        # if fcnt > 100:
        #    sys.exit()

        # Save this stuff because next line will delete it.
        if args.writecatalog:
            # Construct array for Bishop's webpage
            # Things David wants in this file: names (aliases), max mag, max mag
            # date (gregorian), type, redshift (helio), redshift (host), r.a.,
            # dec., # obs., link
            SNE_PAGES.append(
                [entry, ",".join([x['value'] for x in catalog[entry]['alias']]),
                 get_first_value(catalog, entry, 'maxappmag'),
                 get_first_value(catalog, entry, 'maxdate'),
                 get_first_value(catalog, entry, 'claimedtype'),
                 get_first_value(catalog, entry, 'redshift'),
                 get_first_kind(catalog, entry, 'redshift'),
                 get_first_value(catalog, entry, 'ra'),
                 get_first_value(catalog, entry, 'dec'),
                 catalog[entry]['numphoto'], 'https://sne.space/' + plotlink])

            if 'sources' in catalog[entry]:
                lsourcedict = {}
                for sourcerow in catalog[entry]['sources']:
                    if 'name' not in sourcerow:
                        continue
                    strippedname = re.sub('<[^<]+?>', '', sourcerow['name'].encode(
                        'ascii', 'xmlcharrefreplace').decode("utf-8"))
                    alias = sourcerow['alias']
                    if 'bibcode' in sourcerow and 'secondary' not in sourcerow:
                        lsourcedict[alias] = {
                            'bibcode': sourcerow['bibcode'],
                            'count': 0
                        }
                    if strippedname in sourcedict:
                        sourcedict[strippedname] += 1
                    else:
                        sourcedict[strippedname] = 1

                for key in catalog[entry].keys():
                    if isinstance(catalog[entry][key], list):
                        for row in catalog[entry][key]:
                            if 'source' in row:
                                for lsource in lsourcedict:
                                    if lsource in row['source'].split(','):
                                        if key == 'spectra':
                                            lsourcedict[lsource]['count'] += 10
                                        else:
                                            lsourcedict[lsource]['count'] += 1

                ssources = sorted(
                    list(lsourcedict.values()),
                    key=lambda x: x['count'],
                    reverse=True)
                if ssources:
                    catalog[entry]['references'] = ','.join(
                        [y['bibcode'] for y in ssources[:5]])

            lcspye.append(catalog[entry]['numphoto'] >= 5 and
                          catalog[entry]['numspectra'] > 0)
            lconly.append(catalog[entry]['numphoto'] >= 5 and
                          catalog[entry]['numspectra'] == 0)
            sponly.append(catalog[entry]['numphoto'] < 5 and
                          catalog[entry]['numspectra'] > 0)
            lcspno.append(catalog[entry]['numphoto'] < 5 and
                          catalog[entry]['numspectra'] == 0)

            hasalc.append(catalog[entry]['numphoto'] >= 5)
            hasasp.append(catalog[entry]['numspectra'] > 0)

            totalphoto += catalog[entry]['numphoto']
            totalspectra += catalog[entry]['numspectra']

            # Delete unneeded data from catalog, add blank entries when data
            # missing.
            catalogcopy[entry] = OrderedDict()
            for col in COLUMN_KEYS:
                if col in catalog[entry]:
                    catalogcopy[entry][col] = deepcopy(catalog[entry][col])
                else:
                    catalogcopy[entry][col] = None

        del catalog[entry]

        if args.test and spectraavail and photoavail:
            break

    # Write it all out at the end
    if args.writecatalog and not args.eventlist:
        catalog = deepcopy(catalogcopy)

        # Write the MD5 checksums
        jsonstring = json.dumps(md5_dict, indent='\t', separators=(',', ':'))
        with open(paths.md5_file + testsuffix, 'w') as f:
            f.write(jsonstring)

        # Write the host image info
        if args.collecthosts:
            jsonstring = json.dumps(
                host_img_dict, indent='\t', separators=(',', ':'))
            with open(paths.host_imgs_file + testsuffix, 'w') as f:
                f.write(jsonstring)

        if not args.boneyard:
            # Things David wants in this file: names (aliases), max mag, max mag
            # date (gregorian), type, redshift, r.a., dec., # obs., link
            with open(paths.output_html + 'SNE_PAGES.csv' + testsuffix, 'w') as f:
                csvout = csv.writer(f, quotechar='"', quoting=csv.QUOTE_ALL)
                for row in SNE_PAGES:
                    csvout.writerow(row)

            # Make a few small files for generating charts
            with open(paths.output_html + 'sources.csv' + testsuffix, 'w') as f:
                sortedsources = sorted(
                    list(sourcedict.items()),
                    key=operator.itemgetter(1),
                    reverse=True)
                csvout = csv.writer(f)
                csvout.writerow(['Source', 'Number'])
                for source in sortedsources:
                    csvout.writerow(source)

            with open(paths.output_html + 'pie.csv' + testsuffix, 'w') as f:
                csvout = csv.writer(f)
                csvout.writerow(['Category', 'Number'])
                csvout.writerow(['Has light curve and spectra', sum(lcspye)])
                csvout.writerow(['Has light curve only', sum(lconly)])
                csvout.writerow(['Has spectra only', sum(sponly)])
                csvout.writerow(['No light curve or spectra', sum(lcspno)])

            with open(
                    paths.output_html + 'info-snippets/hasphoto.html' + testsuffix,
                    'w') as f:
                f.write("{:,}".format(sum(hasalc)))
            with open(paths.output_html + 'info-snippets/hasspectra.html' +
                      testsuffix, 'w') as f:
                f.write("{:,}".format(sum(hasasp)))
            with open(
                    paths.output_html + 'info-snippets/snecount.html' + testsuffix,
                    'w') as f:
                f.write("{:,}".format(len(catalog)))
            with open(paths.output_html + 'info-snippets/photocount.html' +
                      testsuffix, 'w') as f:
                f.write("{:,}".format(totalphoto))
            with open(paths.output_html + 'info-snippets/spectracount.html' +
                      testsuffix, 'w') as f:
                f.write("{:,}".format(totalspectra))

            ctypedict = dict()
            for entry in catalog:
                cleanedtype = ''
                if 'claimedtype' in catalog[entry] and catalog[entry][
                        'claimedtype']:
                    maxsources = 0
                    for ct in catalog[entry]['claimedtype']:
                        sourcecount = len(ct['source'].split(','))
                        if sourcecount > maxsources:
                            maxsources = sourcecount
                            cleanedtype = ct['value'].strip('?* ')
                if not cleanedtype:
                    cleanedtype = 'Unknown'
                if cleanedtype in ctypedict:
                    ctypedict[cleanedtype] += 1
                else:
                    ctypedict[cleanedtype] = 1
            sortedctypes = sorted(
                list(ctypedict.items()), key=operator.itemgetter(1), reverse=True)
            with open(paths.output_html + 'types.csv' + testsuffix, 'w') as f:
                csvout = csv.writer(f)
                csvout.writerow(['Type', 'Number'])
                for ctype in sortedctypes:
                    csvout.writerow(ctype)

            with open(paths.output_html + 'sitemap.xml', 'w') as f:
                sitemapxml = SITEMAP_TEMPLATE
                sitemaplocs = ''
                for key in catalog.keys():
                    sitemaplocs = sitemaplocs + "  <url>\n    <loc>https://sne.space/sne/" + \
                        key + "</loc>\n  </url>\n"
                sitemapxml = sitemapxml.replace('{0}', sitemaplocs)
                f.write(sitemapxml)

            # Ping Google to let them know sitemap has been updated
            response = urllib.request.urlopen(GOOGLE_PING_URL)

        # Prune extraneous fields not required for main catalog file
        catalogcopy = OrderedDict()
        for entry in catalog:
            catalogcopy[entry] = OrderedDict()
            for col in catalog[entry]:
                catalogcopy[entry][col] = deepcopy(catalog[entry][col])
                if catalogcopy[entry][col]:
                    for row in catalogcopy[entry][col]:
                        if 'source' in row:
                            del row['source']
                        if 'u_value' in row:
                            del row['u_value']
        catalog = deepcopy(catalogcopy)

        # Convert to array since that's what datatables expects
        catalog = list(catalog.values())

        if args.boneyard:
            catprefix = 'bones'
        else:
            catprefix = 'catalog'

        jsonstring = json.dumps(catalog, separators=(',', ':'))
        with open(paths.output + catprefix + '.min.json' + testsuffix, 'w') as f:
            f.write(jsonstring)

        jsonstring = json.dumps(catalog, indent='\t', separators=(',', ':'))
        with open(paths.output + catprefix + '.json' + testsuffix, 'w') as f:
            f.write(jsonstring)

        with open(paths.output_html + 'table-templates/' + catprefix + '.html' +
                  testsuffix, 'w') as f:
            f.write(
                '<table id="example" class="display" cellspacing="0" width="100%">\n')
            f.write('\t<thead>\n')
            f.write('\t\t<tr>\n')
            for h in HEADER:
                f.write('\t\t\t<th class="' + h + '" title="' + DEF_TITLES[h] + '">' +
                        HEADER[h] + '</th>\n')
            f.write('\t\t</tr>\n')
            f.write('\t</thead>\n')
            f.write('\t<tfoot>\n')
            f.write('\t\t<tr>\n')
            for h in HEADER:
                f.write('\t\t\t<th class="' + h + '" title="' + DEF_TITLES[h] + '">' +
                        HEADER[h] + '</th>\n')
            f.write('\t\t</tr>\n')
            f.write('\t</tfoot>\n')
            f.write('</table>\n')

        with open(paths.output + catprefix + '.min.json', 'rb') as f_in, gzip.open(
                paths.output + catprefix + '.min.json.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        if not args.boneyard:
            names = OrderedDict()
            for ev in catalog:
                names[ev['name']] = [x['value'] for x in ev['alias']]
            jsonstring = json.dumps(names, separators=(',', ':'))
            with open(paths.output + 'names.min.json' + testsuffix, 'w') as f:
                f.write(jsonstring)

        if args.deleteorphans and not args.boneyard:

            safefiles = [os.path.basename(x) for x in files]
            safefiles += SAFE_FILES
            for myfile in glob(paths.output_json + '*.json'):
                if not os.path.basename(myfile) in safefiles:
                    print('Deleting orphan ' + myfile)
                    # os.remove(myfile)
    return


def load_dict_file(fname, log):
    log_str = "File '{}'".format(fname)
    if os.path.isfile(fname):
        fdict = json.load(open(fname, 'r'))
        log_str += " loaded."
        log.debug(log_str)
    else:
        fdict = {}
        log_str += " does not exist.  Initializing empty dictionary."
        log.info(log_str)
    return fdict
