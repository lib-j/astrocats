"""X-ray Plotting Methods associatd with WebCat Script.
"""
from math import pi
from statistics import mean

from bokeh.models import HoverTool, DatetimeAxis, Range1d, LinearAxis, ColumnDataSource
from bokeh.plotting import Figure
from astropy.time import Time as astrotime
from astropy import units as un

from astrocats.catalog.utils import round_sig, get_sig_digits, xraycolorf

from .constants import TOOLS_LIST, RADIO_SIGMA


def plot_xray(catalog, entry, p1, dayframe, distancemod, mjdmax,
              min_x_range, max_x_range, photoavail, redshiftfactor):

    # FIX:
    event_name = entry

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
        tt += [(yaxis + " (" + photoufl[0].replace(
                "ergs/s/cm^2", "ergs s⁻¹ cm⁻²") + ")",
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

    return p4
