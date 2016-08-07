"""Spectra Plotting Methods associatd with WebCat Script.
"""

from copy import deepcopy
from math import ceil, isnan
import operator
import numpy

from bokeh.models import HoverTool, Range1d, LinearAxis, ColumnDataSource, CustomJS, Slider
from bokeh.plotting import Figure

from astrocats.catalog.utils import is_number
from .constants import TOOLS_LIST, DEF_COLORS
from .utils import label_format


def plot_spectra(catalog, entry, mjdmax, redshiftfactor):
    # FIX:
    event_name = entry

    spectrumwave = []
    spectrumflux = []
    spectrumerrs = []
    spectrummjdmax = []
    hasepoch = True
    if 'redshift' in catalog[entry]:
        z = float(catalog[entry]['redshift'][0]['value'])
    catalog[entry]['spectra'] = list(
        filter(None, [x if 'data' in x else None
                      for x in catalog[entry]['spectra']]))
    for spectrum in catalog[entry]['spectra']:
        spectrumdata = deepcopy(spectrum['data'])
        specslice = ceil(float(len(spectrumdata)) / 10000)
        spectrumdata = spectrumdata[::specslice]
        spectrumdata = [x for x in spectrumdata
                        if is_number(x[1]) and not isnan(float(x[1]))]
        specrange = range(len(spectrumdata))

        if 'deredshifted' in spectrum and spectrum[
                'deredshifted'] and 'redshift' in catalog[entry]:
            spectrumwave.append(
                [float(spectrumdata[x][0]) * (1.0 + z) for x in specrange])
        else:
            spectrumwave.append([float(spectrumdata[x][0])
                                 for x in specrange])

        spectrumflux.append([float(spectrumdata[x][1]) for x in specrange])
        if 'errorunit' in spectrum:
            spectrumerrs.append([float(spectrumdata[x][2])
                                 for x in specrange])
            spectrumerrs[-1] = [x if is_number(x) and not isnan(float(x))
                                else 0. for x in spectrumerrs[-1]]

        if 'u_time' not in spectrum or 'time' not in spectrum:
            hasepoch = False

        if 'u_time' in spectrum and 'time' in spectrum and spectrum[
                'u_time'] == 'MJD' and 'redshift' in catalog[
                    entry] and mjdmax:
            specmjd = (float(spectrum['time']) - mjdmax) * redshiftfactor
            spectrummjdmax.append(specmjd)

    nspec = len(catalog[entry]['spectra'])

    prunedwave = []
    prunedflux = []
    for i in reversed(range(nspec)):
        ri = nspec - i - 1
        prunedwave.append([])
        prunedflux.append([])
        for wi, wave in enumerate(spectrumwave[i]):
            exclude = False
            if 'exclude' in catalog[entry]['spectra'][i]:
                for exclusion in catalog[entry]['spectra'][i]['exclude']:
                    if 'below' in exclusion:
                        if wave <= float(exclusion['below']):
                            exclude = True
                    elif 'above' in exclusion:
                        if wave >= float(exclusion['above']):
                            exclude = True
            if not exclude:
                prunedwave[ri].append(wave)
                prunedflux[ri].append(spectrumflux[i][wi])

    prunedwave = list(reversed(prunedwave))
    prunedflux = list(reversed(prunedflux))

    prunedscaled = deepcopy(prunedflux)
    for f, flux in enumerate(prunedscaled):
        std = numpy.std(flux)
        prunedscaled[f] = [x / std for x in flux]

    y_height = 0.
    y_offsets = [0. for x in range(nspec)]
    for i in reversed(range(nspec)):
        y_offsets[i] = y_height
        if (i - 1 >= 0 and 'time' in catalog[entry]['spectra'][i] and
                'time' in catalog[entry]['spectra'][i - 1] and
                catalog[entry]['spectra'][i]['time'] ==
                catalog[entry]['spectra'][i - 1]['time']):
            ydiff = 0
        else:
            ydiff = 0.8 * (max(prunedscaled[i]) - min(prunedscaled[i]))
        prunedscaled[i] = [j + y_height for j in prunedscaled[i]]
        y_height += ydiff

    maxsw = max(list(map(max, prunedwave)))
    minsw = min(list(map(min, prunedwave)))
    maxfl = max(list(map(max, prunedscaled)))
    minfl = min(list(map(min, prunedscaled)))
    maxfldiff = max(
        map(operator.sub, list(map(max, prunedscaled)), list(
            map(min, prunedscaled))))
    x_buffer = 0.0  # 0.1*(maxsw - minsw)
    x_range = [-x_buffer + minsw, x_buffer + maxsw]
    y_buffer = 0.1 * maxfldiff
    y_range = [-y_buffer + minfl, y_buffer + maxfl]

    for f, flux in enumerate(prunedscaled):
        prunedscaled[f] = [x - y_offsets[f] for x in flux]

    tt2 = [("Source ID(s)", "@src")]
    if 'redshift' in catalog[entry]:
        tt2 += [("λ (rest)", "@xrest{1.1} Å")]
    tt2 += [
        ("λ (obs)", "@x{1.1} Å"), ("Flux", "@yorig"),
        ("Flux unit", "@fluxunit")
    ]

    if hasepoch:
        tt2 += [("Epoch (" + spectrum['u_time'] + ")", "@epoch{1.11}")]

    if mjdmax:
        tt2 += [("Rest days to max", "@mjdmax{1.11}")]

    hover = HoverTool(tooltips=tt2)

    p2 = Figure(title='Spectra for ' + event_name, x_axis_label=label_format('Observed Wavelength (Å)'), active_drag='box_zoom',
                y_axis_label=label_format('Flux (scaled)' + (' + offset'
                                                             if (nspec > 1) else '')), x_range=x_range, tools=TOOLS_LIST,  # sizing_mode = "scale_width",
                plot_width=485, plot_height=485, y_range=y_range, toolbar_location='above', toolbar_sticky=False)
    p2.xaxis.axis_label_text_font = 'futura'
    p2.yaxis.axis_label_text_font = 'futura'
    p2.xaxis.major_label_text_font = 'futura'
    p2.yaxis.major_label_text_font = 'futura'
    p2.xaxis.axis_label_text_font_size = '11pt'
    p2.yaxis.axis_label_text_font_size = '11pt'
    p2.xaxis.major_label_text_font_size = '8pt'
    p2.yaxis.major_label_text_font_size = '8pt'
    p2.title.align = 'center'
    p2.title.text_font_size = '16pt'
    p2.title.text_font = 'futura'
    p2.add_tools(hover)

    sources = []
    for i in range(len(prunedwave)):
        sl = len(prunedscaled[i])
        fluxunit = catalog[entry]['spectra'][i][
            'u_fluxes'] if 'u_fluxes' in catalog[entry]['spectra'][
                i] else ''

        data = dict(
            x0=prunedwave[i],
            y0=prunedscaled[i],
            yorig=spectrumflux[i],
            fluxunit=[label_format(fluxunit)] * sl,
            x=prunedwave[i],
            y=[y_offsets[i] + j for j in prunedscaled[i]],
            src=[catalog[entry]['spectra'][i]['source']] * sl)
        if 'redshift' in catalog[entry]:
            data['xrest'] = [x / (1.0 + z) for x in prunedwave[i]]
        if hasepoch:
            data['epoch'] = [catalog[entry]['spectra'][i]['time']
                             for j in prunedscaled[i]]
            if mjdmax and spectrummjdmax:
                data['mjdmax'] = [spectrummjdmax[i]
                                  for j in prunedscaled[i]]
        sources.append(ColumnDataSource(data))
        p2.line(
            'x',
            'y',
            source=sources[i],
            color=DEF_COLORS[i % len(DEF_COLORS)],
            line_width=2,
            line_join='round')

    if 'redshift' in catalog[entry]:
        minredw = minsw / (1.0 + z)
        maxredw = maxsw / (1.0 + z)
        p2.extra_x_ranges = {"other wavelength": Range1d(
            start=minredw, end=maxredw)}
        p2.add_layout(
            LinearAxis(
                axis_label="Restframe Wavelength (Å)",
                x_range_name="other wavelength",
                axis_label_text_font_size='11pt',
                axis_label_text_font='futura',
                major_label_text_font_size='8pt',
                major_label_text_font='futura'),
            'above')

    sdicts = dict(
        zip(['s' + str(x) for x in range(len(sources))], sources))
    callback = CustomJS(
        args=sdicts,
        code="""
        var yoffs = [""" + ','.join([str(x) for x in y_offsets]) + """];
        for (s = 0; s < """ + str(len(sources)) + """; s++) {
            var data = eval('s'+s).get('data');
            var redshift = """ +
        str(z if 'redshift' in catalog[entry] else 0.) + """;
            if (!('binsize' in data)) {
                data['binsize'] = 1.0
            }
            if (!('spacing' in data)) {
                data['spacing'] = 1.0
            }
            if (cb_obj.get('title') == 'Spacing') {
                data['spacing'] = cb_obj.get('value');
            } else {
                data['binsize'] = cb_obj.get('value');
            }
            var f = data['binsize']
            var space = data['spacing']
            var x0 = data['x0'];
            var y0 = data['y0'];
            var dx0 = x0[1] - x0[0];
            var yoff = space*yoffs[s];
            data['x'] = [x0[0] - 0.5*Math.max(0., f - dx0)];
            data['xrest'] = [(x0[0] - 0.5*Math.max(0., f - dx0))/(1.0 + redshift)];
            data['y'] = [y0[0] + yoff];
            var xaccum = 0.;
            var yaccum = 0.;
            for (i = 0; i < x0.length; i++) {
                var dx;
                if (i == 0) {
                    dx = x0[i+1] - x0[i];
                } else {
                    dx = x0[i] - x0[i-1];
                }
                xaccum += dx;
                yaccum += y0[i]*dx;
                if (xaccum >= f) {
                    data['x'].push(data['x'][data['x'].length-1] + xaccum);
                    data['xrest'].push(data['x'][data['x'].length-1]/(1.0 + redshift));
                    data['y'].push(yaccum/xaccum + yoff);
                    xaccum = 0.;
                    yaccum = 0.;
                }
            }
            eval('s'+s).trigger('change');
        }
    """)

    binslider = Slider(
        start=0,
        end=20,
        value=1,
        step=0.5,
        width=230,
        title=label_format("Bin size (Angstrom)"),
        callback=callback)
    spacingslider = Slider(
        start=0,
        end=2,
        value=1,
        step=0.02,
        width=230,
        title=label_format("Spacing"),
        callback=callback)

    return p2, binslider, spacingslider
