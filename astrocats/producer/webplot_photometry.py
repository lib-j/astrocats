"""Photometry Plotting Methods associatd with WebCat Script.
"""
from bokeh.models import HoverTool, DatetimeAxis, Range1d, LinearAxis, ColumnDataSource, CustomJS
from bokeh.plotting import Figure
from bokeh.models.widgets import Select
from astropy.time import Time as astrotime


from astrocats.catalog.utils import (bandaliasf, bandcolorf,
                                     bandgroupf, bandshortaliasf, bandwavef)

from .constants import TOOLS_LIST


def plot_photo(catalog, entry, dayframe, distancemod, redshiftfactor, mjdmax, min_x_range, max_x_range):
    """Plot Photometry Bokeh Figure
    """

    # FIX: is this needed for anything?
    event_name = entry

    phototime = [float(x['time']) for x in catalog[entry]['photometry']
                 if 'magnitude' in x]
    phototimelowererrs = [float(x['e_lower_time'])
                          if ('e_lower_time' in x and 'e_upper_time' in x)
                          else (float(x['e_time']) if 'e_time' in x else 0.)
                          for x in catalog[entry]['photometry'] if 'magnitude' in x]
    phototimeuppererrs = [float(x['e_upper_time'])
                          if ('e_lower_time' in x and 'e_upper_time' in x)
                          else (float(x['e_time']) if 'e_time' in x else 0.)
                          for x in catalog[entry]['photometry'] if 'magnitude' in x]
    photoAB = [float(x['magnitude'])
               for x in catalog[entry]['photometry']
               if 'magnitude' in x]
    photoABlowererrs = [float(x['e_lower_magnitude'])
                        if ('e_lower_magnitude' in x)
                        else (float(x['e_magnitude']) if 'e_magnitude' in x else 0.)
                        for x in catalog[entry]['photometry'] if 'magnitude' in x]
    photoABuppererrs = [float(x['e_upper_magnitude'])
                        if ('e_upper_magnitude' in x)
                        else (float(x['e_magnitude']) if 'e_magnitude' in x else 0.)
                        for x in catalog[entry]['photometry'] if 'magnitude' in x]
    photoband = [(x['band'] if 'band' in x else '')
                 for x in catalog[entry]['photometry'] if 'magnitude' in x]
    photoinstru = [(x['instrument'] if 'instrument' in x else '')
                   for x in catalog[entry]['photometry']
                   if 'magnitude' in x]
    photosource = [', '.join(str(j)
                             for j in sorted(int(i) for i in x['source'].split(',')))
                   for x in catalog[entry]['photometry']
                   if 'magnitude' in x]
    phototype = [(x['upperlimit'] if 'upperlimit' in x else False)
                 for x in catalog[entry]['photometry'] if 'magnitude' in x]
    photocorr = [('k' if 'kcorrected' in x else 'raw')
                 for x in catalog[entry]['photometry'] if 'magnitude' in x]

    photoutime = catalog[entry]['photometry'][0][
        'u_time'] if 'u_time' in catalog[entry]['photometry'][0] else 'MJD'
    hastimeerrs = (len(list(filter(None, phototimelowererrs))) and
                   len(list(filter(None, phototimeuppererrs))))
    hasABerrs = (len(list(filter(None, photoABlowererrs))) and
                 len(list(filter(None, photoABuppererrs))))
    tt = [
        ("Source ID(s)", "@src"),
        ("Epoch (" + photoutime + ")",
         "@x{1.11}" + ("<sub>-@xle{1}</sub><sup>+@xue{1}</sup>"
                       if hastimeerrs else ""))
    ]
    tt += [("Apparent Magnitude",
            "@y{1.111}" + ("<sub>-@lerr{1.11}</sub><sup>+@uerr{1.11}</sup>"
                           if hasABerrs else ""))]
    if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
        tt += [("Absolute Magnitude", "@yabs{1.111}" +
                ("<sub>-@lerr{1.11}</sub><sup>+@uerr{1.11}</sup>"
                 if hasABerrs else ""))]
    if len(list(filter(None, photoband))):
        tt += [("Band", "@desc")]
    if len(list(filter(None, photoinstru))):
        tt += [("Instrument", "@instr")]
    hover = HoverTool(tooltips=tt)

    min_y_range = 0.5 + max([
        (x + y) if not z else x
        for x, y, z in list(zip(photoAB, photoABuppererrs, phototype))
    ])
    max_y_range = -0.5 + min([
        (x - y) if not z else x
        for x, y, z in list(zip(photoAB, photoABlowererrs, phototype))
    ])

    p1 = Figure(
        title='Photometry for ' + event_name,
        active_drag='box_zoom',
        # sizing_mode = "scale_width",
        y_axis_label='Apparent Magnitude',
        tools=TOOLS_LIST,
        plot_width=485,
        plot_height=485,
        x_range=(min_x_range, max_x_range),
        y_range=(min_y_range, max_y_range),
        toolbar_location='above',
        toolbar_sticky=False)
    p1.xaxis.axis_label_text_font = 'futura'
    p1.yaxis.axis_label_text_font = 'futura'
    p1.xaxis.major_label_text_font = 'futura'
    p1.yaxis.major_label_text_font = 'futura'
    p1.xaxis.axis_label_text_font_size = '11pt'
    p1.yaxis.axis_label_text_font_size = '11pt'
    p1.xaxis.major_label_text_font_size = '8pt'
    p1.yaxis.major_label_text_font_size = '8pt'
    p1.title.align = 'center'
    p1.title.text_font_size = '16pt'
    p1.title.text_font = 'futura'

    min_x_date = astrotime(min_x_range, format='mjd').datetime
    max_x_date = astrotime(max_x_range, format='mjd').datetime

    p1.extra_x_ranges = {"gregorian date": Range1d(
        start=min_x_date, end=max_x_date)}
    p1.add_layout(
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
        p1.extra_x_ranges["time since max"] = Range1d(
            start=min_xm_range, end=max_xm_range)
        p1.add_layout(
            LinearAxis(
                axis_label="Time since max (" + dayframe + ")",
                major_label_text_font_size='8pt',
                major_label_text_font='futura',
                axis_label_text_font='futura',
                x_range_name="time since max",
                axis_label_text_font_size='11pt'),
            'above')

    if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
        min_y_absmag = min_y_range - distancemod
        max_y_absmag = max_y_range - distancemod
        p1.extra_y_ranges = {"abs mag": Range1d(
            start=min_y_absmag, end=max_y_absmag)}
        p1.add_layout(
            LinearAxis(
                axis_label="Absolute Magnitude",
                major_label_text_font_size='8pt',
                major_label_text_font='futura',
                axis_label_text_font='futura',
                y_range_name="abs mag",
                axis_label_text_font_size='11pt'),
            'right')
    p1.add_tools(hover)

    xs = []
    ys = []
    err_xs = []
    err_ys = []

    for x, y, xlowerr, xupperr, ylowerr, yupperr in list(
            zip(phototime, photoAB, phototimelowererrs, phototimeuppererrs,
                photoABlowererrs, photoABuppererrs)):
        xs.append(x)
        ys.append(y)
        err_xs.append((x - xlowerr, x + xupperr))
        err_ys.append((y - ylowerr, y + yupperr))

    bandset = list(set(photoband))
    bandsortlists = sorted(
        list(
            zip(
                list(map(bandgroupf, bandset)), list(
                    map(bandwavef, bandset)), list(
                        map(bandaliasf, bandset)))))
    bandset = list(filter(None, [i for (k, j, i) in bandsortlists]))

    sources = []
    corrects = ['raw', 'k']
    glyphs = [[] for x in range(len(corrects))]
    for ci, corr in enumerate(corrects):
        for band in bandset:
            bandname = bandaliasf(band)
            indb = [i for i, j in enumerate(photoband) if j == band]
            indt = [i for i, j in enumerate(phototype) if not j]
            indnex = [i for i, j in enumerate(phototimelowererrs)
                      if j == 0.]
            indyex = [i for i, j in enumerate(phototimelowererrs)
                      if j > 0.]
            indney = [i for i, j in enumerate(photoABuppererrs) if j == 0.]
            indyey = [i for i, j in enumerate(photoABuppererrs) if j > 0.]
            indc = [i for i, j in enumerate(photocorr) if j == corr]
            indne = set(indb).intersection(indt).intersection(
                indc).intersection(indney).intersection(indnex)
            indye = set(indb).intersection(indt).intersection(
                indc).intersection(set(indyey).union(indyex))

            noerrorlegend = bandname if len(indne) == 0 else ''

            data = dict(
                x=[phototime[i] for i in indne],
                y=[photoAB[i] for i in indne],
                lerr=[photoABlowererrs[i] for i in indne],
                uerr=[photoABuppererrs[i] for i in indne],
                desc=[photoband[i] for i in indne],
                instr=[photoinstru[i] for i in indne],
                src=[photosource[i] for i in indne])
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[
                    entry]:
                data['yabs'] = [photoAB[i] - distancemod for i in indne]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in indne]
                data['xue'] = [phototimeuppererrs[i] for i in indne]

            sources.append(ColumnDataSource(data))
            glyphs[ci].append(
                p1.circle(
                    'x',
                    'y',
                    source=sources[-1],
                    color=bandcolorf(band),
                    fill_color="white",
                    legend=noerrorlegend,
                    size=4).glyph)

            data = dict(
                x=[phototime[i] for i in indye],
                y=[photoAB[i] for i in indye],
                lerr=[photoABlowererrs[i] for i in indye],
                uerr=[photoABuppererrs[i] for i in indye],
                desc=[photoband[i] for i in indye],
                instr=[photoinstru[i] for i in indye],
                src=[photosource[i] for i in indye])
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[
                    entry]:
                data['yabs'] = [photoAB[i] - distancemod for i in indye]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in indye]
                data['xue'] = [phototimeuppererrs[i] for i in indye]

            sources.append(ColumnDataSource(data))
            glyphs[ci].append(
                p1.multi_line(
                    [err_xs[x] for x in indye], [
                        [ys[x], ys[x]] for x in indye
                    ],
                    color=bandcolorf(band)).glyph)
            glyphs[ci].append(
                p1.multi_line(
                    [[xs[x], xs[x]] for x in indye], [
                        err_ys[x] for x in indye
                    ],
                    color=bandcolorf(band)).glyph)
            glyphs[ci].append(
                p1.circle(
                    'x',
                    'y',
                    source=sources[-1],
                    color=bandcolorf(band),
                    legend=bandname,
                    size=4).glyph)

            upplimlegend = bandname if len(indye) == 0 and len(
                indne) == 0 else ''

            indt = [i for i, j in enumerate(phototype) if j]
            ind = set(indb).intersection(indt)
            data = dict(
                x=[phototime[i] for i in ind],
                y=[photoAB[i] for i in ind],
                lerr=[photoABlowererrs[i] for i in ind],
                uerr=[photoABuppererrs[i] for i in ind],
                desc=[photoband[i] for i in ind],
                instr=[photoinstru[i] for i in ind],
                src=[photosource[i] for i in ind])
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[
                    entry]:
                data['yabs'] = [photoAB[i] - distancemod for i in ind]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in ind]
                data['xue'] = [phototimeuppererrs[i] for i in ind]

            sources.append(ColumnDataSource(data))
            # Currently Bokeh doesn't support tooltips for
            # inverted_triangle, so hide an invisible circle behind for the
            # tooltip
            glyphs[ci].append(
                p1.circle(
                    'x', 'y', source=sources[-1], alpha=0.0, size=7).glyph)
            glyphs[ci].append(
                p1.inverted_triangle(
                    'x',
                    'y',
                    source=sources[-1],
                    color=bandcolorf(band),
                    legend=upplimlegend,
                    size=7).glyph)

            for gi, gly in enumerate(glyphs[ci]):
                if corr != 'raw':
                    glyphs[ci][gi].visible = False

    p1.legend.label_text_font = 'futura'
    p1.legend.label_text_font_size = '8pt'
    p1.legend.label_width = 20
    p1.legend.label_height = 14
    p1.legend.glyph_height = 14

    if any([x != 'raw' for x in photocorr]):
        photodicts = {}
        for ci, corr in enumerate(corrects):
            for gi, gly in enumerate(glyphs[ci]):
                photodicts[corr + str(gi)] = gly
        sdicts = dict(
            zip(['s' + str(x) for x in range(len(sources))], sources))
        photodicts.update(sdicts)
        photocallback = CustomJS(
            args=photodicts,
            code="""
            var show = 'all';
            if (cb_obj.get('value') == 'Raw') {
                show = 'raw';
            } else if (cb_obj.get('value') == 'K-Corrected') {
                show = 'k';
            }
            var viz = (show == 'all') ? true : false;
            var corrects = ["raw", "k"];
            for (c = 0; c < """ + str(len(corrects)) + """; c++) {
                for (g = 0; g < """ + str(len(glyphs[0])) + """; g++) {
                    if (show == 'all' || corrects[c] != show) {
                        eval(corrects[c] + g).attributes.visible = viz;
                    } else if (show != 'all' || corrects[c] == show) {
                        eval(corrects[c] + g).attributes.visible = !viz;
                    }
                }
            }
            for (s = 0; s < """ + str(len(sources)) + """; s++) {
                eval('s'+s).trigger('change');
            }
        """)
        photochecks = Select(
            title="Photometry to show:",
            value="Raw",
            options=["Raw", "K-Corrected", "All"],
            callback=photocallback)
    else:
        photochecks = ''

    return
