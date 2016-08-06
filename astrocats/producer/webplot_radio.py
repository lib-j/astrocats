"""Radio Plotting Methods associatd with WebCat Script.
"""


def plot_radio(catalog, entry):
    phototime = [float(x['time']) for x in catalog[entry]['photometry']
                 if 'fluxdensity' in x]
    phototimelowererrs = [float(x['e_lower_time'])
                          if ('e_lower_time' in x and 'e_upper_time' in x)
                          else (float(x['e_time'])
                                if 'e_time' in x else 0.)
                          for x in catalog[entry]['photometry']
                          if 'fluxdensity' in x]
    phototimeuppererrs = [
        float(x['e_upper_time'])
        if ('e_lower_time' in x and 'e_upper_time' in x) in x else
        (float(x['e_time']) if 'e_time' in x else 0.)
        for x in catalog[entry]['photometry'] if 'fluxdensity' in x
    ]
    photofd = [
        float(x['fluxdensity']) if 'e_fluxdensity' not in x else
        (float(x['fluxdensity']) if
         (float(x['fluxdensity']) > RADIO_SIGMA * float(x['e_fluxdensity']))
         else round_sig(
             RADIO_SIGMA * float(x['e_fluxdensity']),
             sig=get_sig_digits(x['e_fluxdensity'])))
        for x in catalog[entry]['photometry'] if 'fluxdensity' in x
    ]
    photofderrs = [(float(x['e_fluxdensity'])
                    if 'e_fluxdensity' in x else 0.)
                   for x in catalog[entry]['photometry']
                   if 'fluxdensity' in x]
    photoufd = [(x['u_fluxdensity'] if 'fluxdensity' in x else '')
                for x in catalog[entry]['photometry']
                if 'fluxdensity' in x]
    photofreq = [(x['frequency'] if 'fluxdensity' in x else '')
                 for x in catalog[entry]['photometry']
                 if 'fluxdensity' in x]
    photoufreq = [(x['u_frequency'] if 'fluxdensity' in x else '')
                  for x in catalog[entry]['photometry']
                  if 'fluxdensity' in x]
    photoinstru = [(x['instrument'] if 'instrument' in x else '')
                   for x in catalog[entry]['photometry']
                   if 'fluxdensity' in x]
    photosource = [', '.join(
        str(j)
        for j in sorted(
            int(i)
            for i in catalog[entry]['photometry'][x]['source'].split(',')))
                   for x, y in enumerate(catalog[entry]['photometry'])
                   if 'fluxdensity' in y]
    phototype = [
        (True if 'upperlimit' in x or
         RADIO_SIGMA * float(x['e_fluxdensity']) >= float(x['fluxdensity'])
         else False) for x in catalog[entry]['photometry']
        if 'fluxdensity' in x
    ]

    photoutime = catalog[entry]['photometry'][0][
        'u_time'] if 'u_time' in catalog[entry]['photometry'][0] else 'MJD'
    if distancemod:
        dist = (10.0**(1.0 + 0.2 * distancemod) * un.pc).cgs.value
        areacorr = 4.0 * pi * dist**2.0 * ((1.0e-6 * un.jansky).cgs.value)

    x_buffer = 0.1 * (
        max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0

    hastimeerrs = (len(list(filter(None, phototimelowererrs))) and
                   len(list(filter(None, phototimeuppererrs))))
    hasfderrs = len(list(filter(None, photofderrs)))
    tt = [
        ("Source ID(s)", "@src"),
        ("Epoch (" + photoutime + ")",
         "@x{1.11}" + ("<sub>-@xle{1}</sub><sup>+@xue{1}</sup>"
                       if hastimeerrs else ""))
    ]
    tt += [("Flux Density (" + photoufd[0] + ")",
            "@y{1.11}" + ("&nbsp;±&nbsp;@err{1.11}" if hasfderrs else ""))]
    if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
        tt += [("Iso. Lum. (ergs s⁻¹)", "@yabs" + ("&nbsp;±&nbsp;@abserr"
                                                   if hasfderrs else ""))]
    if len(list(filter(None, photofreq))):
        tt += [("Frequency (" + photoufreq[0] + ")", "@desc")]
    if len(list(filter(None, photoinstru))):
        tt += [("Instrument", "@instr")]
    hover = HoverTool(tooltips=tt)

    if photoavail:
        x_range = p1.x_range
    else:
        x_range = (min_x_range, max_x_range)
    min_y_range = min([x - y for x, y in list(zip(photofd, photofderrs))])
    max_y_range = max([x + y for x, y in list(zip(photofd, photofderrs))])
    [min_y_range, max_y_range] = [
        min_y_range - 0.1 * (max_y_range - min_y_range),
        max_y_range + 0.1 * (max_y_range - min_y_range)
    ]

    p3 = Figure(
        title='Radio Observations of ' + event_name,
        active_drag='box_zoom',
        # sizing_mode = "scale_width",
        y_axis_label='Flux Density (µJy)',
        tools=TOOLS_LIST,
        plot_width=485,
        plot_height=485,
        x_range=x_range,
        y_range=(min_y_range, max_y_range),
        toolbar_location='above',
        toolbar_sticky=False)
    p3.xaxis.axis_label_text_font = 'futura'
    p3.yaxis.axis_label_text_font = 'futura'
    p3.xaxis.major_label_text_font = 'futura'
    p3.yaxis.major_label_text_font = 'futura'
    p3.xaxis.axis_label_text_font_size = '11pt'
    p3.yaxis.axis_label_text_font_size = '11pt'
    p3.xaxis.major_label_text_font_size = '8pt'
    p3.yaxis.major_label_text_font_size = '8pt'
    p3.title.align = 'center'
    p3.title.text_font_size = '16pt'
    p3.title.text_font = 'futura'

    min_x_date = astrotime(min_x_range, format='mjd').datetime
    max_x_date = astrotime(max_x_range, format='mjd').datetime

    p3.extra_x_ranges = {"gregorian date": Range1d(
        start=min_x_date, end=max_x_date)}
    p3.add_layout(
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
        p3.extra_x_ranges["time since max"] = Range1d(
            start=min_xm_range, end=max_xm_range)
        p3.add_layout(
            LinearAxis(
                axis_label="Time since max (" + dayframe + ")",
                major_label_text_font_size='8pt',
                major_label_text_font='futura',
                axis_label_text_font='futura',
                x_range_name="time since max",
                axis_label_text_font_size='11pt'),
            'above')

    if distancemod:
        min_y_absmag = min_y_range * areacorr * (1.0 * un.GHz).cgs.value
        max_y_absmag = max_y_range * areacorr * (1.0 * un.GHz).cgs.value
        # [min_y_absmag, max_y_absmag] = [min_y_absmag - 0.1 *
        #                                 (max_y_absmag - min_y_absmag), max_y_absmag + 0.1 * (max_y_absmag - min_y_absmag)]
        p3.extra_y_ranges = {"abs mag": Range1d(
            start=min_y_absmag, end=max_y_absmag)}
        p3.add_layout(
            LinearAxis(
                axis_label="Isotropic Luminosity at 1 GHz (ergs s⁻¹)",
                major_label_text_font_size='8pt',
                major_label_text_font='futura',
                axis_label_text_font='futura',
                y_range_name="abs mag",
                axis_label_text_font_size='11pt'),
            'right')
        p3.yaxis[1].formatter.precision = 1
    p3.add_tools(hover)

    xs = []
    ys = []
    err_xs = []
    err_ys = []

    for x, y, xlowerr, xupperr, yerr in list(
            zip(phototime, photofd, phototimelowererrs, phototimeuppererrs,
                photofderrs)):
        xs.append(x)
        ys.append(y)
        err_xs.append((x - xlowerr, x + xupperr))
        err_ys.append((y - yerr, y + yerr))

    freqset = [str(y) for y in sorted([float(x) for x in set(photofreq)])]
    frequnit = photoufreq[0] if photoufreq else ''

    for freq in freqset:
        indb = [i for i, j in enumerate(photofreq) if j == freq]
        indt = [i for i, j in enumerate(phototype) if not j]
        # Should always have upper error if have lower error.
        indnex = [i for i, j in enumerate(phototimelowererrs) if j == 0.]
        indyex = [i for i, j in enumerate(phototimelowererrs) if j > 0.]
        indney = [i for i, j in enumerate(photofderrs) if j == 0.]
        indyey = [i for i, j in enumerate(photofderrs) if j > 0.]
        indne = set(indb).intersection(indt).intersection(
            indney).intersection(indnex)
        indye = set(indb).intersection(indt).intersection(
            set(indyey).union(indyex))

        freqlabel = str(freq) + " " + frequnit

        noerrorlegend = freqlabel if len(indye) == 0 and len(
            indne) > 0 else ''

        data = dict(
            x=[phototime[i] for i in indne],
            y=[photofd[i] for i in indne],
            err=[photofderrs[i] for i in indne],
            desc=[photofreq[i] for i in indne],
            instr=[photoinstru[i] for i in indne],
            src=[photosource[i] for i in indne])
        if distancemod:
            data['yabs'] = [str(
                round_sig(
                    photofd[i] * (areacorr * float(freq) * ((1.0 * un.GHz
                                                             ).cgs.value)),
                    sig=3)) for i in indne]
            data['abserr'] = [str(
                round_sig(
                    photofderrs[i] * (areacorr * float(freq) * ((
                        1.0 * un.GHz).cgs.value)),
                    sig=3)) for i in indne]
        if hastimeerrs:
            data['xle'] = [phototimelowererrs[i] for i in indne]
            data['xue'] = [phototimeuppererrs[i] for i in indne]

        source = ColumnDataSource(data)
        p3.circle(
            'x',
            'y',
            source=source,
            color=radiocolorf(freq),
            fill_color="white",
            legend=noerrorlegend,
            size=4)

        yeserrorlegend = freqlabel if len(indye) > 0 else ''

        data = dict(
            x=[phototime[i] for i in indye],
            y=[photofd[i] for i in indye],
            err=[photofderrs[i] for i in indye],
            desc=[photofreq[i] for i in indye],
            instr=[photoinstru[i] for i in indye],
            src=[photosource[i] for i in indye])
        if distancemod:
            data['yabs'] = [str(
                round_sig(
                    photofd[i] * (areacorr * float(freq) * ((1.0 * un.GHz
                                                             ).cgs.value)),
                    sig=3)) for i in indye]
            data['abserr'] = [str(
                round_sig(
                    photofderrs[i] * (areacorr * float(freq) * ((
                        1.0 * un.GHz).cgs.value)),
                    sig=3)) for i in indye]
        if hastimeerrs:
            data['xle'] = [phototimelowererrs[i] for i in indye]
            data['xue'] = [phototimeuppererrs[i] for i in indye]

        source = ColumnDataSource(data)
        p3.multi_line(
            [err_xs[x] for x in indye], [[ys[x], ys[x]] for x in indye],
            color=radiocolorf(freq))
        p3.multi_line(
            [[xs[x], xs[x]] for x in indye], [err_ys[x] for x in indye],
            color=radiocolorf(freq))
        p3.circle(
            'x',
            'y',
            source=source,
            color=radiocolorf(freq),
            legend=yeserrorlegend,
            size=4)

        upplimlegend = freqlabel if len(indye) == 0 and len(
            indne) == 0 else ''

        indt = [i for i, j in enumerate(phototype) if j]
        ind = set(indb).intersection(indt)
        data = dict(
            x=[phototime[i] for i in ind],
            y=[photofd[i] for i in ind],
            err=[photofderrs[i] for i in ind],
            desc=[photofreq[i] for i in ind],
            instr=[photoinstru[i] for i in ind],
            src=[photosource[i] for i in ind])
        if distancemod:
            data['yabs'] = [str(
                round_sig(
                    photofd[i] * (areacorr * float(freq) * ((1.0 * un.GHz
                                                             ).cgs.value)),
                    sig=3)) for i in ind]
            data['abserr'] = [str(
                round_sig(
                    photofderrs[i] * (areacorr * float(freq) * ((
                        1.0 * un.GHz).cgs.value)),
                    sig=3)) for i in ind]
        if hastimeerrs:
            data['xle'] = [phototimelowererrs[i] for i in ind]
            data['xue'] = [phototimeuppererrs[i] for i in ind]

        source = ColumnDataSource(data)
        # Currently Bokeh doesn't support tooltips for inverted_triangle,
        # so hide an invisible circle behind for the tooltip
        p3.circle('x', 'y', source=source, alpha=0.0, size=7)
        p3.inverted_triangle(
            'x',
            'y',
            source=source,
            color=radiocolorf(freq),
            legend=upplimlegend,
            size=7)

    p3.legend.label_text_font = 'futura'
    p3.legend.label_text_font_size = '8pt'
    p3.legend.label_width = 20
    p3.legend.label_height = 14
    p3.legend.glyph_height = 14
