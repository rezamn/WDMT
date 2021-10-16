from makeGreen import generate_greens
from readGreen import aquireGreens, reorderGreen, upGrStats
from processData import station_list
from delaz import find_coordnate
from obspy.core import Stream, Trace, UTCDateTime
import numpy as np
import sys


def make_synthetic(args, tensor, distances, azimuths):
    if len(tensor) == 5:
        # convert strike, dip and rake to radians and normalize the moment
        strike = tensor[0] * np.pi / 180
        dip = tensor[1] * np.pi / 180
        rake = tensor[2] * np.pi / 180
        m_fault = tensor[3]
        m_iso = tensor[4]

        M0 = m_fault
        M0 /= 1e+20
        m_iso /= 1e+20

        # mxx = -M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) + np.sin(2 * dip) * np.sin(rake) * np.sin(
        #     strike) ** 2)
        # myy = M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) - np.sin(2 * dip) * np.sin(rake) * np.cos(
        #     strike) ** 2)
        # mzz = M0 * (np.sin(2 * dip) * np.sin(rake))
        # myz = -M0 * (np.cos(dip) * np.cos(rake) * np.sin(strike) - np.cos(2 * dip) * np.sin(rake) * np.cos(strike))
        # mxz = - M0 * (np.cos(dip) * np.cos(rake) * np.cos(strike) + np.cos(2 * dip) * np.sin(rake) * np.sin(strike))
        # mxy = M0 * (np.sin(dip) * np.cos(rake) * np.cos(2 * strike) + 0.5 * np.sin(2 * dip) * np.sin(rake) * np.sin(
        #     2 * strike))
    else:
        print("the tensor format is not correct")
        sys.exit()

    st = Stream()
    for i in range(len(distances)):
        tr = Trace(np.arange(int(args.npts)))
        tr.stats.network = "XX"
        tr.stats.station = "station%2d" % i
        tr.stats.dist = distances[i]
        tr.stats.az = azimuths[i]
        tr.stats.starttime = UTCDateTime(args.ori)
        st.append(tr)

    greens = generate_greens(args, st)
    greens = aquireGreens(args, greens)
    grsta = station_list(greens)
    st = reorderGreen(greens, grsta)

    #st.filter(type="bandpass", freqmin=.01, freqmax=.05, corners=4)

    synt = Stream()
    staz_list = []

    # npts
    npts = st[0].stats.npts
    delta = st[0].stats.delta

    k = 1
    for dist, az in zip(distances, azimuths):
        print(dist, az)
        az = az * np.pi / 180
        str = az - strike
        print(str * 180 / np.pi)
        a = [0, 0, 0, 0, 0]
        a[0] = np.sin(2.0 * str) * np.cos(rake) * np.sin(dip) + 0.5 * np.cos(2.0 * str) * np.sin(rake) * np.sin(
            2.0 * dip)
        a[1] = np.cos(str) * np.cos(rake) * np.cos(dip) - np.sin(str) * np.sin(rake) * np.cos(2.0 * dip)
        a[2] = 0.5 * np.sin(rake) * np.sin(2.0 * dip)
        a[3] = np.cos(2.0 * str) * np.cos(rake) * np.sin(dip) - 0.5 * np.sin(2.0 * str) * np.sin(rake) * np.sin(
            2.0 * dip)
        a[4] = np.sin(strike) * np.cos(rake) * np.cos(dip) + np.cos(strike) * np.sin(rake) * np.cos(2.0 * dip)
        a[4] *= -1.0
        # a = [((mxx - myy) * np.cos(2 * az) / 2) + (mxy * np.sin(2 * az)), mxz * np.cos(az) + myz * np.sin(az),
        #      -(mxx + myy - 2 * mzz) / 6, ((mxx - myy) * np.sin(2 * az) / 2) - (mxy * np.cos(2 * az)),
        #      mxz * np.sin(az) - myz * np.cos(az)]

        # create a new stream of 10 greens function for given distance and azimuth
        tmp_greens = Stream()
        for tr in st:
            if abs(dist - tr.stats.dist) < .01:
                tmp_greens.append(tr)
                print(tr.stats.channel, tr.stats.dist)


        data = np.arange(npts)
        syn = Trace(data)
        lat_e, lon_e, depth = args.epi.split()
        lat_e = eval(lat_e)
        lon_e = eval(lon_e)
        lon_s, lat_s, baz = find_coordnate(lon_e, lat_e, az * 180 / np.pi, dist*1000.0)

        ###############
        # TAN Component
        # print(tmp_greens[0].stats.channel, tmp_greens[1].stats.channel)
        syn.data = M0 * (a[3] * tmp_greens[0].data + a[4] * tmp_greens[1].data)

        # apply Mo
        # syn.data = syn.data * (-1)

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'T'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz


        # add to synt stream
        synt.append(syn)

        ###############
        # RAD Component
        syn = Trace(data)
        # print(tmp_greens[2].stats.channel, tmp_greens[3].stats.channel, tmp_greens[4].stats.channel, tmp_greens[
        # 8].stats.channel)
        syn.data = M0 * (a[0] * tmp_greens[2].data + a[1] * tmp_greens[3].data + a[2] * tmp_greens[4].data) + m_iso * \
                   tmp_greens[8].data

        # apply Mo
        # syn.data = syn.data * (-1)

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'R'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz

        # add to synt stream
        synt.append(syn)

        ###############
        # VER Component
        syn = Trace(data)
        # print(tmp_greens[5].stats.channel, tmp_greens[6].stats.channel, tmp_greens[7].stats.channel,
        #       tmp_greens[9].stats.channel)
        syn.data = -1.0 * M0 * (
                a[0] * tmp_greens[5].data + a[1] * tmp_greens[6].data + a[2] * tmp_greens[7].data) + m_iso * \
                   tmp_greens[9].data

        # apply Mo
        syn.data = syn.data

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'Z'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz

        # add to synt stream
        synt.append(syn)
        k += 1

    #### for check

    # for i in range(len(synt)):
    #     if synt[i].stats.channel == 'R':
    #         staz_list.append(synt[i].stats.station)
    # k = 1
    # for l in range(len(staz_list)):
    #
    #     dati = np.arange(npts) * 0.0  # !!! initialize with floats
    #     az = azimuths[l] * np.pi / 180
    #     dist = distances[l]
    #
    #     ###############
    #     # TAN Component
    #     syn = Trace(dati)
    #
    #     syn.data = mxx * 0.5 * st[l * 10 + k * 0].data * np.sin(2 * az) \
    #                   - myy * 0.5 * st[l * 10 + k * 0].data * np.sin(2 * az) \
    #                   - mxy * 1.0 * st[l * 10 + k * 0].data * np.cos(2 * az) \
    #                   - mxz * 1.0 * st[l * 10 + k * 1].data * np.sin(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 1].data * np.cos(1 * az)
    #
    #         # apply Mo
    #     syn.data = syn.data * (-1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'T'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    #
    #     ###############
    #     # RAD Component
    #     syn = Trace(dati)
    #
    #
    #     syn.data = mxx * 1 / 6 * st[l * 10 + k * 4].data * (+1) \
    #                   - mxx * 0.5 * st[l * 10 + k * 2].data * np.cos(2 * az) \
    #                   + mxx * 1 / 3 * st[l * 10 + k * 8].data \
    #                   + myy * 1 / 6 * st[l * 10 + k * 4].data * (+1) \
    #                   + myy * 0.5 * st[l * 10 + k * 2].data * np.cos(2 * az) \
    #                   + myy * 1 / 3 * st[l * 10 + k * 8].data \
    #                   + mzz * 1 / 3 * st[l * 10 + k * 8].data \
    #                   - mzz * 1 / 3 * st[l * 10 + k * 4].data * (+1) \
    #                   - mxy * 1.0 * st[l * 10 + k * 2].data * np.sin(2 * az) \
    #                   + mxz * 1.0 * st[l * 10 + k * 3].data * np.cos(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 3].data * np.sin(1 * az)
    #
    #     # apply Mo
    #     syn.data = syn.data * (-1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'R'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    #
    #     ###############
    #     # VER Component
    #     syn = Trace(dati)
    #
    #     syn.data = mxx * 1 / 6 * st[l * 10 + k * 7].data \
    #                   - mxx * 0.5 * st[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
    #                   + mxx * 1 / 3 * st[l * 10 + k * 9].data \
    #                   + myy * 1 / 6 * st[l * 10 + k * 7].data \
    #                   + myy * 0.5 * st[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
    #                   + myy * 1 / 3 * st[l * 10 + k * 9].data \
    #                   + mzz * 1 / 3 * st[l * 10 + k * 9].data \
    #                   - mzz * 1 / 3 * st[l * 10 + k * 7].data \
    #                   - mxy * 1.0 * st[l * 10 + k * 5].data * (+1) * np.sin(2 * az) \
    #                   + mxz * 1.0 * st[l * 10 + k * 6].data * (+1) * np.cos(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 6].data * (+1) * np.sin(1 * az)
    #
    #     # apply Mo
    #     syn.data = syn.data * (+1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'Z'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    st = upGrStats(st, synt)

    return st, synt


def make_synthetic_moment(moment, greens, scale):
    gfscale = 1e+20
    if len(moment) == 6:
        mxx = moment[0]/scale
        myy = moment[1]/scale
        mxy = moment[2]/scale
        mxz = moment[3]/scale
        myz = moment[4]/scale
        mzz = moment[5]/scale
    else:
        print("the tensor format is not correct")
        sys.exit()
        # initialize
    synthetics = Stream()
    station_list = []

    # npts
    number_of_samples = greens[0].stats.npts
    delta = greens[0].stats.delta

    # aquire station list
    # print("here here", len(greens))
    for i in range(len(greens)):
        # print(i, greens[i].stats.channel, greens[i].stats.station)
        if greens[i].stats.channel == 'tss':
            station_list.append(greens[i].stats.station)

        # make synthetics Repeated 3 time the component loops over npts on syn.data
    k = 1
    for l in range(len(station_list)):
        tmp_data = np.arange(number_of_samples) * 0.0  # !!! initialize with floats
        az = greens[l * 10 + k].stats.az * np.pi / 180

        ###############
        # TAN Component
        syn = Trace(tmp_data)
        syn.data = mxx * 0.5 * greens[l * 10 + k * 0].data * np.sin(2 * az) \
                   - myy * 0.5 * greens[l * 10 + k * 0].data * np.sin(2 * az) \
                   - mxy * 1.0 * greens[l * 10 + k * 0].data * np.cos(2 * az) \
                   - mxz * 1.0 * greens[l * 10 + k * 1].data * np.sin(1 * az) \
                   + myz * 1.0 * greens[l * 10 + k * 1].data * np.cos(1 * az)

        # apply Mo
        syn.data = syn.data * (-1)

        # update stats
        syn.stats.station = greens[l * 10].stats.station
        syn.stats.channel = 'T'
        syn.stats.az = greens[l * 10].stats.az
        syn.stats.baz = greens[l * 10].stats.baz
        syn.stats.dist = greens[l * 10].stats.dist
        # syn.stats.gcarc = greens[l * 10].stats.gcarc
        syn.stats.evla = greens[l * 10].stats.evla
        syn.stats.evlo = greens[l * 10].stats.evlo
        syn.stats.stlo = greens[l * 10].stats.stlo
        syn.stats.stla = greens[l * 10].stats.stla
        syn.stats.delta = delta

        # add to synt stream
        synthetics.append(syn)

        ###############
        # RAD Component
        syn = Trace(tmp_data)
        syn.data = mxx * 1 / 6 * greens[l * 10 + k * 4].data * (+1) \
                   - mxx * 0.5 * greens[l * 10 + k * 2].data * np.cos(2 * az) \
                   + mxx * 1 / 3 * greens[l * 10 + k * 8].data \
                   + myy * 1 / 6 * greens[l * 10 + k * 4].data * (+1) \
                   + myy * 0.5 * greens[l * 10 + k * 2].data * np.cos(2 * az) \
                   + myy * 1 / 3 * greens[l * 10 + k * 8].data \
                   + mzz * 1 / 3 * greens[l * 10 + k * 8].data \
                   - mzz * 1 / 3 * greens[l * 10 + k * 4].data * (+1) \
                   - mxy * 1.0 * greens[l * 10 + k * 2].data * np.sin(2 * az) \
                   + mxz * 1.0 * greens[l * 10 + k * 3].data * np.cos(1 * az) \
                   + myz * 1.0 * greens[l * 10 + k * 3].data * np.sin(1 * az)

        # apply Mo
        syn.data = syn.data * (-1)

        # update stats
        syn.stats.station = greens[l * 10].stats.station
        syn.stats.channel = 'R'
        syn.stats.az = greens[l * 10].stats.az
        syn.stats.baz = greens[l * 10].stats.baz
        syn.stats.dist = greens[l * 10].stats.dist
        # syn.stats.gcarc = greens[l * 10].stats.gcarc
        syn.stats.evla = greens[l * 10].stats.evla
        syn.stats.evlo = greens[l * 10].stats.evlo
        syn.stats.stlo = greens[l * 10].stats.stlo
        syn.stats.stla = greens[l * 10].stats.stla
        syn.stats.delta = delta

        # add to synt stream
        synthetics.append(syn)

        ###############
        # VER Component
        syn = Trace(tmp_data)
        syn.data = mxx * 1 / 6 * greens[l * 10 + k * 7].data \
                   - mxx * 0.5 * greens[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
                   + mxx * 1 / 3 * greens[l * 10 + k * 9].data \
                   + myy * 1 / 6 * greens[l * 10 + k * 7].data \
                   + myy * 0.5 * greens[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
                   + myy * 1 / 3 * greens[l * 10 + k * 9].data \
                   + mzz * 1 / 3 * greens[l * 10 + k * 9].data \
                   - mzz * 1 / 3 * greens[l * 10 + k * 7].data \
                   - mxy * 1.0 * greens[l * 10 + k * 5].data * (+1) * np.sin(2 * az) \
                   + mxz * 1.0 * greens[l * 10 + k * 6].data * (+1) * np.cos(1 * az) \
                   + myz * 1.0 * greens[l * 10 + k * 6].data * (+1) * np.sin(1 * az)

        # apply Mo
        syn.data = syn.data * (+1)

        # update stats
        syn.stats.station = greens[l * 10].stats.station
        syn.stats.channel = 'Z'
        syn.stats.az = greens[l * 10].stats.az
        syn.stats.baz = greens[l * 10].stats.baz
        syn.stats.dist = greens[l * 10].stats.dist
        # syn.stats.gcarc = greens[l * 10].stats.gcarc
        syn.stats.evla = greens[l * 10].stats.evla
        syn.stats.evlo = greens[l * 10].stats.evlo
        syn.stats.stlo = greens[l * 10].stats.stlo
        syn.stats.stla = greens[l * 10].stats.stla
        syn.stats.delta = delta

        # add to synt stream
        synthetics.append(syn)

    return synthetics

def new_make_synthetic(args, tensor, distances, azimuths):
    if len(tensor) == 5:
        # convert strike, dip and rake to radians and normalize the moment
        strike = tensor[0] * np.pi / 180
        dip = tensor[1] * np.pi / 180
        rake = tensor[2] * np.pi / 180
        m_fault = tensor[3]
        m_iso = tensor[4]

        M0 = m_fault
        M0 /= 1e+20
        m_iso /= 1e+20

        mxx = -M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) + np.sin(2 * dip) * np.sin(rake) * np.sin(
            strike) ** 2) + m_iso
        myy = M0 * (np.sin(dip) * np.cos(rake) * np.sin(2 * strike) - np.sin(2 * dip) * np.sin(rake) * np.cos(
            strike) ** 2) + m_iso
        mzz = M0 * (np.sin(2 * dip) * np.sin(rake)) + m_iso
        myz = -M0 * (np.cos(dip) * np.cos(rake) * np.sin(strike) - np.cos(2 * dip) * np.sin(rake) * np.cos(strike))
        mxz = - M0 * (np.cos(dip) * np.cos(rake) * np.cos(strike) + np.cos(2 * dip) * np.sin(rake) * np.sin(strike))
        mxy = M0 * (np.sin(dip) * np.cos(rake) * np.cos(2 * strike) + 0.5 * np.sin(2 * dip) * np.sin(rake) * np.sin(
            2 * strike))
    else:
        print("the tensor format is not correct")
        sys.exit()

    st = Stream()
    for i in range(len(distances)):
        tr = Trace(np.arange(int(args.npts)))
        tr.stats.network = "XX"
        tr.stats.station = "station%2d" % i
        tr.stats.dist = distances[i]
        tr.stats.az = azimuths[i]
        tr.stats.starttime = UTCDateTime(args.ori)
        st.append(tr)

    greens = generate_greens(args, st)
    greens = aquireGreens(args, greens)
    grsta = station_list(greens)
    st = reorderGreen(greens, grsta)



    #st.filter(type="bandpass", freqmin=.01, freqmax=.05, corners=4)

    synt = Stream()
    staz_list = []

    # npts
    npts = st[0].stats.npts
    delta = st[0].stats.delta

    k = 1
    for dist, az in zip(distances, azimuths):
        print(dist, az)
        az = az * np.pi / 180
        str = az - strike
        # print(str * 180 / np.pi)
        # a = [0, 0, 0, 0, 0]
        # a[0] = np.sin(2.0 * str) * np.cos(rake) * np.sin(dip) + 0.5 * np.cos(2.0 * str) * np.sin(rake) * np.sin(
        #     2.0 * dip)
        # a[1] = np.cos(str) * np.cos(rake) * np.cos(dip) - np.sin(str) * np.sin(rake) * np.cos(2.0 * dip)
        # a[2] = 0.5 * np.sin(rake) * np.sin(2.0 * dip)
        # a[3] = np.cos(2.0 * str) * np.cos(rake) * np.sin(dip) - 0.5 * np.sin(2.0 * str) * np.sin(rake) * np.sin(
        #     2.0 * dip)
        # a[4] = np.sin(strike) * np.cos(rake) * np.cos(dip) + np.cos(strike) * np.sin(rake) * np.cos(2.0 * dip)
        # a[4] *= -1.0
        # a = [((mxx - myy) * np.cos(2 * az) / 2) + (mxy * np.sin(2 * az)), mxz * np.cos(az) + myz * np.sin(az),
        #      -(mxx + myy - 2 * mzz) / 6, ((mxx - myy) * np.sin(2 * az) / 2) - (mxy * np.cos(2 * az)),
        #      mxz * np.sin(az) - myz * np.cos(az)]

        # create a new stream of 10 greens function for given distance and azimuth
        tmp_greens = Stream()
        for tr in st:
            if abs(dist - tr.stats.dist) < .01:
                tmp_greens.append(tr)
                print(tr.stats.channel, tr.stats.dist)

        data = np.arange(npts)
        lat_e, lon_e, depth = args.epi.split()
        lat_e = eval(lat_e)
        lon_e = eval(lon_e)
        lon_s, lat_s, baz = find_coordnate(lon_e, lat_e, az * 180 / np.pi, dist)

        ###############
        # TAN Component
        # print(tmp_greens[0].stats.channel, tmp_greens[1].stats.channel)
        syn = Trace(data)
        syn.data = mxx * 0.5 * tmp_greens[0].data * np.sin(2 * az) \
                   - myy * 0.5 * tmp_greens[0].data * np.sin(2 * az) \
                   - mxy * 1.0 * tmp_greens[0].data * np.cos(2 * az) \
                   - mxz * 1.0 * tmp_greens[1].data * np.sin(1 * az) \
                   + myz * 1.0 * tmp_greens[1].data * np.cos(1 * az)

        # apply Mo
        syn.data = syn.data * (-1)

        # syn.data = M0 * (a[3] * tmp_greens[0].data + a[4] * tmp_greens[1].data)

        # apply Mo
        # syn.data = syn.data * (-1)

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'T'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz


        # add to synt stream
        synt.append(syn)

        ###############
        # RAD Component
        syn = Trace(data)
        # print(tmp_greens[2].stats.channel, tmp_greens[3].stats.channel, tmp_greens[4].stats.channel, tmp_greens[
        # 8].stats.channel)
        syn.data = mxx * 1 / 6 * tmp_greens[4].data * (+1) \
                   - mxx * 0.5 * tmp_greens[2].data * np.cos(2 * az) \
                   + mxx * 1 / 3 * tmp_greens[8].data \
                   + myy * 1 / 6 * tmp_greens[4].data * (+1) \
                   + myy * 0.5 * tmp_greens[2].data * np.cos(2 * az) \
                   + myy * 1 / 3 * tmp_greens[8].data \
                   + mzz * 1 / 3 * tmp_greens[8].data \
                   - mzz * 1 / 3 * tmp_greens[4].data * (+1) \
                   - mxy * 1.0 * tmp_greens[2].data * np.sin(2 * az) \
                   + mxz * 1.0 * tmp_greens[3].data * np.cos(1 * az) \
                   + myz * 1.0 * tmp_greens[3].data * np.sin(1 * az)

        # apply Mo
        syn.data = syn.data * (-1)

        # apply Mo
        # syn.data = syn.data * (-1)

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'R'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz

        # add to synt stream
        synt.append(syn)

        ###############
        # VER Component
        syn = Trace(data)
        # print(tmp_greens[5].stats.channel, tmp_greens[6].stats.channel, tmp_greens[7].stats.channel,
        #       tmp_greens[9].stats.channel)
        syn.data = mxx * 1 / 6 * tmp_greens[7].data \
                   - mxx * 0.5 * tmp_greens[5].data * (+1) * np.cos(2 * az) \
                   + mxx * 1 / 3 * tmp_greens[9].data \
                   + myy * 1 / 6 * tmp_greens[7].data \
                   + myy * 0.5 * tmp_greens[5].data * (+1) * np.cos(2 * az) \
                   + myy * 1 / 3 * tmp_greens[9].data \
                   + mzz * 1 / 3 * tmp_greens[9].data \
                   - mzz * 1 / 3 * tmp_greens[7].data \
                   - mxy * 1.0 * tmp_greens[5].data * (+1) * np.sin(2 * az) \
                   + mxz * 1.0 * tmp_greens[6].data * (+1) * np.cos(1 * az) \
                   + myz * 1.0 * tmp_greens[6].data * (+1) * np.sin(1 * az)

        # apply Mo
        syn.data = syn.data * (+1)

        # apply Mo
        syn.data = syn.data

        # update stats
        syn.stats.network = "XX"
        syn.stats.station = 'STA%02d' % k
        syn.stats.channel = 'Z'
        syn.stats['az'] = az * 180 / np.pi
        syn.stats['dist'] = dist
        syn.stats.delta = delta
        syn.stats['stla'] = lat_s
        syn.stats['stlo'] = lon_s
        syn.stats['evla'] = lat_e
        syn.stats['evlo'] = lon_e
        syn.stats['baz'] = baz

        # add to synt stream
        synt.append(syn)
        k += 1

    #### for check

    # for i in range(len(synt)):
    #     if synt[i].stats.channel == 'R':
    #         staz_list.append(synt[i].stats.station)
    # k = 1
    # for l in range(len(staz_list)):
    #
    #     dati = np.arange(npts) * 0.0  # !!! initialize with floats
    #     az = azimuths[l] * np.pi / 180
    #     dist = distances[l]
    #
    #     ###############
    #     # TAN Component
    #     syn = Trace(dati)
    #
    #     syn.data = mxx * 0.5 * st[l * 10 + k * 0].data * np.sin(2 * az) \
    #                   - myy * 0.5 * st[l * 10 + k * 0].data * np.sin(2 * az) \
    #                   - mxy * 1.0 * st[l * 10 + k * 0].data * np.cos(2 * az) \
    #                   - mxz * 1.0 * st[l * 10 + k * 1].data * np.sin(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 1].data * np.cos(1 * az)
    #
    #         # apply Mo
    #     syn.data = syn.data * (-1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'T'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    #
    #     ###############
    #     # RAD Component
    #     syn = Trace(dati)
    #
    #
    #     syn.data = mxx * 1 / 6 * st[l * 10 + k * 4].data * (+1) \
    #                   - mxx * 0.5 * st[l * 10 + k * 2].data * np.cos(2 * az) \
    #                   + mxx * 1 / 3 * st[l * 10 + k * 8].data \
    #                   + myy * 1 / 6 * st[l * 10 + k * 4].data * (+1) \
    #                   + myy * 0.5 * st[l * 10 + k * 2].data * np.cos(2 * az) \
    #                   + myy * 1 / 3 * st[l * 10 + k * 8].data \
    #                   + mzz * 1 / 3 * st[l * 10 + k * 8].data \
    #                   - mzz * 1 / 3 * st[l * 10 + k * 4].data * (+1) \
    #                   - mxy * 1.0 * st[l * 10 + k * 2].data * np.sin(2 * az) \
    #                   + mxz * 1.0 * st[l * 10 + k * 3].data * np.cos(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 3].data * np.sin(1 * az)
    #
    #     # apply Mo
    #     syn.data = syn.data * (-1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'R'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    #
    #     ###############
    #     # VER Component
    #     syn = Trace(dati)
    #
    #     syn.data = mxx * 1 / 6 * st[l * 10 + k * 7].data \
    #                   - mxx * 0.5 * st[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
    #                   + mxx * 1 / 3 * st[l * 10 + k * 9].data \
    #                   + myy * 1 / 6 * st[l * 10 + k * 7].data \
    #                   + myy * 0.5 * st[l * 10 + k * 5].data * (+1) * np.cos(2 * az) \
    #                   + myy * 1 / 3 * st[l * 10 + k * 9].data \
    #                   + mzz * 1 / 3 * st[l * 10 + k * 9].data \
    #                   - mzz * 1 / 3 * st[l * 10 + k * 7].data \
    #                   - mxy * 1.0 * st[l * 10 + k * 5].data * (+1) * np.sin(2 * az) \
    #                   + mxz * 1.0 * st[l * 10 + k * 6].data * (+1) * np.cos(1 * az) \
    #                   + myz * 1.0 * st[l * 10 + k * 6].data * (+1) * np.sin(1 * az)
    #
    #     # apply Mo
    #     syn.data = syn.data * (+1)
    #
    #     # update stats
    #     syn.stats.station = st[l * 10].stats.station
    #     syn.stats.channel = 'Z'
    #     syn.stats.az = az * 180 / np.pi
    #     # syn.stats.baz = st[l * 10].stats.baz
    #     syn.stats.dist = dist
    #     # syn.stats.gcarc = st[l * 10].stats.gcarc
    #     # syn.stats.evla = st[l * 10].stats.evla
    #     # syn.stats.evlo = st[l * 10].stats.evlo
    #     # syn.stats.stlo = st[l * 10].stats.stlo
    #     # syn.stats.stla = st[l * 10].stats.stla
    #     syn.stats.delta = delta
    #
    #     # add to synt stream
    #     synt.append(syn)
    st = upGrStats(st, synt)

    return st, synt
