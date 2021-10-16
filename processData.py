from obspy.core import UTCDateTime, Stream
from obspy.signal.rotate import rotate_ne_rt
from delaz import delaz


def station_list(data):
    out = []
    for tr in data:
        out.append(tr.stats.station)
    return sorted(list(set(out)))


def cut_window(args, data):
    tb = UTCDateTime(args.ori)
    tb = tb - float(args.pre)
    te = tb + float(args.len)
    data.trim(starttime=tb, endtime=te, pad=True, fill_value=0.0)
    return tb, te, data


def del_az_to_stream(data):
    for tr in data:
        # value not set as initialized
        if tr.stats.evla != -1000.0:
            dist, az, baz = delaz(tr.stats.evla, tr.stats.evlo, tr.stats.stla, tr.stats.stlo)
            tr.stats.dist = dist
            if baz < 0:
                baz += 360
            tr.stats.baz = baz
            if az < 0:
                az += 360
            tr.stats.az = az
        else:
            print("No station coordinates found for station ", tr.stats.station, "! Check station file.")
    return data


def remove_short_traces(data, tolerance, start_time, end_time):
    duration = end_time - start_time
    no_gap = Stream()

    for tr in data:
        observed = tr.stats.npts * tr.stats.delta * 100 / duration
        if observed >= tolerance:
            no_gap.append(tr)
    return no_gap


def sort_stream(data, keys):
    for i in range(len(keys)):
        if keys[i] == 'A':
            keys[i] = 'az'
        if keys[i] == 'D':
            keys[i] = 'dist'
        if keys[i] == 'C':
            keys[i] = 'station'

        for _i in keys[::-1]:
            data.traces.sort(key=lambda x: x.stats[_i], reverse=False)

    return data


def remove_21_component(data):
    tmp = Stream()
    stations = station_list(data)
    new_list = []
    for i in range(len(stations)):
        c = 0
        for j in range(len(data)):
            if data[j].stats.station == stations[i]:
                c += 1
            if c == 3:
                new_list.append(data[j].stats.station)
                break

    for i in range(len(new_list)):
        for j in range(len(data)):
            if data[j].stats.station == new_list[i]:
                tmp.append(data[j])

    return tmp


def remove_mean_trend(data):
    data.detrend()
    data.taper(type="cosine", max_percentage=.1)
    return data


def remove_trace(args, data, key):
    tmp = Stream()
    if key == "s" and args.set != "ALL":
        unavailable = []
        stations = args.set.split()
        for i in range(len(stations)):
            c = 0
            for tr in data:
                if tr.stats.station == stations[i]:
                    tmp.append(tr)
                    c += 1
                    if c == 3:
                        break
            if c == 0:
                unavailable.append(stations[i])
        if len(unavailable) != 0:
            print("Stations", unavailable, "not in stream")

    if key == "d":
        ran = args.range.split()
        for tr in data:
            if eval(ran[0]) <= tr.stats.dist <= eval(ran[1]):
                tmp.append(tr)

    if key == "a":
        ran = args.azi.split()
        for tr in data:
            if eval(ran[0]) <= tr.stats.az <= eval(ran[1]):
                tmp.append(tr)

    if key == 'r':
        for tr in data:
            be = True
            for station in args:
                if tr.stats.station == station:
                    be = False
            if be:
                tmp.append(tr)

    return tmp


def rotate_to_gcp(data):
    # begin loop over data stream
    for i in range(len(data) - 1):
        # split channel
        li0 = list(data[i + 0].stats['channel'])
        li1 = list(data[i + 1].stats['channel'])

        # check if station and part 1 of channel is identical and location
        if li0[0] == li1[0] and li0[1] == li1[1] \
                and data[i + 0].stats['station'] == data[i + 1].stats['station'] \
                and data[i + 0].stats['location'] == data[i + 1].stats['location']:

            rch = li0[0] + li0[1] + 'R'
            tch = li0[0] + li0[1] + 'T'

            # if yes 3 possibility: EN, NE , pass
            if li0[2] == "E" and li1[2] == "N":
                # baz
                baz = data[i].stats['baz']
                if data[i + 0].stats['npts'] == data[i + 1].stats['npts']:
                    # rotate 0-1
                    (data[i + 1].data, data[i + 0].data) = rotate_ne_rt(data[i + 1].data, data[i + 0].data, baz)
                    data[i + 0].stats['channel'] = tch
                    data[i + 1].stats['channel'] = rch
                    i += 1
                else:
                    print("Can't rotate ", data[i + 0].stats['station'], data[i + 0].stats['channel'], " and ",
                          data[i + 1].stats['station'], data[i + 1].stats['channel'])

            elif li0[2] == "N" and li1[2] == "E":
                # baz
                # baz = data[i].stats['baz']
                if data[i + 0].stats['npts'] == data[i + 1].stats['npts']:
                    (data[i + 0].data, data[i + 1].data) = rotate_ne_rt(data[i + 0].data, data[i + 1].data, baz)
                    data[i + 1].stats['channel'] = tch
                    data[i + 0].stats['channel'] = rch
                    i += 1
                else:
                    print("Can't rotate ", data[i + 0].stats['station'], data[i + 0].stats['channel'], " and ",
                          data[i + 1].stats['station'], data[i + 1].stats['channel'])

            else:
                pass

    return data


def remove_gap_data(data):
    gaps = data.get_gaps()
    if len(gaps) > 0:
        list_to_remove = [gaps[0]]
        for i in range(len(gaps)):
            ans = True
            for j in range(len(list_to_remove)):
                if gaps[i][0] == list_to_remove[j][0] and gaps[i][1] == list_to_remove[j][1] and gaps[i][2] == \
                        list_to_remove[j][2] and gaps[i][3] == \
                        list_to_remove[j][3]:
                    ans = False
            if ans:
                list_to_remove.append(gaps[i])
        for i in range(len(list_to_remove)):
            for tr in data.select(network=list_to_remove[i][0], station=list_to_remove[i][1],
                                  channel=list_to_remove[i][3]):
                data.remove(tr)
    return data


def clean_stream(args, data):
    data = remove_mean_trend(data)
    tb, te, data = cut_window(args, data)
    data = remove_gap_data(data)
    data = remove_short_traces(data, 100, tb, te)
    data = remove_21_component(data)
    data = del_az_to_stream(data)
    data = sort_stream(data, ['station', 'dist', 'az'])
    data = rotate_to_gcp(data)
    if args.deselect != "None":
        data = remove_trace(args.deselect.split(), data, "r")
    if args.range != "None":
        data = remove_trace(args, data, "d")
    return data
