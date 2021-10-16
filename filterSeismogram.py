import sys
import numpy as np
from obspy.io.sac import attach_paz, attach_resp
from obspy.signal import invsim
from processData import remove_trace


def remove_instrument_response(args, data):
    # prefilters
    if args.inst == 'Y':
        data.detrend()
        data.taper(type="hamming", max_percentage=.1)
        f = args.flim.split()
        f0 = eval(f[0])
        f1 = eval(f[1])
        f2 = eval(f[2])
        f3 = eval(f[3])
        to_remove = []  # station to purge if no Paz found

        for i in range(len(data)):
            # attach poles and zeros instrument
            if args.dva == '1':
                try:
                    if args.sim == 'PAZ':
                        attach_paz(data[i], data[i].stats.PAZ_file, todisp=False)
                    elif args.sim == 'RESP':
                        attach_resp(data[i], data[i].stats.RESP_file, todisp=True)
                except:
                    print("No appropriate PZs file found for station " + data[i].stats.station + data[i].stats.channel +
                          data[i].stats.network)
                    to_remove.append(data[i].stats.station)
            else:
                try:
                    if args.sim == 'PAZ':
                        attach_paz(data[i], data[i].stats.PAZ_file, tovel=True)
                    elif args.sim == 'RESP':
                        attach_resp(data[i], data[i].stats.RESP_file, tovel=False)
                except:
                    print("No appropriate PZs file found for station " + data[i].stats.station, data[i].stats.channel,
                          data[i].stats.network)
                    to_remove.append(data[i].stats.station)

        # remove stations if len(toPurge>0)
        if len(to_remove) > 0:
            data = remove_trace(data, to_remove, 'r')
            print("Check if station/channel/network/location of the PZs files and the same string within loaded "
                  "binary files ")
            print("do correspond. It may occour for instance that the headers strings of the waveform files (e.g. "
                  "sac, fseed) ")
            print("do not agrees with the same strings of the PZs name files. For instance the name of the network. ")
            print("If these strings do not correspond, modify the name of the PZs files or the header values of the "
                  "waveforms")
            print("You may also choose to remove this station using the option --purge (see help for details)")

        # now do remove
        for i in range(len(data)):
            # remove instrument to displacement
            #          st[i].data=detrend(st[i].data)
            data[i].data = invsim.simulate_seismometer(data[i].data, data[i].stats.sampling_rate,
                                                       paz_remove=data[i].stats.paz,
                                                       taper=True, taper_fraction=0.050,
                                                       pre_filt=(f0, f1, f2, f3))  # ,water_level=60.0)

            # from meters to centimeters
            data[i].data = data[i].data * 100

        return data


def find_j(frequency, n, fs):
    duration = n / fs
    if np.log2(3 * duration * frequency) > 0:
        print(np.log2(3 * duration * frequency))
        return np.int32(np.log2(3 * duration * frequency))
    else:
        return 0


def find_frequency(j, n, fs):
    duration = n / fs
    f1 = 2 ** j * 1 / 3 / duration
    f2 = 2 ** j * 4 / 3 / duration
    if f1 < 1 / duration:
        f1 = 1 / duration
    if f2 > fs / 2:
        f2 = fs / 2
    return f1, f2


def filer_data(args, data):
    hip = args.hip
    lop = args.lop
    bdp = args.bdp

    if hip != "0":
        elements = hip.split()
        corner = int(elements[0])
        frequency = eval(elements[1])
        data.filter(type="highpass", freq=frequency, corners=corner, zerophas=args.zerophase)
    if lop != "0":
        elements = lop.split()
        corner = int(elements[0])
        frequency = eval(elements[1])
        data.filter(type="lowpass", freq=frequency, corners=corner, zerophas=args.zerophase)
    if bdp != "0":
        elements = bdp.split()
        corner = int(elements[0])
        frequency1 = eval(elements[1])
        frequency2 = eval(elements[2])
        data.filter(type="bandpass", freqmin=frequency1, freqmax=frequency2, corners=corner, zerophase=args.zerophase)
    return data


def decimate_data(args, data, spin):
    dt_final = eval(args.DeltaInv)
    # if data is green (g), then find here decimation factor
    if spin == 'g':
        frs = 1 / data[0].stats.delta
        fri = 1 / eval(dt_final)
        c = int(frs / fri)
        ny = 1.0 / (2 * dt_final)
        for i in range(len(data)):
            data[i].filter("lowpass", freq=ny, corners=4, zerophase="False")
            data[i].decimate(c, strict_length=False, no_filter=True)

    elif spin == 'd':
        for i in range(len(data)):
            c = int((1 / data[i].stats.delta) / (1 / dt_final))
            ny = 1.0 / (2 * dt_final)
            data[i].filter("lowpass", freq=ny, corners=4, zerophase="False")
            data[i].decimate(c, strict_length=False, no_filter=True)

    else:
        print("Decimation option not recognized. EXIT!")
        sys.exit()

    # check if correct sampling rate, else interpolate
    if args.inter == 'Y' and spin == 'd':
        for i in range(len(data)):
            if data[i].stats.npts != int(eval(args.len)/eval(args.delta)):
                if args.war == 'Y':
                    print("Station " + str(data[i].stats.station) + "." + str(data[i].stats.channel) + " resampled from " +
                          str(data[i].stats.delta) + " to " + args.DeltaInv)
                data[i].resample(1 / dt_final)

    return data


def resample_data(args, data):
    data.resample(sampling_rate=1 / args.deltaInv, no_filter=False)
    return data


def filtering_data(args, data):
    data = remove_instrument_response(args, data)
    data = filer_data(args, data)
    data = decimate_data(args, data, 'd')
    return data
