import os
import sys
import glob
from pandas import read_table
from obspy import read


def load_data(args):
    """
    A wrapper for reading data from various data structure or database
    Parameters
    ----------
    :param args:
        The :class:`argparse.namespace` argument structure for passing the default values
        :type args: :class:`argparse.namespace`

    Returns
    -------
    :return: data
        The data stream
        :rtype data: :class:`obspy.Stream`
    """
    ##################################################################################
    # Check for data structure
    if not (os.path.exists(args.input) or os.path.exists(args.info)):
        print(r"The data structure is not found")
        sys.exit()
    else:
        if not (os.path.exists(os.path.join(args.input, "data"))):
            print(r"No 'data' folder is found in %s" % args.input)
            sys.exit()
        if not (os.path.exists(os.path.join(args.info, "RESP")) or os.path.exists(os.path.join(args.info, "PAZ"))):
            print(r"No response folder is found in %s" % args.input)
            sys.exit()

        ##################################################################################
    # A simple method for loading different data structure into data stream

    # full seed files
    if args.fseed != 'None':
        try:
            data = read(os.path.join(args.input, 'data', '*.seed'))
        except:
            print("The %s is not found" % args.fseed)
            sys.exit()

    # mini seed files
    if args.mseed != 'None':
        try:
            data = read(os.path.join(args.input, 'data', '*.mseed'))
        except:
            print("mini seed files not found in %s" % args.mseed)
            sys.exit()

    # sac binary files
    if args.sac != 'None':
        try:
            data = read(os.path.join(args.input, 'data', '*.sac'))
        except:
            try:
                data = read(os.path.join(args.input, 'data', '*.SAC'))
            except:
                print("sac files not found in %s" % args.sac)
                sys.exit()

    # sac ascii
    if args.asc != 'None':
        try:
            data = read(os.path.join(args.input, 'data', '*.sac'))
        except:
            print("sac ascii files not found in %s" % args.asc)
            sys.exit()

    data = reset_stations(data, args)
    data.sort()
    if args.paz != 'None':
        resp = list_resp(args, "PAZ")
        data = setting_traces(args, 'PAZ', data, resp)
    elif args.resp != 'None':
        resp = list_resp(args, "RESP")
        data = setting_traces(args, 'RESP', data, resp)
    else:
        print("No response file format is declared use --paz or --resp")
        sys.exit()

    return data


def reset_stations(data, args):
    # here already initialize hypo central information
    tmp = args.epi.split(' ')

    for i in range(len(data)):
        data[i].stats.RESP_file = "None"
        data[i].stats.PAZ_file = "None"
        data[i].stats.stla = -1000.0
        data[i].stats.stlo = -1000.0
        data[i].stats.stev = -1000.0
        data[i].stats.Zcor = -1234.5
        data[i].stats.Rcor = -1234.5
        data[i].stats.Tcor = -1234.5
        data[i].stats.Vcor = -1234.5
        data[i].stats.VR = -1234.5
        data[i].stats.evlo = eval(tmp[0])
        data[i].stats.evla = eval(tmp[1])
        data[i].stats.depth = eval(tmp[2])
        stations = read_table(os.path.join(args.info, "STATIONS", args.sta), delimiter=r'\s+', usecols=range(7),
                              names=["station", "network", "longitude", "latitude", "elevation", "code", "location"])
        for i in range(len(data)):
            sta = data[i].stats.station
            data[i].stats.station = sta[:5].replace(" ", "")
            channel = data[i].stats.channel
            data[i].stats.channel = channel.replace(" ", "")
            network = stations[(stations["station"] == sta)].values[0, :][1]
            data[i].stats.network = network

    return data


def setting_traces(args, mode, data, resp):
    for i in range(len(data)):
        name = data[i].stats.station
        loc = data[i].stats.location
        com = data[i].stats.channel
        net = data[i].stats.network
        if mode == "PAZ":
            resp_file = find_response_file(mode, resp, name, loc, net, com)
            data[i].stats.PAZ_file = resp_file
        elif mode == "RESP":
            res_file = find_response_file(mode, resp, name, loc, net, com)
            data[i].stats.RESP_file = res_file

        coordinate = find_station_coordination(args, net, name)
        data[i].stats.stla = coordinate[0]
        data[i].stats.stlo = coordinate[1]
        data[i].stats.stev = coordinate[2]
    return data


def find_station_coordination(args, network, station):
    # station format should be in format of sta net lan lot elv
    #                                   ex.  CHTR IR 123.4 34.7 1200
    tmp = read_table(os.path.join(args.info, "STATIONS", args.sta), delimiter=r'\s+', header=None, usecols=range(5),
                     names=['station', 'network', 'longitude', 'latitude', 'elevation'])
    tmp = tmp[(tmp['station'] == station) & (tmp['network'] == network)]
    if tmp.size != 0:
        return tmp.values[0, 2:].tolist()  # returns longitude, latitude and elevation
    else:
        return 0, 0, 0


def list_resp(args, mode):
    resp = []
    if mode == "RESP":
        for file in glob.glob(os.path.join(args.info, "RESP", "RESP*")):
            resp.append(file)
    elif mode == "PAZ":
        for file in glob.glob(os.path.join(args.info, "PAZ", "SAC_PZs_*")):
            resp.append(file)
        if len(resp) == 0:
            for file in glob.glob(os.path.join(args.info, "*SAC.PZs.*")):
                resp.append(file)

    return resp


def find_response_file(mode, resp, sta, loc, net, comp):
    if mode == "PAZ":
        for i in range(len(resp)):
            tmp = resp[i].split(os.sep)[-1].split("_")
            if len(tmp) > 1 and tmp[0] == "SAC" and tmp[1] == "PZs" and tmp[2] == net and tmp[3] == sta and tmp[
                4] == comp:
                return resp[i]
            else:
                tmp = resp[i].split(os.sep)[-1].split(".")
                if len(tmp) > 1 and tmp[0] == "SAC" and tmp[1] == "PZs" and tmp[2] == net and tmp[3] == sta and tmp[
                    4] == comp:
                    return resp[i]
    if mode == "RESP":
        for i in range(len(resp)):
            tmp = resp[i].split(os.sep)[-1].split(".")
            if tmp[0] == "RESP" and tmp[1] == net and tmp[2] == sta and tmp[3] == loc and tmp[4] == comp:
                return resp[i]


def list_stream(st):
    for i in range(len(st)):
        print(st[i].stats.network, st[i].stats.station, st[i].stats.channel, st[i].stats.location,
              st[i].stats.stla, st[i].stats.stlo, st[i].stats.evla, st[i].stats.evlo)#, st[i].stats.starttime, st[i].stats.endtime, st[i].stats.PAZ_file)
