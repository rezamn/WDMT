import argparse


def parse_my_line():
    """
    Parse the arguments based on argparse

    Returns
    -------
    :return:
        args
    :rtype:
        :class:`argparse.Namespace`
    """
    #############################################################################################
    # -- options
    #

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='WDMT utility\n------------\n')
    ############################################################################################
    # The directories
    parser.add_argument('--info', default='None',
                        help='The input directory for the information structure. Default None')
    parser.add_argument('--input', default='None',
                        help='The input directory for the data structure. Default None')
    parser.add_argument('--output', default='None',
                        help='The output directory for inversion and results. Default None')

    ############################################################################################
    # The data type
    parser.add_argument('--fseed', default='None',
                        help='The input data is full seed format . Default None')
    parser.add_argument('--mseed', default='None',
                        help='the input data are mini seed format. Default None')
    parser.add_argument('--sac', default='None', help='The input data are in sac binary format. Default None')
    parser.add_argument('--asc', default='None', help='The input data are in sac ascii format. Default None')

    #############################################################################################
    # The instrument response type

    parser.add_argument('--paz', default='None',
                        help='The responses are in sac PAZ file format. Default None')
    parser.add_argument('--resp', default='None', help='The responses are in RESP file format. Default None')
    parser.add_argument('--inst', default='Y', help='Deconvolve the instrument response by paz or resp from the '
                                                    'waveforms.')
    ##############################################################################################
    # the station file name
    parser.add_argument('--sta', default='None', help='File for station file coordinates. Default None.')

    ##############################################################################################
    # the earthquake information
    parser.add_argument('--epi', default='None', help='Epicenter coordinate and depth: Lat Lon Depth.')
    parser.add_argument('--ori', default='None', help='The origin time of earthquake (e.g.: 2012-00-19T17:45:08).')

    ##############################################################################################
    # the data type, the start and the length
    parser.add_argument('--len', default='None', help='Length of signal in seconds After event origin Time.')
    parser.add_argument('--pre', default='0', help='Length of signal in seconds Before event origin Time.')
    parser.add_argument('--dva', default='1', help='Data type: 1(displacement); 2(velocity); 3(acceleration);         '
                                                   'Default = 1')
    parser.add_argument('--com', default='Earthquake',
                        help='comment string on earthquake information string. Default ="Earthquake"')

    ##############################################################################################
    # data analyzing parameter
    parser.add_argument('--range', default='None',
                        help='Min and Max distance range for stations to use in km. Default=None')
    parser.add_argument('--deselect', default='None', help='Station list to purege. Default=None')
    parser.add_argument('--set', default='All', help='Station list to use. Default=All')
    parser.add_argument('--azi', default='0 360', help='Azimuth range to select station. Default=0 360')

    ##############################################################################################
    # filter parameter
    parser.add_argument('--bdp', default='0',
                        help='Bandpass filter "corners fimn fmax". No Defaults. E.g.: "2 0.01 0.1"')
    parser.add_argument('--hip', default='0', help='Highpass filter "corners freq". No Defaults. E.g.: "2 0.01"')
    parser.add_argument('--lop', default='0', help='Lowpass filter "corners freq". No Defaults. E.g.: "2 0.1"')
    parser.add_argument('--zerophase', default='False',
                        help='Zero phase for high, low and bandpass. True/False. Defaul:False')
    parser.add_argument('--taper', default='0.1', help='cos Taper. If taper=-1 no taper is applied. Defaults=0.1')
    parser.add_argument('--sim', default='PAZ',
                        help='Remove instrument method: PAZ for poles and zeros, RESP, for RESP_ files. Default=PZs')
    parser.add_argument('--flim', default='0.002 0.005 0.5 1',
                        help='Corner frequency for deconvolution filtering. Defaults 0.002 0.005 0.5 1')
    parser.add_argument('--deci', default='None',
                        help='Decimation factor for sampling rate. Only integer decimation factor allowed. Default=None')
    parser.add_argument('--inter', default='Y',
                        help='Interpolate data to correct samplingrate if sampling not correct after decimation ['
                             'Y]/N. Artifact may occour. Warning about decimation use --war Y. Default=Y')
    ##############################################################################################
    # earth model and greens parameters
    parser.add_argument('--model', default='None',
                        help='Earth model for greenfunctions. Default=None. See README_MODEL.txt for details')
    parser.add_argument('--npts', default='1024', help='Number of points for greens. Power of 2. Default=1024')
    parser.add_argument('--delta', default='0.5', help='Sampling interval in seconds for greens. Default=0.5')
    parser.add_argument('--cpus', default='1',
                        help='Number of CPU available for greens computation. Min=1, Max=4. Default=1')
    parser.add_argument('--rvel', default='8', help='Reduction velocity. Recommendet value = 8km/s. Default=8')
    parser.add_argument('--war', default='Y', help='Warnings Y/[N]. Default=N')

    ##############################################################################################
    # Analysis
    parser.add_argument('--DeltaInv', default='0.5', help='Delta for data and greens for MT inversion. Default= 0.5 Hz')
    parser.add_argument('--iso', default='0', help='iso=0 [isotropic component set to 0, else iso=1. Default=0')
    parser.add_argument('--zcor', default='uniso', help='Fix Zcor for all 3 components (iso) or for each component (uniso), Default=uniso')

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    arg_run = parse_my_line()
    print(arg_run)
