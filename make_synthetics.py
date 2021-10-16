from makeSynteric import make_synthetic_moment, new_make_synthetic
from myParser import parse_my_line
from myPlot import plot_data, plot_stream_wavelets, plot_green, plot_data_vs_synthetics, plot_data_vs_synthetics_wc, create_map
from WDMT_inv import correlate_phase_shift, set_GTG, set_station_weights, righthand_side, extract_parameters_mt, to_rtf, \
    to_xyz, check_fit
from gaussj import gaussian_jordan
from wavelet import wavelet_transform_j, inverse_wavelet_transform_j, next_pow2
from obspy.imaging.beachball import beachball
# import matplotlib.pyplot as plt
import os
import copy
import numpy as np
from scipy.signal import detrend
from obspy.signal.filter import bandpass, highpass
from mpl_toolkits import basemap
import matplotlib.pyplot as plt
import mplstereonet

args = parse_my_line()
print(args)
strike_dip_rake = [45, 45, 45, 1.0e24, 0]  # [0, 90, 0, 1.3e24, 0]#[24, 69, 47, 1.3e24, 0]
azimuths = [0, 45, 90, 135, 180, 225, 270, 315]
distance = [100, 350, 100, 350, 100, 350, 100, 350]
random_noise_amp = 0.05
freq_noise_amp = 0.00
freq_noise_freq = 0.2
###########################
# dont change anything below if you are not skillful


greens, synthetics = new_make_synthetic(args, strike_dip_rake, distance, azimuths)
greens, synthetics = correlate_phase_shift(greens, synthetics)

if not os.path.exists(os.path.join(args.output, "syn_data")):
    os.mkdir(os.path.join(args.output, "syn_data"))
stations = os.path.join(args.output, 'syn_data', 'synth_stations.txt')
fid = open(stations, "w")

for tr in synthetics:
    print(tr.stats.Zcor, tr.stats.Tcor, tr.stats.Rcor)
    tr_name = "%s.%s.%s.SAC" % (tr.stats.network, tr.stats.station, tr.stats.channel)
    tr.write(os.path.join(args.output, "syn_data", tr_name), format="SAC")
    if tr.stats.channel == "Z":
        fid.write("%s %s %f %f %f\n" % (tr.stats.station, tr.stats.network, tr.stats.stlo, tr.stats.stla, 0.0))

fid.close()

fid = open(os.path.join(args.output,"output", "results.txt"), "w")
fid.write("j\tfrqmin\tfreqmax\tMxx\tMyy\tMzz\tMxy\tMxz\tMyz\tMw\tMo\tDC\tCLVD\tISO\tVAR\tQUALITY\tstrike1\tdip1\trake1"
          "\tstrike2\tdip2\trake2\n") 

# synthetics.filter(type="bandpass", freqmin=.01, freqmax=0.05, corners=2)
# greens.filter(type="bandpass", freqmin=.01, freqmax=0.05, corners=2)


gfscale = 1e+20

(greens, synthetics) = set_station_weights(greens, synthetics)

ori_data = copy.deepcopy(synthetics)
ori_gree = copy.deepcopy(greens)

T_strikes = []
P_strikes = []
T_dips = []
P_dips = []
strike1 = []
dip1 = []
rake1 = []
strike2 = []
dip2 = []
rake2 =[]

for i in range(len(ori_data)):
    noise = np.sin(2 * np.pi * freq_noise_freq * np.arange(ori_data[i].stats.npts) * ori_data[i].stats.delta)
    noise = detrend(noise)
    rand_noise = np.random.normal(size=len(ori_data[i].data))
    rand_noise /= np.max(np.abs(rand_noise))
    rand_noise -= np.mean(rand_noise)
    rand_noise = detrend(rand_noise)
    # rand_noise = highpass(rand_noise, freq=.5, df=1 / ori_data[i].stats.delta, corners=4)

    maxi = np.max(ori_data[i].data)
    ori_data[i].data = ori_data[i].data + freq_noise_amp * maxi * noise + random_noise_amp * maxi * rand_noise
plot_green(args, greens)
plot_data(args, ori_data)
plot_stream_wavelets(args, ori_data)
create_map(args, "map.eps", ori_data)

for j in range(int(next_pow2(ori_data[0].stats.npts))):
    print("j = %d freqmin = %.3f freqmax = %.3f" % (
        j, 2 ** j / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta,
        2 ** (j + 2) / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta))

    for k in range(len(ori_data)):
        synthetics[k].data = wavelet_transform_j(ori_data[k].data, ori_data[k].stats.npts, 1 / ori_data[k].stats.delta,
                                                 j)
        synthetics[k].stats.npts = 2 ** j
    for k in range(len(ori_gree)):
        greens[k].data = wavelet_transform_j(ori_gree[k].data, ori_gree[k].stats.npts, 1 / ori_gree[k].stats.delta, j)
        greens[k].stats.npts = 2 ** j

    (AIV, B, AJ, W) = set_GTG(args, synthetics, greens)
    # Calculate Righthand Side
    B = righthand_side(args, synthetics, B, AJ, W)
    # OK as B value

    # Solve for MT
    M = gaussian_jordan(AIV, B)

    # if isotropic constrain, set Mzz
    if args.iso == "0":
        zz = -1 * (M[0] + M[1])
        M[5] = zz

    M *= -1.

    # *Convert deviatoric moment tensor to AKI convention*
    MTx = to_xyz(M, gfscale)
    MTr = to_rtf(M, gfscale)

    # sy = make_synthetic_moment(M * gfscale, greens)

    # synthetics, sy, variance, quality = check_fit(synthetics, sy)
    # variance = 200
    # quality = 100
    # Here compute Planes, Axes, Mo, Mw
    (Mo, Mw, Pdc, Pclvd, Pciso, EigVal, EigVec, T, N, P, np1, np2) = extract_parameters_mt(MTx, gfscale)

    M = [MTx[0][0], MTx[1][1], MTx[0][1], MTx[0][2], MTx[1][2], MTx[2][2]]
    sy = make_synthetic_moment(M, greens, gfscale)
    sy_t = copy.deepcopy(sy)
    synthetics_t = copy.deepcopy(synthetics)

    for n in range(len(sy_t)):
        sy_t[n].data = inverse_wavelet_transform_j(sy[n].data, ori_data[n].stats.npts, 1 / ori_data[n].stats.delta, j)
        synthetics_t[n].data = inverse_wavelet_transform_j(synthetics[n].data, ori_data[n].stats.npts,
                                                           1 / ori_data[n].stats.delta, j)
    synthetics_t = synthetics_t.taper(max_percentage=.1)
    name_t = "dataVSsynthetcs%d.pdf" % j
    name_wc = "dataVSsynthetcs_wc%d.pdf" % j
    plot_data_vs_synthetics(args, name_t, synthetics_t, sy_t)
    plot_data_vs_synthetics_wc(args, name_wc, synthetics, sy)

    # for iii in range(len(synthetics)):
    #     print(np.max(sy[iii].data), np.max(synthetics[iii].data))

    st, sy, variance, quality = check_fit(synthetics, sy)

    print("Mw = %.1f Mo = %.2e DC = %.1f%% CLVD = %.1f%% ISO = %.1f%% VAR = %.1f QUALITY = %d" % (
        Mw, Mo, Pdc, Pclvd, Pciso, variance, quality))
    print("strike1 = %.1f dip1 = %.1f rake1 = %.1f" % (np1[0], np1[1], np1[2]))
    print("strike2 = %.1f dip2 = %.1f rake2 = %.1f" % (np2[0], np2[1], np2[2]))
    fid.write("%d\t%.4f\t%.4f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.1f\t%.2e\t%.1f\t%.1f\t%.1f\t%.2f\t%d\t%.1f\t%.1f"
              "\t%.1f\t%.1f\t%.1f\t%.1f\n" %
              (j, 2 ** j / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta,
               2 ** (j + 2) / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta,
               MTx[0][0], MTx[1][1], MTx[2][2], MTx[0][1], MTx[0][2], MTx[2][2], Mw, Mo, Pdc, Pclvd, Pciso, variance,
               quality, np1[0], np1[1], np1[2], np2[0], np2[1], np2[2]))

    for kk in range(int(variance)):
        if quality >= 3:
            strike1.append(np1[0])
            dip1.append(np1[1])
            rake1.append(np1[1])
            strike2.append(np2[0])
            dip2.append(np2[1])
            rake1.append(np1[1])
            T_strikes.append(T.strike)
            T_dips.append(T.dip)
            P_strikes.append(P.strike)
            P_dips.append(P.dip)

fid.close()
fig = plt.figure()

ax = fig.add_subplot(111, projection='stereonet')
ax.pole(T_strikes, T_dips, c='r', label='T Pole of the Planes')
# ax.density_contourf(T_strikes, T_dips, measurement='poles', cmap='Reds')
# ax.set_title('Density coutour of the P', y=1.10, fontsize=15)
ax.grid()

# ax = fig.add_subplot(122, projection='stereonet')
ax.pole(P_strikes, P_dips, c='b', label='P Pole of the Planes')
# ax.density_contourf(P_strikes, P_dips, measurement='poles', cmap='Blues')
# ax.set_title('Density coutour of the P', y=1.10, fontsize=15)
ax.plane(strike1, dip1, c='k')
ax.plane(strike2, dip2, c='k')
ax.grid()
#plt.legend(loc="upper right")

plt.savefig(os.path.join(args.output, "output", "focal.eps"))
plt.show()




