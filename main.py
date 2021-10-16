from myParser import parse_my_line
from loadData import load_data, list_stream
from processData import station_list
from processData import clean_stream
from myPlot import plot_data, plot_green, plot_stream_wavelets, plot_data_vs_synthetics, plot_data_vs_synthetics_wc
from filterSeismogram import remove_instrument_response, decimate_data, filtering_data
from makeGreen import generate_greens
from readGreen import aquireGreens, reorderGreen, upGrStats
from wavelet import next_pow2, wavelet_transform_j, inverse_wavelet_transform_j
from WDMT_inv import new_correlate_phase_shift, set_GTG, set_station_weights, righthand_side, extract_parameters_mt, to_rtf, \
    to_xyz, check_fit, align_phase, new_correlate_shift
from gaussj import gaussian_jordan
from makeSynteric import make_synthetic_moment
import obspy.signal
from pandas import read_table
import mplstereonet
import os
import copy
import matplotlib.pyplot as plt

args = parse_my_line()
print(args)


if not os.path.exists(os.path.join(args.output, "output")):
    os.mkdir(os.path.join(args.output, "output"))
fid = open(os.path.join(args.output,"output", "results.txt"), "w")
data = load_data(args)

data = clean_stream(args, data)
list_stream(data)
data.detrend(type='linear')
data.taper(max_percentage=.1)
data = remove_instrument_response(args, data)
data = decimate_data(args, data, 'd')
# for tr in data:
#     print(tr.stats.npts)

greens = generate_greens(args, data)
greens = aquireGreens(args, greens)
grsta = station_list(greens)
greens = reorderGreen(greens, grsta)
greens = upGrStats(greens, data)
# for tr in greens:
#     print(tr.stats.npts)

(greens, data) = set_station_weights(greens, data)
# greens, data = correlate_phase_shift(greens, data)
# data = align_phase(data, args)
# for tr in data:
#     print(tr.stats.station, tr.stats.channel, tr.stats.Zcor, tr.stats.Tcor, tr.stats.Rcor, tr.stats.Vcor)


ori_data = copy.deepcopy(data)
ori_gree = copy.deepcopy(greens)

plot_green(args, greens)
plot_data(args, ori_data)
# plot_stream_wavelets(args, ori_data)

gfscale = 1e+20
T_strikes = []
P_strikes = []
T_dips = []
P_dips = []
plane1 = []
plane2 = []
for j in range(4, int(next_pow2(ori_data[0].stats.npts))):
    print("j = %d freqmin = %.3f freqmax = %.3f" % (
        j, 2 ** j / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta,
        2 ** (j + 2) / 3 / (2 ** next_pow2(ori_data[0].stats.npts)) / ori_data[0].stats.delta))

    for k in range(len(ori_data)):
        data[k].data = wavelet_transform_j(ori_data[k].data, ori_data[k].stats.npts, 1 / ori_data[k].stats.delta, j)
        data[k].data = inverse_wavelet_transform_j(data[k].data, ori_data[0].stats.npts, 1 / ori_data[0].stats.delta,
                                                   j)
    for k in range(len(ori_gree)):
        greens[k].data = wavelet_transform_j(ori_gree[k].data, ori_gree[k].stats.npts, 1 / ori_gree[k].stats.delta, j)
        greens[k].data = inverse_wavelet_transform_j(greens[k].data, ori_gree[0].stats.npts, 1 / ori_gree[0].stats.delta,
                                                   j)

    greens, data = new_correlate_phase_shift(greens, data)
    data = align_phase(data, args)
    # for tr in data:
    #     print(tr.stats.station, tr.stats.channel, tr.stats.Zcor, tr.stats.Tcor, tr.stats.Rcor, tr.stats.Vcor)

    for k in range(len(ori_data)):
        data[k].data = wavelet_transform_j(data[k].data, data[k].stats.npts, 1 / data[k].stats.delta, j)
        data[k].stats.npts = 2 ** j

    for k in range(len(ori_gree)):
        greens[k].data = wavelet_transform_j(ori_gree[k].data, ori_gree[k].stats.npts, 1 / ori_gree[k].stats.delta, j)
        greens[k].stats.npts = 2 ** j

    (AIV, B, AJ, W) = set_GTG(args, data, greens)
    # Calculate Righthand Side
    B = righthand_side(args, data, B, AJ, W)
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
    data_t = copy.deepcopy(data)

    for n in range(len(sy_t)):
        sy_t[n].data = inverse_wavelet_transform_j(sy[n].data, ori_data[n].stats.npts, 1/ori_data[n].stats.delta, j)
        data_t[n].data = inverse_wavelet_transform_j(data[n].data, ori_data[n].stats.npts, 1/ori_data[n].stats.delta, j)

    name_t = "dataVSsynthetcs%d.pdf" % j
    name_wc= "dataVSsynthetcs_wc%d.pdf" % j
    data_t, sy_t = new_correlate_shift(data_t, sy_t)
    data_t = align_phase(data_t, args)
    plot_data_vs_synthetics(args, name_t, data_t, sy_t, len(data_t)//3)
    plot_data_vs_synthetics_wc(args, name_wc, data, sy, len(data)//3)



    # for iii in range(len(synthetics)):
    #     print(np.max(sy[iii].data), np.max(synthetics[iii].data))

    data_t, sy_t, variance, quality = check_fit(data_t, sy_t)

    print("Mw = %.1f Mo = %.2e DC = %.1f%% CLVD = %.1f%% ISO = %.1f%% VAR = %.2f QUALITY = %d" % (
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
        if variance > 20:
            T_strikes.append(T.strike)
            T_dips.append(T.dip)
            P_strikes.append(P.strike)
            P_dips.append(P.dip)
            plane1.append(np1)
            plane2.append(np2)

fid.close()

fig = plt.figure()

ax = fig.add_subplot(111, projection='stereonet')

ax.pole(T_strikes, T_dips, c='r', label='T Pole of the Planes')
# ax.density_contourf(T_strikes, T_dips, measurement='poles', cmap='Reds')
# ax.set_title('Density coutour of the P', y=1.10, fontsize=15)
ax.grid()

ax.pole(P_strikes, P_dips, c='b', label='P Pole of the Planes')
# ax.plane(plane1[:][0], plane1[:][1], c='k')
# ax.plane(plane2[:][0], plane2[:][1], c='k')

# ax.density_contourf(P_strikes, P_dips, measurement='poles', cmap='Blues')
# ax.set_title('Density coutour of the P', y=1.10, fontsize=15)
ax.grid()

plt.show()