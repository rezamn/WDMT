from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.basemap import Basemap
import numpy as np
from wavelet import wavelet_transform, inverse_wavelet_transform_j
import os


def plot_data(args, st, spp=8):
    with PdfPages(os.path.join(args.output, 'output', 'seis.pdf')) as pdf:
        for k in range(0, len(st), spp * 3):
            fig, axs = plt.subplots(spp, 3, sharex='col', sharey='row')

            for n in range(spp):
                for i in range(3):
                    axs[n, i].spines['right'].set_visible(False)
                    axs[n, i].spines['left'].set_visible(False)
                    axs[n, i].spines['top'].set_visible(False)
                    axs[n, i].spines['bottom'].set_visible(False)
                    axs[n, i].set_yticks([])
                    axs[n, i].set_xticks([])
                    axs[n, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                          labelsize=6)
                    if n == spp - 1:
                        axs[n, i].spines['bottom'].set_visible(True)
                        axs[n, i].tick_params(axis='x', which='both', bottom=True, labelsize=6)
                        axs[n, i].xaxis.set_major_locator(MultipleLocator(100))
                        axs[n, i].xaxis.set_minor_locator(MultipleLocator(50))
                    else:
                        axs[n, i].tick_params(labelsize=6)

            for i in range(spp - 1, -1, -1):
                for j in range(3):
                    if i == 0:
                        axs[0, j].set_title(st[(spp - 1 - (spp - 1)) * 3 + j + k].stats.channel[-1])
                    if (spp - 1 - i) * 3 + j + k < len(st):
                        t = np.arange(st[(spp - 1 - i) * 3 + j + k].stats.npts) * st[
                            (spp - 1 - i) * 3 + j + k].stats.delta
                        axs[i, j].plot(t, st[(spp - 1 - i) * 3 + j + k].data, linewidth=.01, color="blue")
                        plt.text(0.9, 0.9, "max = %.2e" % np.max(np.abs(st[(spp - 1 - i) * 3 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=4, color='red')
                        # print(st[(9 - i) * 3 + j + k].stats.npts)
                        if j == 0:
                            plt.text(0.1, 0.9, "dist = %.1f\naz = %.1f" % (
                                st[(spp - 1 - i) * 3 + j + k].stats.dist, st[(spp - 1 - i) * 3 + j + k].stats.az),
                                     ha='center',
                                     va='center', transform=axs[i, j].transAxes,
                                     fontsize=4)

                            axs[i, j].set_ylabel(
                                st[(spp - 1 - i) * 3 + j + k].stats.network + "." + st[
                                    (spp - 1 - i) * 3 + j + k].stats.station,
                                fontsize=6, rotation=0, labelpad=30)

            pdf.savefig(fig, dpi=300)
            plt.close()


def plot_green(args, st):
    sss = 0
    with PdfPages(os.path.join(args.output, 'output', 'greens.pdf')) as pdf:
        for k in range(0, len(st), 50):
            fig, axs = plt.subplots(5, 10, sharex='col', sharey='row')
            sss += 1
            print(sss)

            for n in range(5):
                for i in range(10):
                    axs[n, i].spines['right'].set_visible(False)
                    axs[n, i].spines['left'].set_visible(False)
                    axs[n, i].spines['top'].set_visible(False)
                    axs[n, i].spines['bottom'].set_visible(False)
                    axs[n, i].set_yticks([])
                    axs[n, i].set_xticks([])
                    axs[n, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                          labelsize=4)
                    if n == 4:
                        axs[n, i].spines['bottom'].set_visible(True)
                        axs[n, i].tick_params(axis='x', which='both', bottom=True, labelsize=4)
                        axs[n, i].xaxis.set_major_locator(MultipleLocator(100))
                        axs[n, i].xaxis.set_minor_locator(MultipleLocator(50))
                    else:
                        axs[n, i].tick_params(labelsize=4)

            for i in range(4, -1, -1):
                for j in range(10):
                    if i == 0:
                        axs[0, j].set_title(st[(4 - 4) * 10 + j + k].stats.channel)
                    if (4 - i) * 10 + j + k < len(st):
                        t = np.arange(st[(4 - i) * 10 + j + k].stats.npts) * st[(4 - i) * 10 + j + k].stats.delta
                        axs[i, j].plot(t, st[(4 - i) * 10 + j + k].data, linewidth=.01, color="gray")
                        plt.text(0.9, 0.9, "max = %.1e" % np.max(np.abs(st[(4 - i) * 10 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=3, color='red')
                        # print(st[(9 - i) * 3 + j + k].stats.npts)
                        if j == 0:
                            plt.text(0.1, 0.9, "dist = %.1f\naz = %.1f" % (
                                st[(4 - i) * 10 + j + k].stats.dist, st[(4 - i) * 10 + j + k].stats.az), ha='center',
                                     va='center', transform=axs[i, j].transAxes,
                                     fontsize=3)

                            axs[i, j].set_ylabel("%s.%s" % (st[(4 - i) * 10 + j + k].stats.network,
                                                            st[(4 - i) * 10 + j + k].stats.station),
                                                 fontsize=6, rotation=0, labelpad=30)

            pdf.savefig(fig, dpi=300)
            plt.close()


def plot_stream_wavelets(args, data):
    with PdfPages(os.path.join(args.output, 'output', 'wavelets.pdf')) as pdf:
        for tr in data:
            fig, axs = plt.subplots(7, 2, sharex='col')
            fig.suptitle("%s.%s Channel:%s" % (tr.stats.network, tr.stats.station, tr.stats.channel))

            # max_amplitude = np.max(np.abs(tr.data))
            # ylim = [-max_amplitude * 1.1, max_amplitude * 1.1]

            t = np.arange(tr.stats.npts) * tr.stats.delta
            axs[0, 0].plot(t, tr.data, linewidth=.01, color="gray")
            # axs[0, 0].set_ylim(ylim)

            wc = wavelet_transform(tr.data, tr.stats.npts, 1 / tr.stats.delta)
            vmax = []
            vmin = []

            for c in wc:
                vmax.append(np.max(np.log10(np.abs(c))))
                vmin.append(np.min(np.log10(np.abs(c))))
            vmax = max(vmax)
            vmin = min(vmin)
            i = 0
            for c in wc:
                out = axs[0, 1].imshow(np.log10(np.abs(c.reshape(1, -1))), extent=[0, max(t), i - 0.5, i + 0.5],
                                       vmin=vmin, vmax=vmax, cmap='rainbow', aspect='auto', interpolation='nearest',
                                       origin='lower')
                i += 1

            axs[0, 1].set_ylim([-.5, i - .5])
            fig.colorbar(out, ax=axs[0, 1])

            n = 0
            for j in range(2):
                for i in range(1, 7):
                    if n < len(wc):
                        xt = inverse_wavelet_transform_j(wc[n], tr.stats.npts, 1 / tr.stats.delta, n)
                        axs[i, j].plot(t, xt, linewidth=.01, color="gray")
                        # axs[i, j].set_ylim(ylim)
                        n += 1
                    else:
                        pass
            pdf.savefig(fig, dpi=300)
            plt.close()


def plot_data_vs_synthetics(args, name, st, sy, spp=8):
    with PdfPages(os.path.join(args.output, 'output', name)) as pdf:
        for k in range(0, len(st), 3 * spp):
            fig, axs = plt.subplots(spp, 3, sharex='col', sharey='row')
            # j = np.log2(st[0].stats.npts)
            # T = st[0].stats.npts * st[0].stats.delta
            # plt.suptitle("scale: j = %d , Frequency: f = [%.4f, %.4f]" % (j, 2**j/3/T, 2**(j+2)/3/T))

            for n in range(spp):
                for i in range(3):
                    axs[n, i].spines['right'].set_visible(False)
                    axs[n, i].spines['left'].set_visible(False)
                    axs[n, i].spines['top'].set_visible(False)
                    axs[n, i].spines['bottom'].set_visible(False)
                    axs[n, i].set_yticks([])
                    axs[n, i].set_xticks([])
                    axs[n, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                          labelsize=6)
                    if n == spp - 1:
                        axs[n, i].spines['bottom'].set_visible(True)
                        axs[n, i].tick_params(axis='x', which='both', bottom=True, labelsize=6)
                        axs[n, i].xaxis.set_major_locator(MultipleLocator(100))
                        axs[n, i].xaxis.set_minor_locator(MultipleLocator(50))
                    else:
                        axs[n, i].tick_params(labelsize=6)

            for i in range(spp - 1, -1, -1):
                for j in range(3):
                    if i == 0:
                        axs[0, j].set_title(st[(spp - 1 - (spp - 1)) * 3 + j + k].stats.channel[-1])
                    if (spp - 1 - i) * 3 + j + k < len(st):
                        t = np.arange(st[(spp - 1 - i) * 3 + j + k].stats.npts) * st[
                            (spp - 1 - i) * 3 + j + k].stats.delta
                        axs[i, j].plot(t, st[(spp - 1 - i) * 3 + j + k].data / np.max(
                            np.abs(st[(spp - 1 - i) * 3 + j + k].data)), linewidth=.15, color="blue")
                        axs[i, j].plot(t, sy[(spp - 1 - i) * 3 + j + k].data / np.max(
                            np.abs(sy[(spp - 1 - i) * 3 + j + k].data)), linewidth=.05, color="red")#, linestyle='dotted')
                        plt.text(0.9, 0.9, "max = %.2e" % np.max(np.abs(st[(spp - 1 - i) * 3 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=4, color='blue')
                        plt.text(0.9, 0.9, "\n\nmax = %.2e" % np.max(np.abs(sy[(spp - 1 - i) * 3 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=4, color='red')
                        # print(st[(9 - i) * 3 + j + k].stats.npts)
                        if j == 0:
                            plt.text(0.1, 0.9, "dist = %.1f\naz = %.1f" % (
                                st[(spp - 1 - i) * 3 + j + k].stats.dist, st[(spp - 1 - i) * 3 + j + k].stats.az),
                                     ha='center',
                                     va='center', transform=axs[i, j].transAxes,
                                     fontsize=4)

                            axs[i, j].set_ylabel(
                                st[(spp - 1 - i) * 3 + j + k].stats.network + "." + st[
                                    (spp - 1 - i) * 3 + j + k].stats.station,
                                fontsize=6, rotation=0, labelpad=30)

            pdf.savefig(fig, dpi=300)
            plt.close()


def plot_data_vs_synthetics_wc(args, name, st, sy, spp=8):
    with PdfPages(os.path.join(args.output, 'output', name)) as pdf:
        for k in range(0, len(st), spp * 3):
            fig, axs = plt.subplots(spp, 3, sharex='col', sharey='row')
            # j = np.log2(st[0].stats.npts)
            # T = st[0].stats.npts * st[0].stats.delta
            # plt.suptitle("j = %d , f = [%.4f, %.4f]" % (j, 2**j/3/T, 2**(j+2)/3/T))

            for n in range(spp):
                for i in range(3):
                    axs[n, i].spines['right'].set_visible(False)
                    axs[n, i].spines['left'].set_visible(True)
                    axs[n, i].spines['top'].set_visible(False)
                    axs[n, i].spines['bottom'].set_visible(True)
                    axs[n, i].set_yticks([])
                    axs[n, i].set_xticks([])
                    axs[n, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                          color='gray', labelsize=6)
                    axs[n, i].spines['bottom'].set_position('center')
                    axs[n, i].spines['bottom'].set_color('gray')
                    axs[n, i].spines['left'].set_color('gray')

                    if n == spp - 1:
                        axs[n, i].spines['bottom'].set_visible(True)
                        axs[n, i].tick_params(axis='x', which='both', bottom=True, labelsize=4)
                        axs[n, i].xaxis.set_major_locator(MultipleLocator(10))
                        axs[n, i].xaxis.set_minor_locator(MultipleLocator(1))
                    else:
                        axs[n, i].tick_params(labelsize=4)

            for i in range(spp - 1, -1, -1):
                for j in range(3):
                    if i == 0:
                        axs[0, j].set_title(st[(spp - 1 - (spp - 1)) * 3 + j + k].stats.channel[-1])
                    if (spp - 1 - i) * 3 + j + k < len(st):
                        n = np.arange(st[(spp - 1 - i) * 3 + j + k].stats.npts)
                        axs[i, j].scatter(n, st[(spp - 1 - i) * 3 + j + k].data, marker="o", color="blue", s=6,
                                          facecolors='none', linewidths=.3)
                        axs[i, j].scatter(n, sy[(spp - 1 - i) * 3 + j + k].data, marker="x", color="red", s=6,
                                          linewidths=.3)
                        plt.text(0.9, 0.9, "max = %.2e" % np.max(np.abs(st[(spp - 1 - i) * 3 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=4, color='gray')
                        plt.text(0.9, 0.9, "\n\nmax = %.2e" % np.max(np.abs(sy[(spp - 1 - i) * 3 + j + k].data)),
                                 ha='center',
                                 va='center', transform=axs[i, j].transAxes,
                                 fontsize=4, color='red')
                        # print(st[(9 - i) * 3 + j + k].stats.npts)
                        if j == 0:
                            plt.text(0.1, 0.9, "dist = %.1f\naz = %.1f" % (
                                st[(spp - 1 - i) * 3 + j + k].stats.dist, st[(spp - 1 - i) * 3 + j + k].stats.az),
                                     ha='center',
                                     va='center', transform=axs[i, j].transAxes,
                                     fontsize=4)

                            axs[i, j].set_ylabel(
                                st[(spp - 1 - i) * 3 + j + k].stats.network + "." + st[
                                    (spp - 1 - i) * 3 + j + k].stats.station,
                                fontsize=6, rotation=0, labelpad=30)

            pdf.savefig(fig, dpi=300)
            plt.close()


def create_map(args, name, st):
    s_lats = []
    s_lons = []
    os.path.join(args.output, 'output', name)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    for tr in st:
        if tr.stats["channel"] == "Z":
            s_lons.append(tr.stats.stlo)
            s_lats.append(tr.stats.stla)
    s_lons.append(tr.stats.evlo)
    s_lats.append(tr.stats.evla)
    max_lo = np.max(s_lons)
    min_lo = np.min(s_lons)
    max_la = np.max(s_lats)
    min_la = np.min(s_lats)
    map = Basemap(projection="cyl", llcrnrlat=min_la - .5, llcrnrlon=min_lo - .5, urcrnrlat=max_la + .5,
                  urcrnrlon=max_lo + .5)
    map.bluemarble()
        # map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=1500, verbose=True)
    map.drawparallels(np.arange(np.floor(min_la - .5), np.ceil(max_la - .5)+1, 1.), labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(np.floor(min_lo - .5), np.ceil(max_lo - .5)+1, 1.), labels=[0, 0, 0, 1])
    x, y = map(s_lons[0:-1], s_lats[0:-1])
    map.scatter(x, y, s=80, marker="^", c='k')
    x,y =map(s_lons[-1], s_lats[-1])
    map.scatter(x, y, s=80, marker="*", c='y')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, 'output', name))
    plt.close()
