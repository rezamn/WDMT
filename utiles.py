import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

from wavelet import wavelet_transform, inverse_wavelet_transform_j
from obspy import read
import os


def plot_stream_wavelets(stream, output_path, output_name):
    print(stream)
    stream = read(stream)
    stream.detrend("linear")
    stream.decimate(2)
    stream.decimate(3)
    # stream.decimate(2)
    stream.decimate(5)

    print(stream)
    output_name = os.path.join(output_path, output_name)
    fig, axs = plt.subplots(4, 3, figsize=(16, 9), gridspec_kw={'height_ratios': [1, 3, 1, 1]})
    ax_list = []
    for  i in range(3):
        for j in range(2):
            ax_list.append(axs[i,j])
    ax_list[0].get_shared_x_axes().join(ax_list[0], *ax_list)
    ax1_list = [axs[3,0], axs[3,1], axs[3,2]]
    ax1_list[0].get_shared_y_axes().join(ax1_list[0], *ax1_list)

    t = np.arange(stream[0].stats.npts) * stream[0].stats.delta
    for i in range(3):
        axs[0, i].plot(t, stream[i].data / np.max(np.abs(stream[i].data)))
        if i != 0:
            axs[0, i].spines['right'].set_visible(False)
            axs[0, i].spines['left'].set_visible(False)
            axs[0, i].spines['top'].set_visible(False)
            axs[0, i].spines['bottom'].set_visible(False)
            axs[0, i].set_yticks([])
            axs[0, i].set_xticks([])
            axs[0, i].set_ylim([-1.1, 1.1])
            axs[0, i].set_title(stream[i].stats.channel[-1], fontsize=16)
            axs[0, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                  labelsize=10)
        else:
            axs[0, i].spines['right'].set_visible(False)
            axs[0, i].spines['left'].set_visible(True)
            axs[0, i].spines['top'].set_visible(False)
            axs[0, i].spines['bottom'].set_visible(False)
            axs[0, i].set_xticks([])
            axs[0, i].set_title(stream[i].stats.channel[-1], fontsize=16)
            axs[0, i].set_ylabel("Amplitude", fontsize=12)
            axs[0, i].set_ylim([-1.1, 1.1])
            axs[0, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=True,
                                  labelsize=10)
    wc = []
    for i in range(3):
        if i != 0:
            axs[1, i].spines['right'].set_visible(False)
            axs[1, i].spines['left'].set_visible(False)
            axs[1, i].spines['top'].set_visible(False)
            axs[1, i].spines['bottom'].set_visible(False)
            axs[1, i].set_yticks([])
            axs[1, i].set_xticks([])
            axs[1, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                  labelsize=10)
        else:
            axs[1, i].spines['right'].set_visible(False)
            axs[1, i].spines['left'].set_visible(True)
            axs[1, i].spines['top'].set_visible(False)
            axs[1, i].spines['bottom'].set_visible(False)
            axs[1, i].set_ylabel("Level (j)", fontsize=12)
            axs[1, i].set_xticks([])
            axs[1, i].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                  labelsize=10)
        wc = wavelet_transform(stream[i].data, stream[i].stats.npts, 1 / stream[i].stats.delta)
        vmax = []
        vmin = []
        for c in wc:
            vmax.append(np.max(np.log10(np.abs(c))))
            vmin.append(np.min(np.log10(np.abs(c))))
        # vmax = max(vmax)
        # vmin = min(vmin)
        vmax = 6
        vmin = -6
        print(vmax, vmin)
        cmap = matplotlib.colormaps.get_cmap('jet')
        j = 0
        for c in wc:
            divider = make_axes_locatable(axs[1, i])
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            im = axs[1, i].imshow(np.log10(np.abs(c.reshape(1, -1))), extent=[0, max(t), j - 0.5, j + 0.5],
                                   vmin=vmin, vmax=vmax, cmap=cmap, aspect='auto', interpolation='none',
                                   origin='lower')
            fig.colorbar(im, cax=cax, orientation='horizontal')
            j += 1

        axs[1, i].hlines(y=6, xmin=np.min(t), xmax=np.max(t), color="black", linestyles="dashed")
        axs[1, i].set_ylim([-.5, j - .5])

        xt = inverse_wavelet_transform_j(wc[6], stream[i].stats.npts, 1 / stream[i].stats.delta, 6)
        axs[2, i].plot(t, xt[:617]/np.max(np.abs(xt[:617])), linewidth=.01, color="gray")
        if i != 0:
            axs[2, i].spines['right'].set_visible(False)
            axs[2, i].spines['left'].set_visible(False)
            axs[2, i].spines['top'].set_visible(False)
            axs[2, i].spines['bottom'].set_visible(True)
            axs[2, i].set_yticks([])
            axs[2, i].set_xlabel("Time (sec)", fontsize=14)
            axs[2, i].set_ylim([-1.1, 1.1])
            axs[2, i].tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
                                  labelsize=10)
        else:
            axs[2, i].spines['right'].set_visible(False)
            axs[2, i].spines['left'].set_visible(True)
            axs[2, i].spines['top'].set_visible(False)
            axs[2, i].spines['bottom'].set_visible(True)
            axs[2, i].set_ylabel("Amplitde",fontsize=12)
            axs[2, i].set_xlabel("Time (sec)", fontsize=14)
            axs[2, i].set_ylim([-1.1, 1.1])
            axs[2, i].tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=True,
                                  labelsize=10)

        nn = np.arange(len(wc[6]))
        axs[3, i].scatter(nn, wc[6], s=10, color="blue")
        if i != 0:
            axs[3, i].spines['right'].set_visible(False)
            axs[3, i].spines['left'].set_visible(False)
            axs[3, i].spines['top'].set_visible(False)
            axs[3, i].spines['bottom'].set_visible(True)
            axs[3, i].set_yticks([])
            axs[3, i].set_xlabel("Translation index", fontsize=14)
            # axs[3, i].set_ylim([-1.1, 1.1])
            axs[3, i].tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=False,
                                  labelsize=10)
        else:
            axs[3, i].spines['right'].set_visible(False)
            axs[3, i].spines['left'].set_visible(True)
            axs[3, i].spines['top'].set_visible(False)
            axs[3, i].spines['bottom'].set_visible(True)
            axs[3, i].set_xlabel("Translation index", fontsize=14)
            axs[3, i].set_ylabel("Wavelet\nvalue", fontsize=12)
            # axs[3, i].set_ylim([-1.1, 1.1])
            axs[3, i].ticklabel_format(axis='y', style='sci', scilimits=(1, 4))
            axs[3, i].tick_params(axis='both', which='both', bottom=True, top=False, right=False, left=True,
                                  labelsize=10)


    fig.tight_layout()
    fig.savefig(output_name + ".eps")
    fig.savefig(output_name + ".pdf")
    plt.close()


if __name__ == "__main__":
    data = "D:/WDMT/EX1/TEST/20141230_041934.9/data/I2.AHRM.*.SAC"
    # data = "D:/WDMT/syn_data/XX.STA03.*.SAC"
    out_path = "D:/WDMT/EX1/TEST/output"
    plot_stream_wavelets(data, out_path, "fig05")
