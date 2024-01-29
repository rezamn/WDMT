import numpy as np
import scipy
import matplotlib.pyplot as plt
import os

from matplotlib.ticker import MultipleLocator


def h(omega):
    return np.where(omega > 0, np.exp(-1 / omega ** 2), 0)


def g(omega):
    return h(4 * np.pi / 3 - omega) / (h(omega - 2 * np.pi / 3) + h(4 * np.pi / 3 - omega))


def phi(omega):
    return np.sqrt(g(omega) * g(-omega))


def sy(omega):
    return np.exp(-1j * omega / 2) * np.sqrt(phi(omega / 2) * phi(omega / 2) - phi(omega) * phi(omega))


def next_pow2(n):
    return np.uint32(np.ceil(np.log2(n)))


def wavelet_transform_j(input_signal, n, sampling_frequency, j):
    n = next_pow2(n)
    nn = 2 ** n
    aas = np.fft.fft(input_signal, nn, norm='forward')  # norm='forward' means scaled with 1/NN
    k_max = 2 ** j
    ajk = []
    for si in np.arange(k_max):
        ajk.append(aas[si] * np.conj(sy(2 * np.pi * si / k_max)) + aas[si + 2 ** j] * np.conj(
            sy(2 * np.pi + 2 * np.pi * si / k_max)))
    alpha = 2 * np.real(np.fft.ifft(np.array(ajk), norm='forward')) * np.sqrt(nn / sampling_frequency / k_max)
    return np.array(alpha)


def wavelet_transform(input_signal, n, sampling_frequency):
    n = next_pow2(n)
    nn = 2 ** n
    alpha = []
    # aas = np.fft.fft(input_signal, nn, norm='forward')  # norm='forward' means scaled with 1/NN
    for j in np.arange(n):
        alpha.append(wavelet_transform_j(input_signal, nn, sampling_frequency, j))
        # k_max = 2 ** j ajk = [] for si in np.arange(k_max): ajk.append(aas[si] * np.conj(sy(2 * np.pi * si /
        # k_max)) + aas[si + 2 ** j] * np.conj( sy(2 * np.pi + 2 * np.pi * si / k_max))) alpha.append(2 * np.real(
        # np.fft.ifft(np.array(ajk), norm='forward')) * np.sqrt(nn / sampling_frequency / k_max))
    return np.array(alpha, dtype=object)


def inverse_wavelet_transform_j(alpha, n, sampling_frequency, j):
    n = next_pow2(n)
    nn = 2 ** n
    delta_frequency = sampling_frequency / nn
    frequency = np.arange(-nn / 2, nn / 2) * delta_frequency
    omega = 2 * np.pi * frequency
    sw = 0
    # for j in np.arange(np.size(alpha)):
    for k in np.arange(2 ** j):
        a = (2 ** j) * sampling_frequency / nn
        b = k
        sw = sw + np.sqrt(a) * alpha[k] * (1 / a) * (np.exp(-1j * omega * b / a) * sy(omega / a))

    return np.real(np.fft.ifft(np.fft.ifftshift(sw), norm='backward')) * sampling_frequency


def inverse_wavelet_transform(alpha, sampling_frequency):
    n = 2 ** np.size(alpha)
    # delta_frequency = sampling_frequency / n
    # frequency = np.arange(-n / 2, n / 2) * delta_frequency
    # omega = 2 * np.pi * frequency
    signal = 0
    for j in np.arange(np.size(alpha)):
        signal += inverse_wavelet_transform_j(alpha[j], n, sampling_frequency, j)
        # for k in np.arange(2 ** j):
        #     a = (2 ** j) * sampling_frequency / n
        #     b = k
        #     sw += np.sqrt(a) * alpha[j][k] * (1 / a) * (np.exp(-1j * omega * b / a) * sy(omega / a))

    return signal


def inverse_wavelet_transform_2(alpha, sampling_frequency, j_min=0, j_max=0):
    n = 2 ** np.size(alpha)
    delta_frequency = sampling_frequency / n
    frequency = np.arange(-n / 2, n / 2) * delta_frequency
    omega = 2 * np.pi * frequency
    sw = 0
    if j_max == 0 or j_max > 2 ** np.size(alpha):
        j_max = 2 ** np.size(alpha)
    for j in np.arange(j_min, j_max):
        for k in np.arange(2 ** j):
            a = (2 ** j) * sampling_frequency / n
            b = k
            sw += np.sqrt(a) * alpha[j][k] * (1 / a) * (np.exp(-1j * omega * b / a) * sy(omega / a))

    return np.real(np.fft.ifft(np.fft.ifftshift(sw), norm='backward')) * sampling_frequency


if __name__ == '__main__':
    if not os.path.exists("output"):
        os.mkdir("output")

    N = 2 ** 13
    fs = 100
    df = fs / N
    f = np.arange(-N / 2, N / 2) * df
    w = 2 * np.pi * f

    s = sy(w)

    fig = plt.figure(figsize=(5, 6))

    plt.subplot(2, 1, 1)
    line1, = plt.plot(w, np.real(s), label='real part')
    line2, = plt.plot(w, np.imag(s), label='imaginary part')
    line3, = plt.plot(w, np.abs(s), label='amplitude')

    ax = plt.gca()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')

    plt.xlim([-10, 10])
    plt.ylim([-1.5, 1.5])
    plt.xlabel(r"$\omega$", loc="right")
    plt.xticks(np.linspace(-3*np.pi,3*np.pi,7), ['-3π', '-2π', '-π', '0', 'π', '2π', '3π'])
    plt.ylabel(r"Amplitude", loc="top", rotation=0)
    plt.text(-9, 1.25, "(a)")
    plt.legend(loc="upper right", fontsize=6)

    st = np.fft.ifft(np.fft.ifftshift(s), norm='backward')
    st = np.fft.fftshift(np.real(st))
    t = np.arange(-N / 2, N / 2) / fs

    plt.subplot(2, 1, 2)
    plt.plot(t, st * fs)

    ax = plt.gca()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')

    plt.ylim([-1.5, 1.5])
    plt.xlim([-5, 5])
    plt.xlabel(r"t", loc="right")
    plt.ylabel(r"Amplitude", loc="top", rotation=0)
    plt.text(-4.5, 1.25, "(b)")
    plt.tight_layout()
    plt.savefig(".\\output\\fig1.eps")
    plt.savefig(".\\output\\fig1.pdf")
    plt.show()

    N = 1024
    fs = 10
    T = N / fs
    t = np.arange(N) / fs
    xi = np.sin(2 * np.pi * 1 * t) + np.sin(2 * np.pi * .1 * t) + np.sin(2 * np.pi * .01 * t)
    a_jk = wavelet_transform(xi, N, fs)

    xt = inverse_wavelet_transform(a_jk, fs)

    # plt.figure(3)
    # plt.subplot(211)
    # plt.plot(t, xi)
    # plt.subplot(212)
    # plt.plot(t, xi)
    # plt.plot(t, xt)
    # plt.show()
    # print(np.sum(xt ** 2 - xi ** 2) / np.sum(xi ** 2))

    max_j = np.uint(np.log2(N))
    fig, axs = plt.subplots(max_j + 1, 2, sharex='col', sharey='row', figsize=(12, 8))

    for i in np.arange(max_j):
        xt2 = inverse_wavelet_transform_j(a_jk[i], N, fs, i)
        # xt = inverse_wavelet_transform_2(a_jk, fs, i, i + 1)
        # plt.subplot(max_j, 1, i + 1)
        # axs[i, 1].plot(t, xt)
        axs[i + 1][0].plot(t, xt2, linewidth=.5, color="blue")
        axs[i + 1][0].set_ylim(-2, 2)
        if i != max_j - 1:
            axs[i + 1][0].spines['right'].set_visible(False)
            axs[i + 1][0].spines['left'].set_visible(True)
            axs[i + 1][0].spines['top'].set_visible(False)
            axs[i + 1][0].spines['bottom'].set_visible(False)
            axs[i + 1][0].yaxis.set_major_locator(MultipleLocator(2))
            axs[i + 1][0].yaxis.set_minor_locator(MultipleLocator(1))
            # axs[i].set_yticks([])
            axs[i + 1][0].set_xticks([])
            axs[i + 1][0].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=True,
                                   labelsize=6)
        else:
            axs[i + 1][0].spines['bottom'].set_visible(True)
            axs[i + 1][0].spines['right'].set_visible(False)
            axs[i + 1][0].spines['top'].set_visible(False)
            axs[i + 1][0].tick_params(axis='both', which='both', bottom=True, labelsize=6)
            axs[i + 1][0].xaxis.set_major_locator(MultipleLocator(20))
            axs[i + 1][0].xaxis.set_minor_locator(MultipleLocator(10))
            axs[i + 1][0].yaxis.set_major_locator(MultipleLocator(2))
            axs[i + 1][0].yaxis.set_minor_locator(MultipleLocator(1))
            axs[i + 1][0].set_xlabel("Time (sec)")
        axs[i + 1][0].set_ylabel("Level: j= %d\nfreq: f= [%.2f %.2f]" % (i, 2 ** i / 3 / T, 2 ** i * 4 / 3 / T),
                              fontsize=8, rotation=0, labelpad=40)


    axs[0][0].plot(t, xi, linewidth=.5, color="red")
    axs[0][0].spines['right'].set_visible(False)
    axs[0][0].spines['left'].set_visible(True)
    axs[0][0].spines['top'].set_visible(False)
    axs[0][0].spines['bottom'].set_visible(False)
    axs[0][0].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=True,
                       labelsize=6)
    axs[0][0].set_ylim(-3, 3)
    axs[0][0].yaxis.set_major_locator(MultipleLocator(3))
    axs[0][0].yaxis.set_minor_locator(MultipleLocator(1))
    axs[0][0].set_ylabel("Original signal", fontsize=8, rotation=0, labelpad=30)
    axs[0][0].set_title("(a)", fontsize=10)


    plt.tight_layout()
    # plt.savefig(".\\output\\fig2.eps")
    # plt.show()

    ###############
    N = 1024
    fs = 10
    T = N / fs
    t = np.arange(N) / fs
    xi = np.sin(2 * np.pi * 1 * t) + np.sin(2 * np.pi * .1 * t) + np.sin(2 * np.pi * .01 * t)
    noise = np.random.randn(xi.size) * .15 * np.max(np.abs(xi))
    xi += noise
    win = scipy.signal.windows.tukey(N, alpha=.1)
    xi = scipy.signal.detrend(xi)*win
    a_jk = wavelet_transform(xi, N, fs)
    xt = inverse_wavelet_transform(a_jk, fs)

    # plt.figure(3)
    # plt.subplot(211)
    # plt.plot(t, xi)
    # plt.subplot(212)
    # plt.plot(t, xi)
    # plt.plot(t, xt)
    # plt.show()
    # print(np.sum(xt ** 2 - xi ** 2) / np.sum(xi ** 2))

    max_j = np.uint(np.log2(N))
    # fig, axs = plt.subplots(max_j + 1, 1, sharex='col', sharey='row', figsize=(6, 8))

    for i in np.arange(max_j):
        xt2 = inverse_wavelet_transform_j(a_jk[i], N, fs, i)
        # xt = inverse_wavelet_transform_2(a_jk, fs, i, i + 1)
        # plt.subplot(max_j, 1, i + 1)
        # axs[i, 1].plot(t, xt)
        axs[i + 1][1].plot(t, xt2, linewidth=.5, color="blue")
        axs[i + 1][1].set_ylim(-2, 2)
        if i != max_j - 1:
            axs[i + 1][1].spines['right'].set_visible(False)
            axs[i + 1][1].spines['left'].set_visible(False)
            axs[i + 1][1].spines['top'].set_visible(False)
            axs[i + 1][1].spines['bottom'].set_visible(False)
            # axs[i + 1][1].yaxis.set_major_locator(MultipleLocator(2))
            # axs[i + 1][1].yaxis.set_minor_locator(MultipleLocator(1))
            # axs[i].set_yticks([])
            axs[i + 1][1].set_xticks([])
            axs[i + 1][1].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                                   labelsize=6)
        else:
            axs[i + 1][1].spines['bottom'].set_visible(True)
            axs[i + 1][1].spines['left'].set_visible(False)
            axs[i + 1][1].spines['right'].set_visible(False)
            axs[i + 1][1].spines['top'].set_visible(False)
            axs[i + 1][1].tick_params(axis='both', which='both', bottom=True, left=False, labelsize=6)
            axs[i + 1][1].xaxis.set_major_locator(MultipleLocator(20))
            axs[i + 1][1].xaxis.set_minor_locator(MultipleLocator(10))
            # axs[i + 1][1].yaxis.set_major_locator(MultipleLocator(2))
            # axs[i + 1][1].yaxis.set_minor_locator(MultipleLocator(1))
            axs[i + 1][1].set_xlabel("Time (sec)")
        # axs[i + 1][1].set_ylabel("scale: j= %d\nfreq: f= [%.2f %.2f]" % (i, 2 ** i / 3 / T, 2 ** i * 4 / 3 / T),
        #                       fontsize=6, rotation=0, labelpad=30)

    axs[0][1].plot(t, xi, linewidth=.5, color="red")
    axs[0][1].spines['right'].set_visible(False)
    axs[0][1].spines['left'].set_visible(False)
    axs[0][1].spines['top'].set_visible(False)
    axs[0][1].spines['bottom'].set_visible(False)
    axs[0][1].tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False,
                       labelsize=6)
    axs[0][1].set_ylim(-3, 3)
    axs[0][1].yaxis.set_major_locator(MultipleLocator(3))
    axs[0][1].yaxis.set_minor_locator(MultipleLocator(1))
    # axs[0][1].set_ylabel("Original signal", fontsize=6, rotation=0, labelpad=30)
    axs[0][1].set_title("(b)", fontsize=10)

    plt.tight_layout()
    plt.savefig(".\\output\\fig2.eps")
    plt.savefig(".\\output\\fig2.pdf")
    plt.show()