import numpy as np
import matplotlib.pyplot as plt


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
    N = 2 ** 13
    fs = 100
    df = fs / N
    f = np.arange(-N / 2, N / 2) * df
    w = 2 * np.pi * f

    s = sy(w)

    plt.figure(1)
    line1, = plt.plot(w, np.real(s), label='real part')
    line2, = plt.plot(w, np.imag(s), label='imaginary part')
    line3, = plt.plot(w, np.abs(s), label='amplitude')
    plt.xlim([-10, 10])
    plt.xlabel(r"$\omega$")
    plt.ylabel(r"Amplitude")
    plt.legend()
    plt.show()

    st = np.fft.ifft(np.fft.ifftshift(s), norm='backward')
    st = np.fft.fftshift(np.real(st))
    t = np.arange(-N / 2, N / 2) / fs

    plt.figure(2)
    plt.plot(t, st * fs)
    plt.ylim([-1.5, 1.5])
    plt.xlim([-5, 5])
    plt.show()

    N = 1024
    fs = 10
    t = np.arange(N) / fs
    xi = np.sin(2 * np.pi * 1 * t) + np.sin(2 * np.pi * .1 * t) + np.sin(2 * np.pi * .01 * t)
    a_jk = wavelet_transform(xi, N, fs)

    xt = inverse_wavelet_transform(a_jk, fs)

    plt.figure(3)
    plt.subplot(211)
    plt.plot(t, xi)
    plt.subplot(212)
    plt.plot(t, xi)
    plt.plot(t, xt)
    plt.show()
    print(np.sum(xt ** 2 - xi ** 2) / np.sum(xi ** 2))

    plt.figure(4)
    max_j = int(np.log2(N))
    for i in range(max_j):
        print(i)
        xt2 = inverse_wavelet_transform_j(a_jk[i], N, fs, i)
        xt = inverse_wavelet_transform_2(a_jk, fs, i, i + 1)
        plt.subplot(max_j, 1, i + 1)
        plt.plot(t, xt)
        plt.plot(t, xt2)
        plt.ylim([-2, 2])
    plt.show()
