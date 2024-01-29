import numpy as np
from obspy.signal.cross_correlation import xcorr_max, correlate
from obspy.imaging.beachball import MomentTensor, mt2axes, aux_plane, mt2plane
import copy


def correlate_phase_shift(gr, st):
    # Correlation for Zcor for all and for each component
    # T --> TSS,TDS         | S1 - u1,u2
    # R --> RSS,RDS,RDD     | S2 - u3,u4,u5
    # Z --> ZSS,ZDS,ZDD     | S3 - u6,u7,u8

    # width for time cross correlation
    print("gr", len(gr), "st", len(st))
    wid = int(float(st[0].stats.npts)/10)

    for i in range(int(len(st) / 3)):

        x = np.arange(16.).reshape(8, 2)
        t = np.arange(4.).reshape(2, 2)
        r = np.arange(6.).reshape(3, 2)
        v = np.arange(6.).reshape(3, 2)

        # Tangential
        a, b = xcorr_max(correlate(st[3 * i + 0], gr[10 * i + 0], wid))
        x[0][0] = abs(b)
        x[0][1] = a
        t[0][0] = abs(b)
        t[0][1] = a
        a, b = xcorr_max(correlate(st[3 * i + 0], gr[10 * i + 1], wid))
        x[1][0] = abs(b)
        x[1][1] = a
        t[1][0] = abs(b)
        t[1][1] = a

        # Radial
        a, b = xcorr_max(correlate(st[3 * i + 1], gr[10 * i + 2], wid))
        x[2][0] = abs(b)
        x[2][1] = a
        r[0][0] = abs(b)
        r[0][1] = a
        a, b = xcorr_max(correlate(st[3 * i + 1], gr[10 * i + 3], wid))
        x[3][0] = abs(b)
        x[3][1] = a
        r[1][0] = abs(b)
        r[1][1] = a
        a, b = xcorr_max(correlate(st[3 * i + 1], gr[10 * i + 4], wid))
        x[4][0] = abs(b)
        x[4][1] = a
        r[2][0] = abs(b)
        r[2][1] = a

        # Vertical
        a, b = xcorr_max(correlate(st[3 * i + 2], gr[10 * i + 5], wid))
        x[5][0] = abs(b)
        x[5][1] = a
        v[0][0] = abs(b)
        v[0][1] = a
        a, b = xcorr_max(correlate(st[3 * i + 2], gr[10 * i + 6], wid))
        x[6][0] = abs(b)
        x[6][1] = a
        v[1][0] = abs(b)
        v[1][1] = a
        a, b = xcorr_max(correlate(st[3 * i + 2], gr[10 * i + 7], wid))
        x[7][0] = abs(b)
        x[7][1] = a
        v[2][0] = abs(b)
        v[2][1] = a

        # sort for zcor
        X = np.array(sorted(sorted(x, key=lambda e: e[1]), key=lambda e: e[0]))
        T = np.array(sorted(sorted(t, key=lambda e: e[1]), key=lambda e: e[0]))
        R = np.array(sorted(sorted(r, key=lambda e: e[1]), key=lambda e: e[0]))
        V = np.array(sorted(sorted(v, key=lambda e: e[1]), key=lambda e: e[0]))
        Zco = X[-1][1]
        Tco = T[-1][1]
        Rco = R[-1][1]
        Vco = V[-1][1]

        # Update stats
        for l in range(0, 3):
            st[3 * i + l].stats.Zcor = Zco
            st[3 * i + l].stats.Tcor = Tco
            st[3 * i + l].stats.Rcor = Rco
            st[3 * i + l].stats.Vcor = Vco
        for l in range(0, 10):
            gr[10 * i + l].stats.Zcor = Zco
            gr[10 * i + l].stats.Tcor = Tco
            gr[10 * i + l].stats.Rcor = Rco
            gr[10 * i + l].stats.Vcor = Vco

    return gr, st


def new_correlate_phase_shift(gr, st):
    # Correlation for Zcor for all and for each component
    # T --> TSS,TDS         | S1 - u1,u2
    # R --> RSS,RDS,RDD     | S2 - u3,u4,u5
    # Z --> ZSS,ZDS,ZDD     | S3 - u6,u7,u8
    wid = 'full'
    for i in range(int(len(st) / 3)):

        t = np.zeros(2)
        r = np.zeros(3)
        v = np.zeros(3)
        #TODO: normalize the greens functins in correlation
        #TODO: take all shifts from RTZ greens function
        # Tangential
        t[0] = np.argmax(np.abs(np.correlate(st[3 * i + 0].data, gr[10 * i + 0].data, wid))) - gr[10 * i + 0].stats.npts
        t[1] = np.argmax(np.abs(np.correlate(st[3 * i + 0].data, gr[10 * i + 1].data, wid))) - gr[10 * i + 1].stats.npts

        # Radial
        r[0] = np.argmax(np.abs(np.correlate(st[3 * i + 1].data, gr[10 * i + 2].data, wid))) - gr[10 * i + 2].stats.npts
        r[1] = np.argmax(np.abs(np.correlate(st[3 * i + 1].data, gr[10 * i + 3].data, wid))) - gr[10 * i + 3].stats.npts
        r[2] = np.argmax(np.abs(np.correlate(st[3 * i + 1].data, gr[10 * i + 4].data, wid))) - gr[10 * i + 4].stats.npts

        # Vertical
        v[0] = np.argmax(np.abs(np.correlate(st[3 * i + 2].data, gr[10 * i + 5].data, wid))) - gr[10 * i + 5].stats.npts
        v[1] = np.argmax(np.abs(np.correlate(st[3 * i + 2].data, gr[10 * i + 6].data, wid))) - gr[10 * i + 6].stats.npts
        v[2] = np.argmax(np.abs(np.correlate(st[3 * i + 2].data, gr[10 * i + 7].data, wid))) - gr[10 * i + 7].stats.npts

        # sort for zcor
        T = np.min(t)
        R = np.min(r)
        V = np.min(v)
        X = (T + R + V) / 3
        Zco = V
        Tco = T
        Rco = R
        Vco = X

        # Update stats
        for l in range(0, 3):
            st[3 * i + l].stats.Zcor = Zco
            st[3 * i + l].stats.Tcor = Tco
            st[3 * i + l].stats.Rcor = Rco
            st[3 * i + l].stats.Vcor = Vco
        for l in range(0, 10):
            gr[10 * i + l].stats.Zcor = Zco
            gr[10 * i + l].stats.Tcor = Tco
            gr[10 * i + l].stats.Rcor = Rco
            gr[10 * i + l].stats.Vcor = Vco

    return gr, st


def new_correlate_shift(st, syn):
    wid = 'full'
    for i in range(int(len(st) // 3)):
        # print(i, st[3 * i + 0].stats.station)

        t = np.zeros(1)
        r = np.zeros(1)
        v = np.zeros(1)

        # Tangential
        t[0] = np.argmax(np.correlate(st[3 * i + 0].data, syn[3 * i + 0].data, wid)) - syn[3 * i + 0].stats.npts - 1

        # Radials
        r[0] = np.argmax(np.correlate(st[3 * i + 1].data, syn[3 * i + 1].data, wid)) - syn[3 * i + 1].stats.npts - 1

        # Vertical
        v[0] = np.argmax(np.correlate(st[3 * i + 2].data, syn[3 * i + 2].data, wid)) - syn[3 * i + 2].stats.npts - 1

        # sort for zcor
        T = t[0]
        R = r[0]
        V = v[0]
        X = (T + R + V) / 3
        # print(st[3*i+0].stats.station, T, R, V)
        Zco = V
        Tco = T
        Rco = R
        Vco = X

        # Update stats
        for l in range(3):
            st[3 * i + l].stats.Zcor2 = Zco
            st[3 * i + l].stats.Tcor2 = Tco
            st[3 * i + l].stats.Rcor2 = Rco
            st[3 * i + l].stats.Vcor2 = Vco

    return st, syn


def align_phase(st, args):
    rvel = float(args.rvel)
    pre = float(args.pre)
    if args.zcor == "iso":

        for i in range(len(st)):

            ph = int(st[i].stats.Vcor)

            if ph >= 1:
                st[i].data = np.roll(st[i].data, -ph)
                for k in range(1, ph):
                    st[i].data[-k] = 0

            elif ph <= 1:
                st[i].data = np.roll(st[i].data, -ph)
                for k in range(0, ph):
                    st[i].data[k] = 0

            else:
                pass

    else:
        for i in range(len(st)):

            if st[i].stats.channel[-1] == "Z":
                ph = int(st[i].stats.Zcor)
            elif st[i].stats.channel[-1] == "R":
                ph = int(st[i].stats.Rcor)
            else:
                ph = int(st[i].stats.Tcor)

            if ph >= 1:
                st[i].data = np.roll(st[i].data, -ph)
                for k in range(1, ph):
                    st[i].data[-k] = 0

            elif ph <= 1:
                st[i].data = np.roll(st[i].data, -ph)
                for k in range(0, ph):
                    st[i].data[k] = 0

            else:
                pass
    return st

def align_phase2(st, args):

    for i in range(len(st)):

        if st[i].stats.channel[-1] == "Z":
            ph = int(st[i].stats.Zcor2)
        elif st[i].stats.channel[-1] == "R":
            ph = int(st[i].stats.Rcor2)
        elif st[i].stats.channel[-1] == "T":
            ph = int(st[i].stats.Tcor2)

        if ph >= 1:
            st[i].data = np.roll(st[i].data, -ph)
            for k in range(1, ph):
                st[i].data[-k] = 0

        elif ph <= 1:
            st[i].data = np.roll(st[i].data, -ph)
            for k in range(0, ph):
                st[i].data[k] = 0

        else:
            pass

    return st


def trim_streams(gr, st):
    npts = gr[0].stats.npts

    for i in range(len(st)):
        if st[i].stats.npts >= npts:
            st[i].data = st[i].data[:npts]
        else:
            st[i].data = st[i].data.resize[npts]

    return gr, st


def set_station_weights(gr, st):
    # find closer distance
    min_dist = 1000000
    for i in range(len(st)):
        if st[i].stats.dist <= min_dist:
            min_dist = st[i].stats.dist

    # set station weights
    for i in range(len(st)):
        st[i].stats.W = st[i].stats.dist / min_dist
    for i in range(len(gr)):
        gr[i].stats.W = gr[i].stats.dist / min_dist

    return gr, st


def set_GTG(args, data, greens):
    if args.iso == "0":
        iso_flag = 5
    else:
        iso_flag = 6

    # set W matrix
    W = np.zeros(len(data) * data[0].stats.npts)
    n = 0
    for i in range(len(data)):
        for j in range(data[i].stats.npts):
            W[n] = data[i].stats.W
            n += 1

    # set GTG matrix
    AIV = np.zeros(shape=(iso_flag, iso_flag))

    # set GTd matrix
    B = np.zeros(shape=(iso_flag, 2))

    # set AJ matrix
    # -- Initialize
    # number of elements for each row: 3* Nrstaz * npts
    Np = data[0].stats.npts
    nsta = int(len(data) / 3)
    nrElementsInRow = int(3 * Np * nsta)
    AJ = np.zeros(shape=(iso_flag, nrElementsInRow))
    #
    # -- allocate gr and st
    cnt1 = cnt2 = cnt3 = 0
    for i in range(0, nsta):
        Np = data[i].stats.npts
        Z = 0
        cnt1 = cnt2 = cnt3
        #       cnt2 = cnt3
        cnt2 += Np
        cnt3 += 2 * Np
        #       print i,Np,Z,cnt1,cnt2,cnt3
        for j in range(0, Np):

            alpha = data[i * 3].stats.az * np.pi / 180.

            # Mxx term
            AJ[0][cnt1] = 0.5 * np.sin(2 * alpha) * greens[i * 10 + 0].data[j]
            if (iso_flag == 6):
                AJ[0][cnt2] = (1. / 6.) * (-1.0) * greens[i * 10 + 4].data[j] - 0.5 * np.cos(2. * alpha) * \
                              greens[i * 10 + 2].data[j] + (1. / 3.) * (-1.0) * greens[i * 10 + 8].data[j]
                AJ[0][cnt3] = (1. / 6.) * (-1.0) * greens[i * 10 + 7].data[j] - 0.5 * np.cos(2. * alpha) * (-1.0) * \
                              greens[i * 10 + 5].data[j] + (1. / 3.) * (-1.0) * greens[i * 10 + 9].data[j]
            else:
                AJ[0][cnt2] = (1. / 2.) * greens[i * 10 + 4].data[j] - 0.5 * np.cos(2. * alpha) * \
                              greens[i * 10 + 2].data[j]
                AJ[0][cnt3] = (1. / 2.) * (-1.0) * greens[i * 10 + 7].data[j] - 0.5 * np.cos(2. * alpha) * (-1.0) * \
                              greens[i * 10 + 5].data[j]

            # Myy term
            AJ[1][cnt1] = -1.0 * 0.5 * np.sin(2. * alpha) * greens[i * 10 + 0].data[j]
            if (iso_flag == 6):
                AJ[1][cnt2] = (1. / 6.) * greens[i * 10 + 4].data[j] + 0.5 * np.cos(2. * alpha) * \
                              greens[i * 10 + 2].data[
                                  j] + (1. / 3.) * (-1.0) * greens[i * 10 + 8].data[j]
                AJ[1][cnt3] = (1. / 6.) * (-1.0) * greens[i * 10 + 7].data[j] + 0.5 * np.cos(2. * alpha) * (-1.0) * \
                              greens[i * 10 + 5].data[j] + (1. / 3.) * (-1.0) * greens[i * 10 + 9].data[j]
            else:
                AJ[1][cnt2] = (1. / 2.) * greens[i * 10 + 4].data[j] + 0.5 * np.cos(2. * alpha) * \
                              greens[i * 10 + 2].data[j]
                AJ[1][cnt3] = (1. / 2.) * (-1.0) * greens[i * 10 + 7].data[j] + 0.5 * np.cos(2. * alpha) * (-1.0) * \
                              greens[i * 10 + 5].data[j]

            # Mxy term
            AJ[2][cnt1] = (-1.0) * np.cos(2. * alpha) * greens[i * 10 + 0].data[j]
            AJ[2][cnt2] = (-1.0) * np.sin(2. * alpha) * greens[i * 10 + 2].data[j]
            AJ[2][cnt3] = (-1.0) * np.sin(2. * alpha) * (-1.0) * greens[i * 10 + 5].data[j]
            #           print("----   %d %d %g" % (i,j,gr[i*10+5].data[j]))

            # Mxz term
            AJ[3][cnt1] = (-1.0) * np.sin(alpha) * greens[i * 10 + 1].data[j]
            AJ[3][cnt2] = np.cos(alpha) * greens[i * 10 + 3].data[j]
            AJ[3][cnt3] = np.cos(alpha) * (-1.0) * greens[i * 10 + 6].data[j]

            # Myz term*/
            AJ[4][cnt1] = np.cos(alpha) * greens[i * 10 + 1].data[j]
            AJ[4][cnt2] = np.sin(alpha) * greens[i * 10 + 3].data[j]
            AJ[4][cnt3] = np.sin(alpha) * (-1.0) * greens[i * 10 + 6].data[j]

            # Mzz term*/
            if (iso_flag == 6):
                AJ[5][cnt1] = 0.0
                AJ[5][cnt2] = 1 / 3 * (-1.0) * greens[i * 10 + 8].data[j] - 1 / 3 * (1.0) * greens[i * 10 + 4].data[j]
                AJ[5][cnt3] = 1 / 3 * (-1.0) * greens[i * 10 + 9].data[j] - 1 / 3 * (-1.0) * greens[i * 10 + 7].data[j]

            # increment counters
            cnt1 += 1
            cnt2 += 1
            cnt3 += 1

    # Compute GTG
    for i in range(0, iso_flag):
        for j in range(0, iso_flag):
            for k in range(0, cnt3):
                AIV[i][j] += AJ[i][k] * AJ[j][k] * W[k]

    return AIV, B, AJ, W


def righthand_side(args, st, B, AJ, W):
    if args.iso == "0":
        iso_flag = 5
    else:
        iso_flag = 6

    cnt1 = cnt2 = cnt3 = 0
    tmp = np.zeros(3 * len(st) * st[0].stats.npts)
    for i in range(len(st) // 3):
        Np = st[i].stats.npts
        Z = 0
        cnt1 = cnt2 = cnt3
        #       cnt2 = cnt3
        cnt2 += Np
        cnt3 += 2 * Np
        #       for j in range(0,Np):
        for j in range(Z, Np + Z):
            if j < len(st[i * 3 + 0].data):
                # print(i, j, cnt1, cnt2, cnt3)
                tmp[cnt1] = st[i * 3 + 0].data[j]
                tmp[cnt2] = st[i * 3 + 1].data[j]
                tmp[cnt3] = st[i * 3 + 2].data[j]
            #             print("%d %d %g %g %g" % (i,j,st[i*3+0].data[j], st[i*3+1].data[j], st[i*3+2].data[j]))
            else:
                tmp[cnt1] = 0.0
                tmp[cnt2] = 0.0
                tmp[cnt3] = 0.0

            cnt1 += 1
            cnt2 += 1
            cnt3 += 1

    # To build B must realign j to Z
    nn = 0
    for i in range(0, iso_flag):
        for j in range(0 + 0, cnt3 + 0):
            B[i][0] += AJ[i][j - 0] * tmp[j] * W[j]
            nn += 1

    return B


def to_rtf(moment_vector, scale):
    rtf = np.zeros(shape=(3, 3))
    rtf[0][0] = (1.0 * scale * moment_vector[0])
    rtf[0][1] = (-1.0 * scale * moment_vector[2])
    rtf[0][2] = (1.0 * scale * moment_vector[3])
    rtf[1][0] = (-1.0 * scale * moment_vector[2])
    rtf[1][1] = (1.0 * scale * moment_vector[1])
    rtf[1][2] = (-1.0 * scale * moment_vector[4])
    rtf[2][0] = (1.0 * scale * moment_vector[3])
    rtf[2][1] = (-1.0 * scale * moment_vector[4])
    rtf[2][2] = (1.0 * scale * moment_vector[5])

    return rtf


def to_xyz(moment_vector, scale):
    xyz = np.zeros(shape=(3, 3))
    xyz[0][0] = (1.0 * scale * moment_vector[0])
    xyz[0][1] = (1.0 * scale * moment_vector[2])
    xyz[0][2] = (1.0 * scale * moment_vector[3])
    xyz[1][0] = (1.0 * scale * moment_vector[2])
    xyz[1][1] = (1.0 * scale * moment_vector[1])
    xyz[1][2] = (1.0 * scale * moment_vector[4])
    xyz[2][0] = (1.0 * scale * moment_vector[3])
    xyz[2][1] = (1.0 * scale * moment_vector[4])
    xyz[2][2] = (1.0 * scale * moment_vector[5])

    return xyz


def extract_parameters_mt(MTx, scale):
    # Iso Mo
    MoIso = (MTx[0][0] + MTx[1][1] + MTx[2][2]) / 3

    c = MomentTensor(MTx[2][2], MTx[0][0], MTx[1][1], MTx[0][2], -MTx[1][2], -MTx[0][1], scale)

    # Principal axes
    (T, N, P) = mt2axes(c)

    # Nodal planes
    np0 = mt2plane(c)
    np2 = aux_plane(np0.strike, np0.dip, np0.rake)
    # Convention rake: up-down
    if np0.rake > 180:
        np0.rake = np0.rake - 360
    if np0.rake < -180:
        np0.rake = np0.rake + 360

    np1 = [np0.strike, np0.dip, np0.rake]

    # Compute Eigenvectors and Eigenvalues
    # Seismic Moment and Moment Magnitude
    (EigVal, EigVec) = np.linalg.eig(MTx)
    b = copy.deepcopy(EigVal)
    b.sort()
    Mo = (abs(b[0]) + abs(b[2])) / 2.
    Mw = np.log10(Mo) / 1.5 - 10.73

    # Compute double-couple, CLVD & iso
    d = copy.deepcopy(EigVal)
    d[0] = abs(d[0])
    d[1] = abs(d[1])
    d[2] = abs(d[2])
    d.sort()
    eps = abs(d[0]) / abs(d[2])
    pcdc = 100.0 * (1.0 - 2.0 * eps)
    pcclvd = 200.0 * eps
    pcdc = pcdc / 100.0
    pcclvd = pcclvd / 100.0
    pciso = abs(MoIso) / Mo
    pcsum = pcdc + pcclvd + pciso
    pcdc = 100.0 * pcdc / pcsum
    pcclvd = 100.0 * pcclvd / pcsum
    pciso = 100.0 * pciso / pcsum

    Pdc = pcdc
    Pclvd = pcclvd

    return Mo, Mw, Pdc, Pclvd, pciso, EigVal, EigVec, T, N, P, np1, np2


def check_fit1(st, sy):
    weight_sum = 0.0
    variance = 0.0

    for i in range(len(st)):
        st[i].data = st[i].data / np.max(np.abs(st[i].data))
        sy[i].data = sy[i].data / np.max(np.abs(sy[i].data))

    for i in range(len(st)):
        Var = sum(st[i].data*sy[i].data)/np.sqrt(sum(st[i].data**2)*sum(sy[i].data**2))
        st[i].stats.VR = Var
        weight_sum += st[i].stats.W
        variance += st[i].stats.W * Var

    variance /= weight_sum
    variance *= 100

    if variance < 20.0:
        quality = 0
    elif 20.0 <= variance < 40.0:
        quality = 1
    elif 40.0 <= variance < 60.0:
        quality = 2
    elif 60.0 <= variance < 80.0:
        quality = 3
    else:
        quality = 4

    return st, sy, variance, quality



def check_fit(st, sy):
    weight_sum = 0.0
    total_energy = 0.0
    variance = 0.0
    total_variance = 0.0
    total_power = 0.0

    # for i in range(len(st)):
    #     st[i].data = st[i].data / np.max(st[i].data)
    #     sy[i].data = sy[i].data / np.max(sy[i].data)

    for i in range(len(st) // 3):

        power = 0.0
        energy = 0.0

        # For 1 station
        for j in range(len(st[i * 3 + 0].data)):
            tmp_energy = st[i * 3 + 0].data[j] - sy[i * 3 + 0].data[j]
            energy += tmp_energy * tmp_energy
            tmp_energy = st[i * 3 + 1].data[j] - sy[i * 3 + 1].data[j]
            energy += tmp_energy * tmp_energy
            tmp_energy = st[i * 3 + 2].data[j] - sy[i * 3 + 2].data[j]
            energy += tmp_energy * tmp_energy

            power += st[i * 3 + 0].data[j] ** 2
            power += st[i * 3 + 1].data[j] ** 2
            power += st[i * 3 + 2].data[j] ** 2

        # resume 1 station info
        weight_sum += float(st[i * 3 + 0].stats.W)
        total_energy += energy
        variance += float(st[i * 3 + 0].stats.W) * energy
        total_power += power
        total_variance += float(st[i * 3 + 0].stats.W) * power
        energy /= power
        var = (1.0 - energy) * 100.0
        st[i * 3 + 0].stats.VR = var
        st[i * 3 + 1].stats.VR = var
        st[i * 3 + 2].stats.VR = var

    # end single stations checkfit loop
    # begin general solution
    vred = (1.0 - total_energy) * 100.0
    variance /= weight_sum
    total_variance /= weight_sum
    variance /= total_variance

    # THIS is he VR quantity of the original code.
    # Keep this value for compatibility
    # The same if for station VR in the loop above
    variance = (1.0 - variance) * 100.0

    # Set Quality
    if variance < 20.0:
        quality = 0
    elif 20.0 <= variance < 40.0:
        quality = 1
    elif 40.0 <= variance < 60.0:
        quality = 2
    elif 60.0 <= variance < 80.0:
        quality = 3
    else:
        quality = 4

    return st, sy, variance, quality

def pack_meta_moment_tensor(M, Mo, Mw, Pdc, Pclvd, EigVal, EigVec, T, N, P, np1, np2, var, quality, k):
    meta = np.zeros(39)
    meta[0] = Mo
    meta[1] = Mw
    meta[2] = Pdc
    meta[3] = Pclvd

    for i in range(6):
        meta[i + 4] = M[i] * k

    for i in range(3):
        meta[i + 10] = EigVal[i]

    flat = [x for sublist in EigVec for x in sublist]
    for i in range(9):
        meta[i + 13] = flat[i]

    meta[22] = T.val
    meta[23] = T.dip
    meta[24] = T.strike
    meta[25] = N.val
    meta[26] = N.dip
    meta[27] = N.strike
    meta[28] = P.val
    meta[29] = P.dip
    meta[30] = P.strike

    meta[31] = np1[0]
    meta[32] = np1[1]
    meta[33] = np1[2]
    meta[34] = np2[0]
    meta[35] = np2[1]
    meta[36] = np2[2]

    meta[37] = var
    meta[38] = quality

    return meta
