import numpy as np
def get_win(winlen=480, hop_len=160):
    ans = np.sqrt(np.hanning(winlen + 2))
    WA = ans[1:-1]

    NormsWA = np.zero(3 * winlen)
    for i in range(0, winlen * 2, hop_len):
        for j in range(winlen):
            NormsWA[i + 1 + j] += WA[j] ** 2

    NormsWA = NormsWA[winlen + 1 : 2 * winlen + 1]
    NormsWA[NormsWA <= 0] = 1
    WS = WA / NormsWA

    return WA, WS
