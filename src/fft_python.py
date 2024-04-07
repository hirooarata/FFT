# CC0
# A program of FFT
# version 5.0 2022/11/06 19:30 Python, BUG
# version 5.1 2022/11/06 19:30 Python, OK
import math

import numpy as np


def dif_fft(flg: float, x: np.ndarray) -> np.ndarray:
    #      ver 1.0
    #       dif_fft (frequeny-inplace type)
    #         The algorithm of FFT is referred from  "Discrete-Time Signal Processing"
    #         by A.V.Oppenheim, R.W.Schafer, p.599
    #         prentice-hall,Englewood cliffs, New Jersey 07632
    #         ISBN 013216292X
    #       --- flg= -1.0 .. regular FFT transform
    #       ---    = +1.0 .. inverse FFT transform
    #       --- x[n] .. data (complex)
    n: int = len(x)
    nu: int = math.floor(math.log2(n))
    # print("dif_fft:n=", n, "nu=log2(n)=", nu)

    # butterfly
    for ll in range(nu):
        le: int = 2 ** (nu - ll)
        le1: int = le // 2
        u: complex = 1 + 1j * 0
        w: complex = math.cos(math.pi / le1) + 1j * math.sin(flg * math.pi / le1)
        for j in range(le1):
            for i in range(j, n, le):
                ip = i + le1
                xp: complex = x[i] + x[ip]
                x[ip] = (x[i] - x[ip]) * u
                x[i] = xp
            u = u * w

    # bit reverse
    nv2 = n // 2
    j: int = 0
    for i in range(n - 1):
        if i < j:
            x[j], x[i] = x[i], x[j]
        k: int = nv2
        while k <= j:
            j -= k
            k //= 2
        j += k

    # normalize for ifft
    if flg >= 0:
        for i in range(n):
            x[i] /= n
    return x


#
def fft(x: np.ndarray):
    x = dif_fft(-1, x)
    return x


def ifft(x: np.ndarray):
    x = dif_fft(1, x)
    return x


def fft_shift(x: np.ndarray)->np.ndarray:
    nn2 = x.size
    nn = nn2 // 2
    for k in range(nn):
        x[k + nn], x[k] = x[k], x[k + nn]
    return x


def fft_checking(nn2: int):
    x = np.zeros(nn2, dtype=np.complex128)
    y = np.zeros(nn2, dtype=np.complex128)
    z = np.zeros(nn2, dtype=np.complex128)
    x[1] = 1 + 1j
    x[nn2-1] = 1 + 1j
    print("\n\nx=", x)
    x0 = x.copy()
    y = dif_fft(1.0, x).copy()
    p = fft_shift(y)
    y = fft_shift(p)
    print("\n\ny=", y)
    z = dif_fft(-1.0, y).copy()
    print("\nz=", z, "\nx=", x0)
    print("\nz-x", z - x0)
    print("up/down:loopback_error=", f"{np.linalg.norm(z - x0):.1e}")


if __name__ == "__main__":
    fft_checking(8192)
