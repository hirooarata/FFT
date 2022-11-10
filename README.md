# FFT
## FFT programs with Python and C

## dif_fft (frequeny-inplace type)

The algorithm of FFT is referred from  "Discrete-Time Signal Processing" by A.V.Oppenheim, R.W.Schafer, p.599, 
prentice-hall,Englewood cliffs, New Jersey 07632,ISBN 013216292X

You can get a FFT program with FORTRAN in this book.

## Usage(Python)

x = dif_fft(flg, x)

flg= -1.0 .. regular FFT(From Time domain to Frequency domain)

flg= +1.0 .. inverse FFT(From Frequency domain to Time domain)

x[nn2]    .. complex ndarray

## Usage(C)

dif_fft(nn2, flg, x);

nn2       .. size of array

flg= -1.0 .. regular FFT(From Time domain to Frequency domain)

flg= +1.0 .. inverse FFT(From Frequency domain to Time domain)

x[nn2]    .. double complex array

## Sample

x= [0.+0.j 1.+1.j 0.+0.j ... 0.+0.j 0.+0.j 1.+1.j]


y= [0.00024414+0.00024414j 0.00024414+0.00024414j 0.00024414+0.00024414j ...
 0.00024414+0.00024414j 0.00024414+0.00024414j 0.00024414+0.00024414j]

z= [0.00000000e+00+0.00000000e+00j 1.00000000e+00+1.00000000e+00j
 0.00000000e+00+0.00000000e+00j ... 7.49041986e-15-7.26014285e-15j
 0.00000000e+00+0.00000000e+00j 1.00000000e+00+1.00000000e+00j] 
 
x= [0.+0.j 1.+1.j 0.+0.j ... 0.+0.j 0.+0.j 1.+1.j]

z-x [ 0.00000000e+00+0.00000000e+00j  8.88178420e-15+1.95399252e-14j
  0.00000000e+00+0.00000000e+00j ...  7.49041986e-15-7.26014285e-15j
  0.00000000e+00+0.00000000e+00j -4.01900735e-14-4.78506124e-14j]

up/down:loopback_error= 8.9e-14
