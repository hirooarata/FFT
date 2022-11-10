# FFT
## MIT LICENCE

Copyright © 2022 <hirooarata>
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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

