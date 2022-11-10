/*
Copyright © 2022 <hirooarata>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
// fft.c
// A program of FFT
// version 1.0 2022/11/09 20:30 Python, OK;
// version 2.0 2022/11/10 20:30 Python to C, OK;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void dif_fft(int n, double flg, double complex x[]){
    //      ver 1.0;
    //       dif_fft (frequeny-inplace type);
    //         The algorithm of FFT is referred from  "Discrete-Time Signal Processing";
    //         by A.V.Oppenheim, R.W.Schafer, p.599;
    //         prentice-hall,Englewood cliffs, New Jersey 07632;
    //         ISBN 013216292X;
    //       --- flg= -1.0 .. regular FFT;
    //       ---    = +1.0 .. inverse FFT;
    //       --- x[n] .. data (complex);
    int nu = (floor)(log2((double)n));
    // printf("dif_fft:n=%d, nu=log2(n)=%d", n, nu);
    double pi= 3.1415926535897932384626433832795029;
    for ( int ll=0; ll < nu; ++ll ){ // butterfly;
        int  le = (int)pow(2.0, (double)(nu - ll));
        int  le1 = le / 2;
        double complex u  = 1+ 0I;
        double complex w = cos(pi / (double)le1)+1I* sin(flg * pi / (double)le1);
        for ( int j=0; j < le1; ++j ){
            for ( int i = j; i < n; i += le ){
                int ip = i + le1;
                double complex xp = x[i] + x[ip];
                x[ip] = (x[i] - x[ip]) * u;
                x[i] = xp;
                }
            u = u * w;
            }
        }
    
    int nv2 = n / 2; // bit reverse;
    for (int j=0, i=0;  i < (n - 1); ++i ){
        if ( i < j ){
            double complex t=x[j];
            x[j] = x[i];
            x[i]=t;
            }
        int k = nv2;
        while( k <= j ){
            j -= k;
            k /= 2;
            }
        j += k;
        }
    
    if (flg >= 0.0){  // normalize for ifft
        for ( int i = 0; i < n; ++i ){
            x[i] /=(double)n;
            }
        }
}

//
void fft(int nn2, double complex x[]){
    dif_fft(nn2, -1, x);
}

void ifft(int nn2, double complex x[]){
    dif_fft(nn2, 1, x);
}

void fft_shift(int nn2, double complex x[]){
    int nn = nn2 / 2;
    for (int k = 0; k < nn; ++k ){
        // x[k + nn], x[k] = x[k], x[k + nn];
        double complex t=x[k];
        x[k] = x[k + nn];
        x[k + nn]=t;
        }
}

void fft_copy(int nn2, double complex y[], double complex x[]){
    for (int k=0; k < nn2; ++k ){
        y[k] = x[k];
        }
}

void fft_checking(int nn2) {
    // x = np.zeros(nn2, dtype = np.complex128);
    // y = np.zeros(nn2, dtype = np.complex128);
    // z = np.zeros(nn2, dtype = np.complex128);
    double complex x[8192]={0.0+1I*0.0};
    double complex x0[8192]={0.0+1I*0.0};
    double complex y[8192]={0.0+1I*0.0};
    double complex z[8192]={0.0+1I*0.0};

    x[1] = 1 + 1I;
    x[nn2 - 1] = 1 + 1I;
    // printf("\n\nx=", x);
    fft_copy( nn2, x0, x );
    dif_fft( nn2, 1.0, x );
    fft_copy( nn2, y, x );
    fft_shift( nn2,y );
    fft_shift( nn2,y );
    fft_copy( nn2, z, y );            // printf("\n\ny=", y);
    dif_fft( nn2, -1.0, z );
    printf("\nx=[%f,%f]",  creal(x0[1]),cimag(x0[1]));
    printf("\nz=[%f,%f]", creal(z[1]), cimag(z[1]));
    printf("\nz-x=[%f,%f]", creal(z[1] - x0[1]),cimag(z[1] - x0[1]) );
    // printf("up/down:loopback_error=", f "{np.linalg.norm(z - x0):.1e}");
    }

int main(void){
    fft_checking(8192);
    return 0;
}
