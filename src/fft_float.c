/* CC0 */
#include <stdio.h>
#include <math.h>

/*------------------------------------------------------------------------------*/
void dif_fft_float(int n, float flg, float * xr, float * xi) {
    /*     VER 1.0                                                              */
    /*       dif_fft (frequeny-inplace type);                                   */
    /*         The algorithm of FFT is referred from                            */
    /*         "Discrete-Time Signal Processing"                                */
    /*         by A.V.Oppenheim, R.W.Schafer, p.599;                            */
    /*         prentice-hall,Englewood clif (fs, New Jersey 07632;              */
    /*         ISBN 013216292X;                                                 */
    /*         Translated from FORTRAN to C language                            */
    /*         No complex routine rquired                                       */
    /*      --- FLG= -1.0 .. REGULAR FFT TRANSFORM                              */
    /*      ---    = +1.0 .. INVERSE FFT TRANSFORM                              */
    /*      --- xr .. DATA (REAL PART)                                          */
    /*      --- xi .. DATA (IMAG PART)                                          */

    const float pi = 3.1415926535897932384626433832795029;
    int le, le1, k, nu, nv2, i, j, ll, ip;
    float xrp, xip, xrm, xim, ur, ui, wr, wi, tr, ti;

    nu = (int)(log2((float) n));
    /*printf("dif_fft_float:m = %d, n = %d ¥n",m,n);*/

    for (ll = 0; ll < nu; ++ll) { /*butterfly*/
        le = (int) pow(2.0, (float)(nu - ll));
        le1 = le / 2;
        ur = 1.0;
        ui = 0.0;
        wr = cos(pi / (float) le1);
        wi = sin(flg * pi / (float) le1);
        for (j = 0; j < le1; ++j) {
            for (i = j; i < n; i += le) {
                ip = i + le1;
                xrp = xr[i] + xr[ip];
                xip = xi[i] + xi[ip];
                xrm = xr[i] - xr[ip];
                xim = xi[i] - xi[ip];
                xr[ip] = xrm * ur - xim * ui;
                xi[ip] = xrm * ui + xim * ur;
                xr[i] = xrp;
                xi[i] = xip;
            }
            tr = ur * wr - ui * wi;
            ui = ur * wi + ui * wr;
            ur = tr;
        }
    }

    nv2 = n / 2; /*bit reverse*/
    for (j = 0, i = 0; i < (n - 1); ++i) {
        if (i < j) {
            tr = xr[j];
            ti = xi[j];
            xr[j] = xr[i];
            xi[j] = xi[i];
            xr[i] = tr;
            xi[i] = ti;
        }
        k = nv2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    if (flg >= 0.0) { /*normalize for ifft*/
        for (i = 0; i < n; ++i) {
            xr[i] /= (float)n;
            xi[i] /= (float)n;
        }
    }
}

int main(void) {
    int n;
    float xr[1024] = {
        0.0
    };
    float xi[1024] = {
        0.0
    };
    float xr0,xi0;
    n=1024;
    xr[1] = 1.0;
    xi[1] = 1.0;
    xr0 = xr[1];
    xi0 = xi[1];
    dif_fft_float(n, 1.0, xr, xi);
    dif_fft_float(n, -1.0, xr, xi);
    printf("\nx=[%f,%f]", xr0, xi0);
    printf("\nz=[%f,%f]", xr[1], xi[1]);
    printf("\nz-x=[%f,%f]\n", xr[1] - xr0, xi[1] - xi0);
    return 0;
}
/*----------------------------------------------------------------------------*/
