#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
using namespace std;

#define M_PI 3.14159265358979323846
#define Del 0.00000001

void F(double x, double y, double px, double py, double alpha, double& a, double& b, double& pa, double& pb, double& B0)
{
    a = y;
    b = py - x * exp(-alpha * x);
    pa = py * (1 - alpha * x) * exp(-alpha * x);
    pb = -px;
    B0 = py * py;
}

double Change_h(double h, double tol, double err)
{
    double facmax = 1.5;
    double facmin = 0.8;
    if (err < Del) return h * facmax;
    else return h * min(facmax, max(facmin, 0.9 * pow(tol / err, 1.0 / 6)));
}

void SOL(double X0, double Y0, double PX0, double PY0, double T, double alpha, double tol, FILE* file, double& a, double& b, double& pa, double& pb, int q)
{
    //printf("1\n");
    double h = 0.01;
    double err;
    double Y = Y0, YY = Y0, X = X0, XX = X0, PX = PX0, PXX = PX0, PY = PY0, PYY = PY0;
    double kx1, ky1, kpx1, kpy1,kb01, kx2, ky2, kpx2, kpy2, kb02, kx3, ky3, kpx3, kpy3, kb03, kx4, ky4, kpx4, kpy4, kb04, kx5, ky5, kpx5, kpy5, kb05, kx6, ky6, kpx6, kpy6, kb06, kx7, ky7, kpx7, kpy7, kb07;
    double y = Y0;
    double x = X0;
    double px = PX0;
    double py = PY0;
    double t = -h;
    double b0 = X0 * X0;
    double B0 = X0 * X0;
    double B00 = X0 * X0;
    double xk, yk, pxk, pyk;
    while (t < T)
    {
        err = sqrt((X - XX) * (X - XX) + (Y - YY) * (Y - YY)+(PX-PXX)*(PX-PXX)+(PY-PYY)*(PY-PYY));
        if (err > tol)
        {
           // printf("%.15lf  >  %.15lf\n", err, tol);
            h = Change_h(h, tol, err);
        }
        else
        { 
            y = Y;
            x = X;
            px = PX;
            py = PY;
            b0 = B0;
            t += h;
            h = Change_h(h, tol, err);
            if (t + h > T)
            {
                h = T - t;
            }
            if (q) fprintf(file, "%.10f %.10f %.10f %.10f %.10f %.10f\n", t, x, y, px, py, B0);
        }
       // printf("%.15lf \n", err);
      //  printf("%.5lf    %.5lf\n", x, y);
        xk = x;
        yk = y;
        pxk = px;
        pyk = py;
        F(xk, yk, pxk, pyk, alpha, kx1, ky1, kpx1, kpy1, kb01);
        yk = y + h * (ky1 * 1.0 / 5);
        xk = x + h * (kx1 * 1.0 / 5);
        pyk = py + h * (kpy1 * 1.0 / 5);
        pxk = px + h * (kpx1 * 1.0 / 5);
        F(xk, yk, pxk, pyk, alpha, kx2, ky2, kpx2, kpy2, kb02);
        yk = y + h * (ky1 * 3.0 / 40 + ky2 * 9.0 / 40);
        xk = x + h * (kx1 * 3.0 / 40 + kx2 * 9.0 / 40);
        pyk = py + h * (kpy1 * 3.0 / 40 + kpy2 * 9.0 / 40);
        pxk = px + h * (kpx1 * 3.0 / 40 + kpx2 * 9.0 / 40);
        F(xk, yk, pxk, pyk, alpha, kx3, ky3, kpx3, kpy3, kb03);
        yk = y + h * (ky1 * 44.0 / 45 - ky2 * 56.0 / 15 + ky3 * 32.0 / 9);
        xk = x + h * (kx1 * 44.0 / 45 - kx2 * 56.0 / 15 + kx3 * 32.0 / 9);
        pyk = py + h * (kpy1 * 44.0 / 45 - kpy2 * 56.0 / 15 + kpy3 * 32.0 / 9);
        pxk = px + h * (kpx1 * 44.0 / 45 - kpx2 * 56.0 / 15 + kpx3 * 32.0 / 9);
        F(xk, yk, pxk, pyk, alpha, kx4, ky4, kpx4, kpy4, kb04);
        yk = y + h * (ky1 * 19372.0 / 6561 - ky2 * 25360.0 / 2187 + ky3 * 64448.0 / 6561 - ky4 * 212.0 / 729);
        xk = x + h * (kx1 * 19372.0 / 6561 - kx2 * 25360.0 / 2187 + kx3 * 64448.0 / 6561 - kx4 * 212.0 / 729);
        pyk = py + h * (kpy1 * 19372.0 / 6561 - kpy2 * 25360.0 / 2187 + kpy3 * 64448.0 / 6561 - kpy4 * 212.0 / 729);
        pxk = px + h * (kpx1 * 19372.0 / 6561 - kpx2 * 25360.0 / 2187 + kpx3 * 64448.0 / 6561 - kpx4 * 212.0 / 729);
        F(xk, yk, pxk, pyk, alpha, kx5, ky5, kpx5, kpy5, kb05);
        yk = y + h * (ky1 * 9017.0 / 3168 - ky2 * 355.0 / 33 + ky3 * 46732.0 / 5247 + ky4 * 49.0 / 176 - ky5 * 5103.0 / 18656);
        xk = x + h * (kx1 * 9017.0 / 3168 - kx2 * 355.0 / 33 + kx3 * 46732.0 / 5247 + kx4 * 49.0 / 176 - kx5 * 5103.0 / 18656);
        pyk = py + h * (kpy1 * 9017.0 / 3168 - kpy2 * 355.0 / 33 + kpy3 * 46732.0 / 5247 + kpy4 * 49.0 / 176 - kpy5 * 5103.0 / 18656);
        pxk = px + h * (kpx1 * 9017.0 / 3168 - kpx2 * 355.0 / 33 + kpx3 * 46732.0 / 5247 + kpx4 * 49.0 / 176 - kpx5 * 5103.0 / 18656);
        F(xk, yk, pxk, pyk, alpha, kx6, ky6, kpx6, kpy6, kb06);
        yk = y + h * (ky1 * 35.0 / 384 + ky3 * 500.0 / 1113 + ky4 * 125.0 / 192 - ky5 * 2187.0 / 6784 + ky6 * 11.0 / 84);
        xk = x + h * (kx1 * 35.0 / 384 + kx3 * 500.0 / 1113 + kx4 * 125.0 / 192 - kx5 * 2187.0 / 6784 + kx6 * 11.0 / 84);
        pyk = py + h * (kpy1 * 35.0 / 384 + kpy3 * 500.0 / 1113 + kpy4 * 125.0 / 192 - kpy5 * 2187.0 / 6784 + kpy6 * 11.0 / 84);
        pxk = px + h * (kpx1 * 35.0 / 384 + kpx3 * 500.0 / 1113 + kpx4 * 125.0 / 192 - kpx5 * 2187.0 / 6784 + kpx6 * 11.0 / 84);
        F(xk, yk, pxk, pyk, alpha, kx7, ky7, kpx7, kpy7, kb07);
        Y = y + h * (ky1 * 35.0 / 384 + ky3 * 500.0 / 1113 + ky4 * 125.0 / 192 - ky5 * 2187.0 / 6784 + ky6 * 11.0 / 84);
        YY = y + h * (ky1 * 5179.0 / 57600 + ky3 * 7571.0 / 16695 + ky4 * 393.0 / 640 - ky5 * 92097.0 / 339200 + ky6 * 187.0 / 2100 + ky7 * 1.0 / 40);
        X = x + h * (kx1 * 35.0 / 384 + kx3 * 500.0 / 1113 + kx4 * 125.0 / 192 - kx5 * 2187.0 / 6784 + kx6 * 11.0 / 84);
        XX = x + h * (kx1 * 5179.0 / 57600 + kx3 * 7571.0 / 16695 + kx4 * 393.0 / 640 - kx5 * 92097.0 / 339200 + kx6 * 187.0 / 2100 + kx7 * 1.0 / 40);
        PY = py + h * (kpy1 * 35.0 / 384 + kpy3 * 500.0 / 1113 + kpy4 * 125.0 / 192 - kpy5 * 2187.0 / 6784 + kpy6 * 11.0 / 84);
        PYY = py + h * (kpy1 * 5179.0 / 57600 + kpy3 * 7571.0 / 16695 + kpy4 * 393.0 / 640 - kpy5 * 92097.0 / 339200 + kpy6 * 187.0 / 2100 + kpy7 * 1.0 / 40);
        PX = px + h * (kpx1 * 35.0 / 384 + kpx3 * 500.0 / 1113 + kpx4 * 125.0 / 192 - kpx5 * 2187.0 / 6784 + kpx6 * 11.0 / 84);
        PXX = px + h * (kpx1 * 5179.0 / 57600 + kpx3 * 7571.0 / 16695 + kpx4 * 393.0 / 640 - kpx5 * 92097.0 / 339200 + kpx6 * 187.0 / 2100 + kpx7 * 1.0 / 40);
        B0 = b0 + h * (kb01 * 35.0 / 384 + kb03 * 500.0 / 1113 + kb04 * 125.0 / 192 - kb05 * 2187.0 / 6784 + kb06 * 11.0 / 84);
        B00 = b0 + h * (kb01 * 5179.0 / 57600 + kb03 * 7571.0 / 16695 + kb04 * 393.0 / 640 - kb05 * 92097.0 / 339200 + kb06 * 187.0 / 2100 + kb07 * 1.0 / 40);
    }
    //printf("2\n");
    a = Y;
    b = X;
    pa = PY;
    pb = PX;
    printf("%.10f\n", b0);
}

int Matrix(double a11, double a12, double a21, double a22, double b1, double b2, double& x1, double& x2)
{
    double D;
    D = a11 * a22 - a12 * a21;
    if (fabs(D) < Del)
    {
        return 0;
    }
    else
    {
        x1 = (a22 * b1 - a12 * b2) / D;
        x2 = (a11 * b2 - a21 * b1) / D;
        return 1;
    }
}

void Method(double X0, double Y0, double PX0, double PY0, double T0, double alpha, double tol, FILE* file)
{
    double a, b, pa, pb, a1, b1, pa1, pb1;
    double a11, a12, a21, a22;
    double PX = PX0, PY = PY0;
    double PXX, PYY;
    double g = 1.0;
    int i = 1;
    double toll = 0.00000001;
    int n = 0;
    double h;
   // printf("1\n");
   // printf("%lf\n", T0);
    SOL(PX, Y0, PX, PY, T0, alpha, toll, file, a, b, pa, pb, 0);
    h = (a - 1) * (a - 1) + b * b;
   // printf("3\n");
    while (h > tol)
    {
        printf( "Norm=%.10lf\n", h);
       // fprintf(stdout, "n=%d T=%.10lf X=%.10lf Y=%.10lf\n", n, T, b, a);
        n++;
        SOL(PX+Del, Y0, PX + Del, PY, T0, alpha, toll, file, a1, b1, pa1, pb1, 0);
        a11 = -(b1 - b) / Del;
        a21 = -(a1 - a) / Del;
        SOL(PX, Y0, PX, PY + Del, T0, alpha, toll, file, a1, b1, pa1, pb1, 0);
        a12 = -(b1 - b) / Del;
        a22 = -(a1 - a) / Del;
       // fprintf(stdout, "\n%.5lf  %.5lf \n%.5lf  %.5lf\n\n", a11, a21, a12, a22);
        if (!Matrix(a11, a12, a21, a22, b, a - 1, PXX, PYY))
        {
            PXX = 1;
            PYY = 1;
        }
        printf("X=%.10lf    Y=%.10lf\n", b, a);
        i = 1;
        do
        {
            g = 1.0 / i;
            SOL(PX, Y0, PX + g * PXX, PY + g * PYY, T0, alpha, toll, file, a1, b1, pa1, pb1, 0);
            i++;
            if (i == 30) break;
        } while (b1 * b1 - b * b + (a1 - 1) * (a1 - 1) - (a - 1) * (a - 1) > 0);
        PX += g * PXX;
        PY += g * PYY;
        SOL(PX, Y0, PX, PY, T0, alpha, toll, file, a, b, pa, pb, 0);
        h = b * b + (a - 1) * (a - 1);
    }
    SOL(PX, Y0, PX, PY, T0, alpha, toll, file, a, b, pa, pb, 1);
    printf( "%.10lf %.10lf\n", PX, PY);

}

int main(void)
{
    FILE* file1 = fopen("Zn25.txt", "w+");
    Method(2, 0, 2, 1, M_PI / 2.0, 0.0, 0.00001, file1);
    fclose(file1);
    return 0;
}