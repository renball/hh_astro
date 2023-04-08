#pragma once
#include<math.h>
class Neyron
{
private:
	#define N 4

	#define V f[0]
	#define m f[1]
	#define n f[2]
	#define h f[3]

	double f[N];
	
	
	const double C = 1; // muF/cm^2

	const double g_K = 35; // mS/cm^2
	const double g_Na = 40; // mS/cm^2
	const double g_L = 0.3; // mS/cm^2

	const double E_K = -77; // mV
	const double E_Na = 55; // mV
	const double E_L = -65; // mV

public:
    double alpha_n(double f[N])
    {
        return 0.02 * (V - 25) / (1 - exp(-(V - 25) / 9));
    }

    double beta_n(double f[N])
    {
        return -0.002 * (V - 25) / (1 - exp((V - 25) / 9));
    }

    double alpha_m(double f[N])
    {
        return 0.182 * (V + 35) / (1 - exp(-(V + 35) / 9));
    }

    double beta_m(double f[N])
    {
        return -0.124 * (V + 35) / (1 - exp((V + 35) / 9));
    }

    double alpha_h(double f[N])
    {
        return 0.25 * exp(-(V + 90) / 12);
    }

    double beta_h(double f[N])
    {
        return 0.25 * exp((V + 62) / 6) / exp((V + 90) / 12);
    }

    double HodgkinHuxley(int i, double f[N], double I)
    {
        switch (i)
        {
        case 0:
            return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * n * (E_K - V) + g_L * (E_L - V) + I) / C);

        case 1:
            return 1000 * (alpha_m(f) * (1 - m) - beta_m(f) * m);

        case 2:
            return 1000 * (alpha_n(f) * (1 - n) - beta_n(f) * n);

        case 3:
            return 1000 * (alpha_h(f) * (1 - h) - beta_h(f) * h);
        }
        return 0;
    }

    void RungeKutta(double dt, double f[N], double f_next[N], double I, double& cmax, double& lv, double& T_period, double& lt, double ct)
    {
        double k[N][4];

        if (HodgkinHuxley(0, f, I) <= 0 && lv > 0 && V > 0) {
            cmax++;
            if (cmax > 1) {
                T_period = ct - lt;
            }
            lt = ct;
        }

        // k1
        for (int i = 0; i < N; i++)
            k[i][0] = HodgkinHuxley(i, f, I) * dt;

        double phi_k1[N];
        for (int i = 0; i < N; i++)
            phi_k1[i] = f[i] + k[i][0] / 2;

        // k2
        for (int i = 0; i < N; i++)
            k[i][1] = HodgkinHuxley(i, phi_k1, I) * dt;

        double phi_k2[N];
        for (int i = 0; i < N; i++)
            phi_k2[i] = f[i] + k[i][1] / 2;

        // k3
        for (int i = 0; i < N; i++)
            k[i][2] = HodgkinHuxley(i, phi_k2, I) * dt;

        double phi_k3[N];
        for (int i = 0; i < N; i++)
            phi_k3[i] = f[i] + k[i][2] / 2;

        // k4
        for (int i = 0; i < N; i++)
            k[i][3] = HodgkinHuxley(i, phi_k3, I) * dt;

        for (int i = 0; i < N; i++)
            f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;


        lv = HodgkinHuxley(0, f, I);
    }

    void CopyArray(double source[N], double target[N])
    {
        for (int i = 0; i < N; i++)
            target[i] = source[i];
    }
};

