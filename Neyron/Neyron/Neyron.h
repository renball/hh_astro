#pragma once
#include<math.h>

#define Nn 4

#define V f[0]
#define m f[1]
#define n f[2]
#define h f[3]

class Neyron
{
private:	
	const double C = 1; // muF/cm^2

	const double g_K = 35; // mS/cm^2
	const double g_Na = 40; // mS/cm^2
	const double g_L = 0.3; // mS/cm^2

	const double E_K = -77; // mV
	const double E_Na = 55; // mV
	const double E_L = -65; // mV

    double Isyn;
    double E_syn = 0;
    double Theta_syn = 0;
    double k_syn = 0.2;
    double g_syn = 0.0;


public:
    double f[Nn];

    double alpha_n(double f[Nn])
    {
        return 0.02 * (V - 25) / (1 - exp(-(V - 25) / 9));
    }

    double beta_n(double f[Nn])
    {
        return -0.002 * (V - 25) / (1 - exp((V - 25) / 9));
    }

    double alpha_m(double f[Nn])
    {
        return 0.182 * (V + 35) / (1 - exp(-(V + 35) / 9));
    }

    double beta_m(double f[Nn])
    {
        return -0.124 * (V + 35) / (1 - exp((V + 35) / 9));
    }

    double alpha_h(double f[Nn])
    {
        return 0.25 * exp(-(V + 90) / 12);
    }

    double beta_h(double f[Nn])
    {
        return 0.25 * exp((V + 62) / 6) / exp((V + 90) / 12);
    }

    double HodgkinHuxley(int i, double f[Nn], double I,
                         int n_ind, double V_1)
    {
        switch (i)
        {
        case 0:
            if ((n_ind % 2) == 0) {
                return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * n * (E_K - V) + g_L * (E_L - V) + I) / C);
            }
            else
                Isyn = g_syn * (E_syn - V) / (1 + exp(-((V_1)-Theta_syn) / k_syn));
                return 1000 * ((g_Na * pow(m, 3) * h * (E_Na - V) + g_K * n * (E_K - V) + g_L * (E_L - V) + I+Isyn) / C);
        case 1:
            return 1000 * (alpha_m(f) * (1 - m) - beta_m(f) * m);

        case 2:
            return 1000 * (alpha_n(f) * (1 - n) - beta_n(f) * n);

        case 3:
            return 1000 * (alpha_h(f) * (1 - h) - beta_h(f) * h);
        }
        return 0;
    }

    void RungeKutta(double dt, double f[Nn], double f_next[Nn], double I,
                    double& cmax, double& lv, double& T_period, double& lt, double ct,
                    int n_ind, double V_1)
    {
        double k[Nn][4];

        if (HodgkinHuxley(0, f, I,n_ind,V_1) <= 0 && lv > 0 && V > 0) {
            cmax++;
            if (cmax > 1) {
                T_period = ct - lt;
            }
            lt = ct;
        }

        // k1
        for (int i = 0; i < Nn; i++)
            k[i][0] = HodgkinHuxley(i, f, I, n_ind, V_1) * dt;

        double phi_k1[Nn];
        for (int i = 0; i < Nn; i++)
            phi_k1[i] = f[i] + k[i][0] / 2;

        // k2
        for (int i = 0; i < Nn; i++)
            k[i][1] = HodgkinHuxley(i, phi_k1, I, n_ind, V_1) * dt;

        double phi_k2[Nn];
        for (int i = 0; i < Nn; i++)
            phi_k2[i] = f[i] + k[i][1] / 2;

        // k3
        for (int i = 0; i < Nn; i++)
            k[i][2] = HodgkinHuxley(i, phi_k2, I, n_ind, V_1) * dt;

        double phi_k3[Nn];
        for (int i = 0; i < Nn; i++)
            phi_k3[i] = f[i] + k[i][2] / 2;

        // k4
        for (int i = 0; i < Nn; i++)
            k[i][3] = HodgkinHuxley(i, phi_k3, I, n_ind, V_1) * dt;

        for (int i = 0; i < Nn; i++)
            f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;


        lv = HodgkinHuxley(0, f, I, n_ind, V_1);
    }

    void CopyArray(double source[Nn], double target[Nn])
    {
        for (int i = 0; i < Nn; i++)
            target[i] = source[i];
    }
};

