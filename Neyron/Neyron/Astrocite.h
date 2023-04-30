#pragma once
#include<math.h>

#define Na 3

#define Ca  f[0]
#define IP3 f[1]
#define z   f[2]

class Astrocite
{
private:
	double c_0 = 2;
	double c_1 = 0.185;
	double v_1 = 6;
	double v_2 = 0.11;
	double v_3 = 2.2;
	double v_5 = 0.025;
	double v_6 = 0.2;
	double k_1 = 0.5;
	double k_2 = 1;
	double k_3 = 0.1;
	double k_4 = 1.1;
	double a_2 = 0.14;
	double d_1 = 0.13;
	double d_2 = 1.049;
	double d_3 = 0.9434;
	double d_5 = 0.082;
	double alpha_G = 25;
	double beta_G = 500;
	double alpha = 0.8;
	double tau_IP3 = 7.143;
	double IP3_star = 0.16;
	//double d_Ca = 0.001;
	//double d_IP3 = 0.12;

	double v_4=0.4;

public:

	double f[Na];

	double UllahJung(int i, double f[Na])
	{
		switch (i)
		{
		case 0: // Ca
			return (c_1 * v_1 * pow(IP3, 3) * pow(Ca, 3) * pow(z, 3) * (c_0 / c_1 - (1 + 1 / c_1) * Ca) / pow(((IP3 + d_1) * (Ca + d_5)), 3)) - (v_3 * pow(Ca, 2) / (pow(k_3, 2) + pow(Ca, 2))) + (c_1 * v_2 * (c_0 / c_1 - (1 + 1 / c_1) * Ca)) + (v_5 + v_6 * pow(IP3, 2) / (pow(k_2, 2) + pow(IP3, 2))) - k_1 * Ca;

		case 1: // IP3
			return (IP3_star - IP3) / tau_IP3 + v_4 * (Ca + (1 - alpha) * k_4) / (Ca + k_4);

		case 2: // z
			return a_2 * (d_2 * ((IP3 + d_1) / (IP3 + d_3)) * (1 - z) - Ca * z);
		}
		return 0;
	}

	void RungeKutta(double dt, double f[Na], double f_next[Na])
	{
		double k[Na][4];

		// k1
		for (int i = 0; i < Na; i++)
			k[i][0] = UllahJung(i, f) * dt;

		double phi_k1[Na];
		for (int i = 0; i < Na; i++)
			phi_k1[i] = f[i] + k[i][0] / 2;

		// k2
		for (int i = 0; i < Na; i++)
			k[i][1] = UllahJung(i, phi_k1) * dt;

		double phi_k2[Na];
		for (int i = 0; i < Na; i++)
			phi_k2[i] = f[i] + k[i][1] / 2;

		// k3
		for (int i = 0; i < Na; i++)
			k[i][2] = UllahJung(i, phi_k2) * dt;

		double phi_k3[Na];
		for (int i = 0; i < Na; i++)
			phi_k3[i] = f[i] + k[i][2] / 2;

		// k4
		for (int i = 0; i < Na; i++)
			k[i][3] = UllahJung(i, phi_k3) * dt;

		for (int i = 0; i < Na; i++)
			f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
	}

	void CopyArray(double source[Na], double target[Na])
	{
		for (int i = 0; i < Na; i++)
			target[i] = source[i];
	}

};

