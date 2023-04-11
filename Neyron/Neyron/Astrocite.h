#pragma once
#include<math.h>

#define Na 3

#define Ca  f[0]
#define IP3 f[1]
#define z   f[2]

class Astrocite
{
private:
	double c0 = 2;
	double c1 = 0.185;
	double v1 = 6;
	double v2 = 0.11;
	double v3 = 2.2 ;
	double v5 = 0.025 ;
	double v6 = 0.2 ;
	double k1 = 0.5 ;
	double k2 = 1 ;
	double k3 = 0.1;
	double k4 = 1.1 ;
	double a2 = 0.14 ;
	double d1 = 0.13;
	double d2 = 1.049;
	double d3 = 0.9434;
	double d5 = 0.082;
	double aG = 25;
	double bG = 500;
	double a= 0.8;
	double tIP3 = 7.143;
	double IP3X=0.16;
	double dCa = 0.001;
	double dIP3 = 0.12;


	double v4 = 0.5;

public:

	double f[Na];

	double Jchannel() {
		return c1 * v1 * IP3 * IP3 * IP3 * Ca * Ca * Ca * z * z * z * ((c0 / c1) - (1 + 1 / c1) * Ca)/(((IP3+d1)*(Ca+d5))* ((IP3 + d1) * (Ca + d5))* ((IP3 + d1) * (Ca + d5)));
	}

	double Jpump() {
		return v3 * Ca * Ca / (k3 * k3 + Ca * Ca);
	}
	double JLeak() {
		return c1 * v2 * (c0 / c1 - (1 + 1 / c1) * Ca);
	}
	double Jin() {
		return v5 + v6 * IP3 * IP3 / (k2 * k2 + IP3 * IP3);
	}
	double Jout() {
		return k1 * Ca;
	}
	double JPLC() {
		return v4 * (Ca + (1 - a) * k4) / (Ca + k4);
	}

	double U_J(int i, double f[Na]) {
		double Jc; Jc = Jchannel();
		double Jp;  Jp = Jpump();
		double Jl;  Jl = JLeak();
		double Ji; Ji = Jin();
		double Jo;  Jo = Jout();
		double J_plc;  J_plc = JPLC();

		switch (i)
		{
		case 0:
			return Jc - Jp + Jl + Ji + Jo;

		case 1:
			return ((IP3X - IP3) / tIP3) + J_plc;

		case 2:
			return (a2 * (d2 * ((IP3 + d1) / (IP3 + d3))) * (1 - z) - Ca * z);
		}
		return 0;
	}

	void RungeKutta(double dt, double f[Na], double f_next[Na])
	{
		double k[Na][4];

		// k1
		for (int i = 0; i < Na; i++)
			k[i][0] = U_J(i, f) * dt;

		double phi_k1[Na];
		for (int i = 0; i < Na; i++)
			phi_k1[i] = f[i] + k[i][0] / 2;

		// k2
		for (int i = 0; i < Na; i++)
			k[i][1] = U_J(i, phi_k1) * dt;

		double phi_k2[Na];
		for (int i = 0; i < Na; i++)
			phi_k2[i] = f[i] + k[i][1] / 2;

		// k3
		for (int i = 0; i < Na; i++)
			k[i][2] = U_J(i, phi_k2) * dt;

		double phi_k3[Na];
		for (int i = 0; i < Na; i++)
			phi_k3[i] = f[i] + k[i][2] / 2;

		// k4
		for (int i = 0; i < Na; i++)
			k[i][3] = U_J(i, phi_k3) * dt;

		for (int i = 0; i < Na; i++)
			f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
	}

	void CopyArray(double source[Na], double target[Na])
	{
		for (int i = 0; i < Na; i++)
			target[i] = source[i];
	}


};

