#include "Neyron.h"
#include "Neyro_dif.h"
#include <fstream>
#include <iostream>
using namespace std;

void neyro_dif(double t, double dt, double tmax) {
	fstream tx_LV, tx_R;
	double I = 1.4; double dI = 0.001; double Imax = 1.4;
	
	Neyron n1[2];
	for (int i = 0; i < 2; i++) {
		n1[i].V = -58.7085; n1[i].m = 0.0953; n1[i].n = 0.000913; n1[i].h = 0.3662;
	}
	while (I <= Imax) {
		tx_LV.open("last_values_N.txt");
		//-58.7085 0.0953 0.000913 0.3662
		for (int i = 0; i < 2; i++) {
			tx_LV >> n1[i].V >> n1[i].m >> n1[i].n >> n1[i].h;
		}
		tx_LV.close();

		t = 0.0;
		double cmax = 0;
		double lv = 0.0;
		double lt = 0.0;
		double T_period = 0.0;

		tx_R.open("results.txt");
		while (t <= tmax) {
			double fn_next[Nn];
			for (int i = 0; i < 2; i++) {
				double V_1 = 0.0;
				if (i != 0) {V_1 = n1[i - 1].V;}

				n1[i].RungeKutta(dt, n1[i].f, fn_next, I, cmax, lv, T_period, lt, t,i,V_1);
				n1[i].CopyArray(fn_next, n1[i].f);

				cout << i << " " << t << " " << n1[i].V << " " << n1[i].m << " " << n1[i].n << " " << n1[i].h << " " << '\n';

				tx_R << t << " " << n1[i].V << " " << n1[i].m << " " << n1[i].n << " " << n1[i].h << " " << '\n';
			}
			t += dt;
		}
		tx_R.close();

		tx_LV.open("last_values_N.txt");
		for (int i = 0; i < 2; i++) {
			tx_LV << n1[i].V << " " << n1[i].m << " " << n1[i].n << " " << n1[i].h;
		}
		tx_LV.close();
		I += 0.1;
	}
}