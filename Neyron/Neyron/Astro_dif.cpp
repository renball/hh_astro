#include "Astrocite.h"
#include "Astro_dif.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void astro_dif(double t, double dt, double tmax) {
	fstream tx_LV, tx_R;

	const int count = 1;
	Astrocite a1[count];

	tx_LV.open("last_values_A.txt");
	//0.07 0.16 0.67
	for (int i = 0; i < count; i++) {
		tx_LV >> a1[i].Ca >> a1[i].IP3 >> a1[i].z;
	}
	tx_LV.close();

	for (int i = 0; i < count; i++) {
		a1[i].Ca = 0.07; a1[i].IP3 = 0.16; a1[i].z = 0.67;
	}


	tx_R.open("results.txt");
	t = 0.0;
	while (t <= tmax) {
		double fa_next[Na];
		for (int i = 0; i < count; i++) {
			a1[i].RungeKutta(dt, a1[i].f, fa_next);
			a1[i].CopyArray(fa_next, a1[i].f);
			cout << fixed;
			cout << setprecision(6);
			cout << t <<" "<< i << " " << a1[i].Ca << " " << a1[i].IP3 << " " << a1[i].z << '\n';
			tx_R << t << " " << a1[i].Ca << " " << a1[i].IP3 << " " << a1[i].z << '\n';
		}
		t += dt;

	}
	tx_R.close();

	tx_LV.open("last_values_A.txt");
	for (int i = 0; i < count; i++) {
		tx_LV << a1[i].Ca << " " << a1[i].IP3 << " " << a1[i].z << '\n';
	}
	tx_LV.close();
}