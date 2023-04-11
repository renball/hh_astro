#include "Astrocite.h"
#include "Astro_dif.h"
#include <fstream>
#include <iostream>

using namespace std;

void astro_dif(double t, double dt, double tmax) {
	fstream tx_LV, tx_R;
	Astrocite a1;
	a1.Ca = 0.07; a1.IP3 = 0.16; a1.z = 0.67;
	
	tx_LV.open("last_values_A.txt");
	//0.07 0.16 0.67
	tx_LV >> a1.Ca >> a1.IP3 >> a1.z;
	tx_LV.close();

	tx_R.open("results.txt");
	t = 0.0;
	while (t <= tmax) {
		double fa_next[Na];
		a1.RungeKutta(dt, a1.f, fa_next);
		a1.CopyArray(fa_next, a1.f);
		t += dt;
		cout << t << " " << a1.Ca << " " << a1.IP3 << " " << a1.z << '\n';
		tx_R << t << " " << a1.Ca << " " << a1.IP3 << " " << a1.z << '\n';
	}
	tx_R.close();

	tx_LV.open("last_values_A.txt");
	tx_LV << a1.Ca << " " << a1.IP3 << " " << a1.z << '\n';
	tx_LV.close();
}