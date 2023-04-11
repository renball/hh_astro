#include "Neyron.h"
#include "Astrocite.h"
#include <iostream>
#include <fstream>
#include "Astro_dif.h"
#include "Neyro_dif.h"

using namespace std;

void main() {
	int in;
	double t = 0.0; double dt = 0.00005; double tmax = 1.0;
	cin >> in;

	if (in == 1) {
		neyro_dif(t, dt, tmax);
	}
	else if (in == 2) {
		astro_dif(t, dt, tmax);
	}
	cout << "_____________________________________";
}