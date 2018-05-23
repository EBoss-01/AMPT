#include <iostream>
#include <fstream>

void filetest () {

	// Arrays created to hold the values I need for this test.
	float xmomentum[4] = {0.162000, 0.111000, -0.090000, -0.073000};
	float ymomentum[4] = {0.272000, 0.309000, -0.053000, 0.105000};
	float xposition[4] = {0.18, 0.18, 0.18, 0.18};
	float yposition[4] = {-0.20, -0.20, -0.20, -0.20};
	float time[4] = {0.25, 0.27, 0.28, 0.29};
	float partonmass = 0.200;
	float dt = 0.05; // Timestep in fm/c
	float xvelocity[4];
	float yvelocity[4];
	float timing[40];
	float times[40] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39};

	for (int i = 0; i < 40; i++) {
		timing[i] = times[i]/40;

		//cout << timing[i] << endl;
	}

	for (int j = 0; j < 4; j++) {

		xvelocity[j] = xmomentum[j]/partonmass;
		yvelocity[j] = ymomentum[j]/partonmass;

		//cout << xvelocity[j] << endl;
		//cout << yvelocity[j] << endl;
		if (j == 3) {
			for (int k = 0; k < 40; k++) {

				if (timing[k] >= time[j]) {
					float xt2 = xposition[j] + xvelocity[j] * timing[k];
					float yt2 = yposition[j] + yvelocity[j] * timing[k];

					cout << "{" << xt2 << "," << yt2 << "}" << endl;
				}
			}
		}

		for (int i = 0; i < 40; i++) {

			if ((timing[i] <= time[j]) && (timing[i] >= time[j - 1])) {

				float xt = xposition[j] + xvelocity[j] * timing[i];
				float yt = yposition[j] + yvelocity[j] * timing[i];

				cout << "{" << xt << "," << yt << "}" << endl;
			}

			else {
				continue;
			}
		}



	}
}