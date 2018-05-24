#include <iostream>
#include <fstream>

void filetest () {

	// Arrays created to hold the values I need for this test.
	float xmomentum[4] = {0.162000, 0.111000, -0.090000, -0.073000};
	float ymomentum[4] = {0.272000, 0.309000, -0.053000, 0.105000};
	float xposition[4] = {0.18};
	float yposition[4] = {-0.20};
	float time[4] = {0.25, 0.27, 0.28, 0.29};
	float partonmass = 0.200;
	float dt = 0.025; // Timestep in fm/c
	float Nstep = 100;
	float actualtime = 0;
	float xvelocity[4];
	float yvelocity[4];

	for (int i = 0; i < Nstep; i++) {

		actualtime = dt * i;

		for (int j = 0; j < 4; j++) {

			xvelocity[j] = xmomentum[j]/partonmass;
			yvelocity[j] = ymomentum[j]/partonmass;

			if (j > 0) {

				xposition[j] = xposition[j-1] * time[j];

				//cout << "calculated position: " << xposition[j] << endl;
			}

			if ((j == 3) && (actualtime >= time[3])) {
				float xt2 = xposition[3] + xvelocity[3] * actualtime;
				float yt2 = yposition[3] + yvelocity[3] * actualtime;

				cout << "{" << xt2 << "," << yt2 << "}" << "," << endl;
			}


			if ((actualtime < time[j]) && (actualtime >= time[j - 1])) {

				float xt = xposition[j] + xvelocity[j] * actualtime;
				float yt = yposition[j] + yvelocity[j] * actualtime;

				cout << "{" << xt << "," << yt << "}" << "," << endl;
			}

		}

	}
}