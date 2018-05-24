#include <iostream>
#include <fstream>

	// Arrays created to hold the values I need for this test.
	float xmomentum[4] = {0.162000, 0.111000, -0.090000, -0.073000};
	float ymomentum[4] = {0.272000, 0.309000, -0.053000, 0.105000};
	float xposition[4] = {0.18};
	float yposition[4] = {-0.20};
	float times[4] = {0.25, 0.27, 0.28, 0.29};
	float partonmass = 0.200;
	float dt = 0.0005; // Timestep in fm/c
	float Nstep = 600;
	float xvelocity[4];
	float yvelocity[4];

// This function is used to determine the stage that will be used for each part of the calculation.
int getStage(float actualtime) {

	if (actualtime < times[0]) return -999;

	for (int i = 0; i < 4; i++) {

		if (actualtime > times[i] && actualtime <= times[i+1]) {

			return i;
		}

		else if (actualtime > times[3]) {
			return 3;
		}
	}

	return 0;
}

void filetest () {

	for (int i = 0; i < Nstep; i++) {

		float actualtime = i * dt;

		// This is where the stage is determined.
		int stage = getStage(actualtime);

		if (stage < 0) {
			continue;
		}

		// This is to reset the time for when the parton is born. This is also used to reset the time at each scattering.
		actualtime = actualtime - times[stage];

		// Calculation of the velocity at each stage.
		xvelocity[stage] = xmomentum[stage]/partonmass;
		yvelocity[stage] = ymomentum[stage]/partonmass;

		// Determine the initial position of the parton.
		float x0, y0;
		if (stage == 0) {

			x0 = xposition[stage];
			y0 = yposition[stage];
		}
		else {
			x0 = xposition[stage-1] + xvelocity[stage-1] * (times[stage] - times[stage-1]);
			y0 = yposition[stage-1] + yvelocity[stage-1] * (times[stage] - times[stage-1]);

			xposition[stage] = x0;
			yposition[stage] = y0;
		}

		// Calculating the positions as a function of time.
		float xt = xposition[stage] + (xvelocity[stage] * actualtime);
		float yt = yposition[stage] + (yvelocity[stage] * actualtime);

		cout << "{" << xt << "," << yt << "}," << endl;
	}
}