//-----------------------------------------------------------------------------------------
//This code will be used to read in the file outputs from AMPT and parse them together.
//This will allow me to reconstruct the parton evolution and make a scattering animation.
//
//File 1: parton-initial-afterPropagation_1.dat
//File 2: parton-collisionHistory_1.dat
//
//05-21-18
//------------------------------------------------------------------------------------------

#include "TLatex.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TRegexp.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TPaletteAxis.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//---------------------------------
//Structure for the Parton
//---------------------------------

struct parton {
// Parton initial momentum
	float px;
	float py;
	float pz;

// Parton initial location
	float x;
	float y;
	float z;
	float t;

//Parton mass in GeV
	float m;

};

//----------------------------------
// Variables
//----------------------------------

// This is the counter to keep track of which event is being read.
int eventnumber = 0;

// Vector for the scattering times for the specific partons.
vector<parton> ScatteringTime;

// A vector of other vectors hopefully
vector<vector<parton> > EventPartons;

// This is my practice vector to be certain that my filling is working as intended.
vector<parton> PracticeVector;

// This vector will store the coordinates for each parton at every time. 
vector<float> xvals;
vector<float> yvals;
vector<int> goesboom;

// Time step in fm/c
float dt = 0.05;
float NStep = 100;

//----------------------------------
//Functions
//----------------------------------

void myText(Double_t x,Double_t y,Color_t color,const char *text,Double_t tsize = 0.05,double angle=-1) {
	TLatex l;
	l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextColor(color);
	if (angle > 0) l.SetTextAngle(angle);
	l.DrawLatex(x,y,text);
}

// This function is used to determine the stage that will be used for each part of the calculation.
int getStage(float actualtime, vector<parton> v) {

	if (actualtime < v[0].t) return -999;

	for (int i = 0; i < v.size(); i++) {

		if (actualtime > v[i].t && actualtime <= v[i+1].t) {

			return i;
		}
	}

	if (actualtime > v[v.size()-1].t) {
		return v.size()-1;
	}

	return 0;
}

// This function is used to produce all the different histograms for each time step.
void draw(vector<float> x, vector<float> y, int iterate, vector<int> bigboom) {

	cout << "+++++++++++++++++" << iterate << endl;

	TCanvas* c = new TCanvas(Form("c_%i", iterate),Form("c_%i", iterate),480,480);
	gStyle->SetOptStat(0);

	c->SetTickx();
	c->SetTicky();

	TH2F* hTimestep = new TH2F(Form("hTimestep_[%i]", iterate), Form("hTimestep_%i", iterate), 100, -2.5, 2.5, 100, -2.5, 2.5);
	hTimestep->SetTitle("");
	hTimestep->GetXaxis()->SetTitle("x [fm]");
	hTimestep->GetYaxis()->SetTitle("y [fm]");
	hTimestep->Draw();

	for (int j = 0; j < x.size(); j++) {
		TEllipse *tell = new TEllipse(x[j], y[j],0.09,0.09);
		tell->SetFillColor(kBlue);
		tell->Draw("same");

		if (bigboom[j] > 0) {
			tell->SetFillColor(kRed);
			tell->Draw("same");
		}
	}

	if (iterate < 10) {
		c->SaveAs(Form("frame_test/Iteration_00%i.png", iterate));
	}
	else if (iterate >= 10 && iterate < 100) {
		c->SaveAs(Form("frame_test/Iteration_0%i.png", iterate));
	}

}

// This function is used to do the position calculations and will return x,y coordinates at every time.
void calculatePosition (float actualtime, vector<parton> v, float &xt, float &yt, int &boom) {

	float xposition[v.size()];
	float yposition[v.size()];
	float xvelocity[v.size()];
	float yvelocity[v.size()];
	float x0, y0;

	for (int i = 0; i < v.size(); i++) {
		// Calculation of the velocity at each stage.
		float energy = TMath::Sqrt(pow(v[i].px,2) + pow(v[i].py,2) + pow(v[i].pz,2) + pow(v[i].m,2));
		TLorentzVector ev(v[i].px, v[i].py, v[i].pz, energy);
		float beta = ev.Beta();
		float phi = ev.Phi();
		xvelocity[i] = beta * TMath::Cos(phi);
		yvelocity[i] = beta * TMath::Sin(phi);


		// Determine the initial position of the parton.
		if (i == 0) {

			x0 = v[i].x;
			y0 = v[i].y;

			xposition[i] = x0;
			yposition[i] = y0;
		}
		else {
			//cout << "+++ " << stage << "   " << xvelocity[stage-1] << endl;

			x0 = xposition[i-1] + xvelocity[i-1] * (v[i].t - v[i-1].t);
			y0 = yposition[i-1] + yvelocity[i-1] * (v[i].t - v[i-1].t);

			xposition[i] = x0;
			yposition[i] = y0;
		}
	}

	// Call getStage function to determine the stages for every time.
	int stage = getStage(actualtime,v);
	if (stage < 0) {
		xt =-999;
		yt=-999;
		return;
	}

	//cout << stage << endl << endl;

	actualtime = actualtime - v[stage].t;

	//cout << "--- " << stage << "   " << xvelocity[stage] << endl;


	// Calculating the positions as a function of time.
	xt = xposition[stage] + (xvelocity[stage] * actualtime);
	yt = yposition[stage] + (yvelocity[stage] * actualtime);
	boom = stage;
}

// This function will loop over each parton and then calculate the positions at every given moment in time. 
void processEvent() {
	int iterate = 0;

	// This loops through each time.
	for (int i = 0; i < NStep; i++) {

		float actualtime = i * dt;

		// Put Parton Loop here
		
		for (int p = 0; p < EventPartons.size(); p++) {

			std::vector<parton> v = EventPartons[p];

			float xt, yt;
			int boom;
			calculatePosition(actualtime,v,xt,yt,boom);

			cout << "{" << xt << "," << yt << "}," << endl;

			xvals.push_back(xt);
			yvals.push_back(yt);
			goesboom.push_back(boom);
		}

		draw(xvals, yvals, iterate, goesboom);

		cout << "------------" << iterate << endl;

		xvals.clear();
		yvals.clear();
		goesboom.clear();
		iterate++;

	}
}

// This will be my attempt to read in the files that I will be using. 
void EventAnimation(void) {
	ifstream myInitialFileInfo;
	ifstream myEvolutionFile;

	myInitialFileInfo.open("parton-initial-afterPropagation_1.dat");

	if (!myInitialFileInfo) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file parton-initial-afterPropagation_1.dat" << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened parton-initial-afterPropagation_1.dat" << endl;
	}

	myEvolutionFile.open("parton-collisionsHistory_1.dat");

	if (!myEvolutionFile) {
	// As before this will let me know if this specific file fails to open.
		cout << "Unable to open file parton-collisionsHistory_1.dat" << endl;
		return;
	}
	else {
	// This is here simply to confirm that this file opened successfully as well.
		cout << "Successfully opened parton-collisionsHistory_1.dat" << endl;
	}

	while (myInitialFileInfo) {

		if (eventnumber % 100 == 0) {

			cout << "Reading Event Number " << eventnumber << endl;
		}

		// Information needed to read the event header.
		int iterationN;
		int nPartons;
		int nBaryons;
		int nMesons;
		int particleC;
		int particleNC;

		myInitialFileInfo >> eventnumber >> iterationN >> nPartons >> nBaryons >> nMesons >> particleC >> particleNC;

		// This line prevents the file from trying to read the last line multiple times.
		if (!myInitialFileInfo) break;

		// This loop fills in the initial information for each parton.
		for (int i = 0; i < nPartons; i++) {

			int partID;
			float momenta[3];
			float mass;
			double spacetime[4];

			myInitialFileInfo >> partID >> momenta[0] >> momenta[1] >> momenta[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

			parton partinfo;
			partinfo.px = momenta[0];
			partinfo.py = momenta[1];
			partinfo.pz = momenta[2];
			partinfo.m = mass;
			partinfo.x = spacetime[0];
			partinfo.y = spacetime[1];
			partinfo.z = spacetime[2];
			partinfo.t = spacetime[3];

			vector<parton> auxillary;
			auxillary.push_back(partinfo);
			EventPartons.push_back(auxillary);
		}

		//-------------------------------------------
		// Additional variables for evolution
		//-------------------------------------------
		string line;
		int evt;
		int junk1;
		int partonindex1;
		int partonindex2;
		// Initial parton information
		int parton1_id_initial;
		float parton1_momenta_initial[3];
		float parton1_mass_initial;
		double parton1_spacetime_initial[4];
		int parton2_id_initial;
		float parton2_momenta_initial[3];
		float parton2_mass_initial;
		double parton2_spacetime_initial[4];
		// Final parton information
		int parton1_final_id;
		float parton1_final_momenta[3];
		float parton1_final_mass;
		double parton1_final_spacetime[4];
		int parton2_final_id;
		float parton2_final_momenta[3];
		float parton2_final_mass;
		double parton2_final_spacetime[4];

		while (std::getline(myEvolutionFile,line)) {

			if (!myEvolutionFile) break;
			
			stringstream heading(line);
			string description;

			if (heading >> description >> evt >> junk1 >> partonindex1 >> partonindex2) {

				if (evt == eventnumber) {

					myEvolutionFile >> parton1_id_initial >> parton1_momenta_initial[0] >> parton1_momenta_initial[1] >> parton1_momenta_initial[2] >> parton1_mass_initial >> parton1_spacetime_initial[0] >> parton1_spacetime_initial[1] >> parton1_spacetime_initial[2] >> parton1_spacetime_initial[3];

					myEvolutionFile >> parton2_id_initial >> parton2_momenta_initial[0] >> parton2_momenta_initial[1] >> parton2_momenta_initial[2] >> parton2_mass_initial >> parton2_spacetime_initial[0] >> parton2_spacetime_initial[1] >> parton2_spacetime_initial[2] >> parton2_spacetime_initial[3];

					myEvolutionFile >> parton1_final_id >> parton1_final_momenta[0] >> parton1_final_momenta[1] >> parton1_final_momenta[2] >> parton1_final_mass >> parton1_final_spacetime[0] >> parton1_final_spacetime[1] >> parton1_final_spacetime[2] >> parton1_final_spacetime[3];

					myEvolutionFile >> parton2_final_id >> parton2_final_momenta[0] >> parton2_final_momenta[1] >> parton2_final_momenta[2] >> parton2_final_mass >> parton2_final_spacetime[0] >> parton2_final_spacetime[1] >> parton2_final_spacetime[2] >> parton2_final_spacetime[3];

					parton part1;
					part1.px = parton1_final_momenta[0];
					part1.py = parton1_final_momenta[1];
					part1.pz = parton1_final_momenta[2];
					part1.x = parton1_final_spacetime[0];
					part1.y = parton1_final_spacetime[1];
					part1.z = parton1_final_spacetime[2];
					part1.t = parton1_final_spacetime[3];
					part1.m = parton1_final_mass;

					parton part2;
					part2.px = parton2_final_momenta[0];
					part2.py = parton2_final_momenta[1];
					part2.pz = parton2_final_momenta[2];
					part2.x = parton2_final_spacetime[0];
					part2.y = parton2_final_spacetime[1];
					part2.z = parton2_final_spacetime[2];
					part2.t = parton2_final_spacetime[3];
					part2.m = parton2_final_mass;

					EventPartons[partonindex1 - 1].push_back(part1);
					EventPartons[partonindex2 - 1].push_back(part2);
				}
			}
		}

		// Call function that Processes the Event.
		processEvent();

/*
		for(int i=0; i<40; i++)
		{
			std::vector<parton> parton_history = EventPartons[i];

			cout << "---> PARTON " << i+1 << " ---------" << endl;

			for(int j=0; j<parton_history.size(); j++)
			{
				parton p = parton_history[j];

				cout << Form("px %i step = %f", j, p.px) << endl;
				cout << Form("py %i step = %f", j, p.py) << endl;
			}

		}
		*/
	}

	myEvolutionFile.close();
	myInitialFileInfo.close();
	return;
}

