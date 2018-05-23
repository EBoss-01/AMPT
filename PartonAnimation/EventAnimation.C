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

				cout << "Event number = " << evt << endl;

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
	}

	myEvolutionFile.close();
	myInitialFileInfo.close();
	return;
}

