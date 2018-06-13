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
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//---------------------------------
// TTree Creation
//---------------------------------

TFile* f1 = NULL;
TTree *tree = NULL;

//---------------------------------
//Structure for the Parton
//---------------------------------

struct parton {
// Event Number
	int evtN;
// Parton ID
	int pID;
// Parton initial momentum
	float px;
	float py;
	float pz;
// Parton scattering angle
	float sangle;
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
float Pi = 3.14159;

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

// This function calculates the scattering angle and includes the Lorentz transformation into the CoM.
void getAngle(float benergy1, float benergy2, float benergy3, float benergy4, TLorentzVector vb1, TLorentzVector vb2, TLorentzVector vb3, TLorentzVector vb4, TLorentzVector com, float parton1_momenta_initial[], float parton2_momenta_initial[], float parton1_mass_initial, float parton2_mass_initial, float &parton1_angle, float &parton2_angle, double parton1_spacetime_initial[], double parton2_spacetime_initial[], float parton1_final_momenta[], float parton2_final_momenta[], float parton1_final_mass, float parton2_final_mass, double parton1_final_spacetime[], double parton2_final_spacetime[]) 
{
	benergy1 = TMath::Sqrt(pow(parton1_momenta_initial[0],2) + pow(parton1_momenta_initial[1],2) + pow(parton1_momenta_initial[2],2) + pow(parton1_mass_initial,2));
	benergy2 = TMath::Sqrt(pow(parton2_momenta_initial[0],2) + pow(parton2_momenta_initial[1],2) + pow(parton2_momenta_initial[2],2) + pow(parton2_mass_initial,2));
	benergy3 = TMath::Sqrt(pow(parton1_final_momenta[0],2) + pow(parton1_final_momenta[1],2) + pow(parton1_final_momenta[2],2) + pow(parton1_final_mass,2));
	benergy4 = TMath::Sqrt(pow(parton2_final_momenta[0],2) + pow(parton2_final_momenta[1],2) + pow(parton2_final_momenta[2],2) + pow(parton2_final_mass,2));

	vb1.SetPxPyPzE(parton1_momenta_initial[0],parton1_momenta_initial[1],parton1_momenta_initial[2],benergy1);
	vb2.SetPxPyPzE(parton2_momenta_initial[0],parton2_momenta_initial[1],parton2_momenta_initial[2],benergy2);
	vb3.SetPxPyPzE(parton1_final_momenta[0],parton1_final_momenta[1],parton1_final_momenta[2],benergy3);
	vb4.SetPxPyPzE(parton2_final_momenta[0],parton2_final_momenta[1],parton2_final_momenta[2],benergy4);

	com = vb1 + vb2;

	vb1.Boost(-com.BoostVector());
	vb2.Boost(-com.BoostVector());
	vb3.Boost(-com.BoostVector());
	vb4.Boost(-com.BoostVector());


	parton1_angle = vb1.Angle(vb3.Vect());
	parton2_angle = vb2.Angle(vb4.Vect());
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
void draw(vector<float> x, vector<float> y, int iterate, vector<int> bigboom, int &actualevent) {

	// This is where all the canvases that will be used are created.
	TCanvas* c = new TCanvas(Form("c_%i_%i", iterate, eventnumber),Form("c_%i_%i", iterate, eventnumber),480,480);
	gStyle->SetOptStat(0);

	c->SetTickx();
	c->SetTicky();

	TH2F* hTimestep = new TH2F(Form("hTimestep_[%i]_[%i]", iterate, eventnumber), Form("hTimestep_%i_%i", iterate, eventnumber), 200, -2.5, 2.5, 200, -2.5, 2.5);
	hTimestep->SetTitle(Form("Event %i", eventnumber));
	hTimestep->GetXaxis()->SetTitle("x [fm]");
	hTimestep->GetYaxis()->SetTitle("y [fm]");
	hTimestep->Draw();

	// This loop is used to make the ellipse shape for the partons.
	for (int j = 0; j < x.size(); j++) {
		TEllipse *tell = new TEllipse(x[j], y[j],0.09,0.09);
		tell->SetFillColor(kBlue);
		tell->Draw("same");

		if (bigboom[j] > 0) {
			tell->SetFillColor(kRed);
			tell->Draw("same");
		}
	}

	// These if statements are simply to ensure that the images will be saved in order allowing them to be animated.
	if (iterate < 10) {
		c->SaveAs(Form("frame/Event_%i_00%i.png", eventnumber, iterate));
	}
	else if (iterate >= 10 && iterate < 100) {
		c->SaveAs(Form("frame/Event_%i_0%i.png", eventnumber, iterate));
	}
	else if (iterate >= 100) {
		c->SaveAs(Form("frame/Event_%i_%i.png", eventnumber, iterate));
	}

}

// This function is used to make a snapshot of all partons at their initial positions and show their momentum vectors.
void singleFrameDraw(int &actualevent, float &initialpx, float &initialpy, float &initialx, float &initialy, int &scatteringn, TCanvas* c1, int counter1, TH2F* Collisionplots) {

	c1->cd();

	if (counter1 == 1) {
		Collisionplots->Draw();
	}

	else if (counter1 > 1) {
		Collisionplots->Draw("same");
	}

	float phi = TMath::ATan2(initialpy,initialpx);
	float Tmomentum = TMath::Sqrt(pow(initialpx,2)+pow(initialpy,2));
	float x2 = initialx + (Tmomentum*TMath::Cos(phi));
	float y2 = initialy + (Tmomentum*TMath::Sin(phi));
	TArrow* ar1 = new TArrow(initialx,initialy,x2,y2,0.02,">");
	ar1->Draw();
	ar1->SetLineWidth(2);

	for (int i = 0; i < EventPartons.size(); i++) {

		TEllipse *tell2 = new TEllipse(initialx, initialy,0.09,0.09);

		if (scatteringn == 0) {
			tell2->SetFillColor(kBlue);
			tell2->Draw("same");
		}

		if (scatteringn > 0) {
			tell2->SetFillColor(kRed);
			tell2->Draw("same");
		}
	}

	if (counter1 == EventPartons.size()) {
		c1->SaveAs(Form("Single_frame/Event_%i.png", eventnumber));
	}
}

// This function will show what scattering angle applies to what parton.
void angleDraw(int &actualevent, float &initialpx, float &initialpy, float &initialpz, float &formationt, float &initialx, float &initialy, float &initialmass, int &scatteringn, float scatteringpx[], float scatteringpy[], float scatteringpz[], float scatteringt[], int counter1, vector<parton> v, TCanvas* c2, TH2F* Angleplots) {
	c2->cd();

	if (counter1 == 1) {
		Angleplots->Draw();
	}

	else if (counter1 > 1) {
		Angleplots->Draw("same");
	}

	float scatangle = 0;

	TEllipse *tell3 = new TEllipse(initialx, initialy,0.09,0.09);
	tell3->SetFillColor(kBlue);
	tell3->Draw("same");

	if (scatteringn > 0) {

		for (int j = 0; j < v.size(); j++) {

			if (j == 0) {

				scatangle = v[j].sangle;
			}

			if (j > 0) {

				scatangle = v[j].sangle;
			}

			if (v[j].sangle < (Pi/6)) {
				tell3->SetFillColor(kBlue-9);
				tell3->Draw("same");
			}

			if (v[j].sangle > (Pi/6) && scatangle < (Pi/3)) {
				tell3->SetFillColor(kGreen+3);
				tell3->Draw("same");
			}

			if (v[j].sangle > (Pi/3) && scatangle < (Pi/2)) {
				tell3->SetFillColor(kOrange-2);
				tell3->Draw("same");
			}

			if (v[j].sangle > (Pi/2)) {
				tell3->SetFillColor(kPink-9);
				tell3->Draw("same");
			}
		}
	}

	if (counter1 == EventPartons.size()) {
		c2->SaveAs(Form("Angle/AngleEvent_%i.png", eventnumber));
	}
}

// This function is used to do the position calculations and will return x,y coordinates at every time.
void calculatePosition (float actualtime, vector<parton> v, float &xt, float &yt, int &boom, int &actualevent, int &partonid, float &initialpx, float &initialpy, float &initialpz, float &formationt, float &initialx, float &initialy, float &initialz, float &initialmass, int &scatteringn, float scatteringpx[], float scatteringpy[], float scatteringpz[], float scatteringt[], float scatteringx[], float scatteringy[], float scatteringz[], TCanvas* c1, TCanvas* c2) {

	float xposition[v.size()];
	float yposition[v.size()];
	float xvelocity[v.size()];
	float yvelocity[v.size()];
	float x0, y0;

	// The v vector in this loop is looking at the scattering's for each individual parton so the size varies in length.
	for (int i = 0; i < v.size(); i++) {
		// Calculation of the velocity at each stage.
		float energy = TMath::Sqrt(pow(v[i].px,2) + pow(v[i].py,2) + pow(v[i].pz,2) + pow(v[i].m,2));
		TLorentzVector ev(v[i].px, v[i].py, v[i].pz, energy);
		float beta = ev.Beta();
		float phi = ev.Phi();
		xvelocity[i] = beta * TMath::Cos(phi);
		yvelocity[i] = beta * TMath::Sin(phi);

		scatteringn = v.size()-1;


		// Determine the initial position of the parton.
		if (i == 0) {

			x0 = v[i].x;
			y0 = v[i].y;

			xposition[i] = x0;
			yposition[i] = y0;

			// This is filling in the branch for initial conditions.
			actualevent = eventnumber;
			partonid = v[i].pID;
			initialpx = v[i].px;
			initialpy = v[i].py;
			initialpz = v[i].pz;
			formationt = v[i].t;
			initialx = x0;
			initialy = y0;
			initialz = v[i].z;
			initialmass = v[i].m;

		}
		else {
			//cout << "+++ " << stage << "   " << xvelocity[stage-1] << endl;

			x0 = xposition[i-1] + xvelocity[i-1] * (v[i].t - v[i-1].t);
			y0 = yposition[i-1] + yvelocity[i-1] * (v[i].t - v[i-1].t);

			xposition[i] = x0;
			yposition[i] = y0;

			// This is filling in the branches for scattering.
			scatteringpx[i-1] = v[i].px;
			scatteringpy[i-1] = v[i].py;
			scatteringpz[i-1] = v[i].pz;
			scatteringt[i-1] = v[i].t;
			scatteringx[i-1] = xposition[i];
			scatteringy[i-1] = yposition[i];
			scatteringz[i-1] = v[i].z;
		}
	}


	// Call getStage function to determine the stages for every time.
	int stage = getStage(actualtime,v);
	if (stage < 0) {
		xt =-999;
		yt=-999;
		return;
	}

	actualtime = actualtime - v[stage].t;

	// Calculating the positions as a function of time.
	xt = xposition[stage] + (xvelocity[stage] * actualtime);
	yt = yposition[stage] + (yvelocity[stage] * actualtime);
	boom = stage;

}

// This function will loop over each parton and then calculate the positions at every given moment in time. 
void processEvent(int &actualevent, int &partonid, float &initialpx, float &initialpy, float &initialpz, float &formationt, float &initialx, float &initialy, float &initialz, float &initialmass, int &scatteringn, float scatteringpx[], float scatteringpy[], float scatteringpz[], float scatteringt[], float scatteringx[], float scatteringy[], float scatteringz[],TCanvas* c1, TH2F* Collisionplots, TCanvas* c2, TH2F* Angleplots) {
	
	int iterate = 0;
	// This counter used to prevent the tree from being filled more than it should be.
	int counter1 = 0;

	// This loops through each time.
	for (int i = 0; i < NStep; i++) {

		float actualtime = i * dt;

		// Put Parton Loop here. At each timestep each parton is looped over. 
		
		for (int p = 0; p < EventPartons.size(); p++) {

			std::vector<parton> v = EventPartons[p];

			float xt, yt;
			int boom;
			calculatePosition(actualtime,v,xt,yt,boom,actualevent,partonid,initialpx,initialpy,initialpz,formationt,initialx,initialy,initialz,initialmass,scatteringn,scatteringpx,scatteringpy,scatteringpz,scatteringt,scatteringx,scatteringy,scatteringz,c1,c2);

			xvals.push_back(xt);
			yvals.push_back(yt);
			goesboom.push_back(boom);

			counter1 = counter1;
			counter1++;

			// This statement is to ensure that the tree is only filled once for all partons.
			if (counter1 <= EventPartons.size()) {
				tree->Fill();
			}
			else {
				continue;
			}

			angleDraw(actualevent,initialpx,initialpy,initialpz,formationt,initialx,initialy,initialmass,scatteringn,scatteringpx,scatteringpy,scatteringpz,scatteringt,counter1,v,c2,Angleplots);

			singleFrameDraw(actualevent,initialpx,initialpy,initialx,initialy,scatteringn,c1,counter1,Collisionplots);

		}

		//draw(xvals, yvals, iterate, goesboom, actualevent);


		xvals.clear();
		yvals.clear();
		goesboom.clear();
		iterate++;

	}

	EventPartons.clear();
}

// This will be my attempt to read in the files that I will be using. 
void EventAnimation(void) {

	f1 = new TFile("Multievent_AMPT_File_Tree.root", "RECREATE");

	Int_t scatteringn = 1000;
	Int_t actualevent;
	Int_t partonid;
	Float_t initialpx;
	Float_t initialpy;
	Float_t initialpz;
	Float_t formationt;
	Float_t initialx;
	Float_t initialy;
	Float_t initialz;
	Float_t initialmass;
	Float_t scatteringpx[scatteringn];
	Float_t scatteringpy[scatteringn];
	Float_t scatteringpz[scatteringn];
	Float_t scatteringt[scatteringn];
	Float_t scatteringx[scatteringn];
	Float_t scatteringy[scatteringn];
	Float_t scatteringz[scatteringn];

	tree = new TTree("tree", "An Orange Tree");
	tree->Branch("event_number",&actualevent,"event_number/I");
	tree->Branch("parton_id",&partonid,"parton_id/I");
	tree->Branch("initial_px",&initialpx,"initial_px/F");
	tree->Branch("initial_py",&initialpy,"initial_py/F");
	tree->Branch("initial_pz",&initialpz,"initial_pz/F");
	tree->Branch("formation_t",&formationt,"formation_t/F");
	tree->Branch("initial_x",&initialx,"initial_x/F");
	tree->Branch("initial_y",&initialy,"initial_y/F");
	tree->Branch("initial_z",&initialz,"initial_z/F");
	tree->Branch("initial_m",&initialmass,"initial_m/F");
	tree->Branch("scattering_n",&scatteringn,"scattering_n/I");
	tree->Branch("scattering_px",scatteringpx,"scattering_px[scattering_n]/F");
	tree->Branch("scattering_py",scatteringpy,"scattering_py[scattering_n]/F");
	tree->Branch("scattering_pz",scatteringpz,"scattering_pz[scattering_n]/F");
	tree->Branch("scattering_t",scatteringt,"scattering_t[scattering_n]/F");
	tree->Branch("scattering_x",scatteringx,"scattering_x[scattering_n]/F");
	tree->Branch("scattering_y",scatteringy,"scattering_y[scattering_n]/F");
	tree->Branch("scattering_z",scatteringz,"scattering_z[scattering_n]/F");
	ifstream myInitialFileInfo;
	ifstream myEvolutionFile;

	myInitialFileInfo.open("parton-initial-afterPropagation.dat");

	if (!myInitialFileInfo) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file parton-initial-afterPropagation.dat" << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened parton-initial-afterPropagation.dat" << endl;
	}


	while (myInitialFileInfo) {

		/*if (eventnumber % 100 == 0) {

			cout << "Reading Event Number " << eventnumber << endl;
		}*/

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
			partinfo.evtN = eventnumber;
			partinfo.pID = partID;
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
		float parton1_angle;
		float parton1_momenta_initial[3];
		float parton1_mass_initial;
		double parton1_spacetime_initial[4];
		int parton2_id_initial;
		float parton2_angle;
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

		float benergy1;
		float benergy2;
		float benergy3;
		float benergy4;

		TLorentzVector com;
		TLorentzVector vb1;
		TLorentzVector vb2;
		TLorentzVector vb3;
		TLorentzVector vb4;

		// This is where the evolution file is opened each time.
		myEvolutionFile.open("parton-collisionsHistory.dat");

		if (!myEvolutionFile) {
		// As before this will let me know if this specific file fails to open.
			cout << "Unable to open file parton-collisionsHistory.dat" << endl;
			return;
		}
		else {
		// This is here simply to confirm that this file opened successfully as well.
			cout << "Successfully opened parton-collisionsHistory.dat" << endl;
		}

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

					getAngle(benergy1,benergy2,benergy3,benergy4,vb1,vb2,vb3,vb4,com,parton1_momenta_initial,parton2_momenta_initial,parton1_mass_initial,parton2_mass_initial,parton1_angle,parton2_angle,parton1_spacetime_initial,parton2_spacetime_initial,parton1_final_momenta,parton2_final_momenta,parton1_final_mass,parton2_final_mass,parton1_final_spacetime,parton2_final_spacetime);

					parton part1;
					part1.evtN = evt;
					part1.pID = parton1_final_id;
					part1.px = parton1_final_momenta[0];
					part1.py = parton1_final_momenta[1];
					part1.pz = parton1_final_momenta[2];
					part1.sangle = parton1_angle;
					part1.x = parton1_final_spacetime[0];
					part1.y = parton1_final_spacetime[1];
					part1.z = parton1_final_spacetime[2];
					part1.t = parton1_final_spacetime[3];
					part1.m = parton1_final_mass;

					parton part2;
					part2.evtN = evt;
					part2.pID = parton2_final_id;
					part2.px = parton2_final_momenta[0];
					part2.py = parton2_final_momenta[1];
					part2.pz = parton2_final_momenta[2];
					part2.sangle = parton2_angle;
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

		// This is where the evolution file is closed each time.
		myEvolutionFile.close();

		// This is a new canvas used for the single frame images.
		TCanvas* c1 = new TCanvas(Form("c1_%i", eventnumber), Form("c1_%i", eventnumber),480,480);
		gStyle->SetOptStat(0);

		c1->SetTickx();
		c1->SetTicky();

		// Histograms used for the single frame images.
		TH2F* Collisionplots = new TH2F(Form("Collisionplots_%i", eventnumber), Form("Collisionplots_%i", eventnumber),200, -2.5, 2.5, 200, -2.5, 2.5);
		Collisionplots->SetTitle(Form("Event %i", eventnumber));
		Collisionplots->GetXaxis()->SetTitle("x [fm]");
		Collisionplots->GetYaxis()->SetTitle("y [fm]");

		// New Canvas for scattering plots.
		TCanvas* c2 = new TCanvas(Form("c2_%i",eventnumber), Form("c2_%i", eventnumber),480,480);
		gStyle->SetOptStat(0);
		c2->SetTickx();
		c2->SetTicky();

		// Histogram for scattering plots
		TH2F* Angleplots = new TH2F(Form("Angleplots_%i", eventnumber), Form("Angleplots_%i", eventnumber),200,-2.5,2.5,200,-2.5,2.5);
		Angleplots->SetTitle(Form("Event %i", eventnumber));
		Angleplots->GetXaxis()->SetTitle("x [fm]");
		Angleplots->GetYaxis()->SetTitle("y [fm]");

		// Call function that Processes the Event.
		processEvent(actualevent,partonid,initialpx,initialpy,initialpz,formationt,initialx,initialy,initialz,initialmass,scatteringn,scatteringpx,scatteringpy,scatteringpz,scatteringt,scatteringx,scatteringy,scatteringz,c1,Collisionplots,c2,Angleplots);

	}

	f1->Write();
	f1->Close();

	myInitialFileInfo.close();
	return;
}

