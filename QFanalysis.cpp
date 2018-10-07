#define QFanalysis_cxx
#include "QFanalysis.h"
#include <stdlib.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Riostream.h"

using namespace std;

void QFanalysis::Loop()
{
	//   In a ROOT session, you can do:
	//      Root > .L QFanalysis.C
	//      Root > QFanalysis t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	//char filename[14] = "small_all.root";
	char basic[13];
	char LBAname[21];
	char LBCname[21];	char EBCname[21];

	float ETA[4][64][48];	
	float PHI[4][64][48];
	float Rmeters[4][64][48];
	int partner[4][48];

	//Made by Joao Pedro!
	gStyle->SetOptStat(0);
	int QFcut;
	QFcut = 100;
	TString SQFcut;
	SQFcut = Form("%i",QFcut);

	QFcut_OF2 = 20;
	TString SQFcut_OF2;
	SQFcut_OF2 = Form("%i",QFcut_OF2);

	cout << "QF < " << SQFcut << endl;

	int Ratio_cell_all, Ratio_cell_0to05, Ratio_cell_05to15, Ratio_cell_15to3, Ratio_cell_3to10, Ratio_cell_justCOF, Ratio_cell_justOF2;
	int Ratio_tower_all, Ratio_tower_0to05, Ratio_tower_05to15, Ratio_tower_15to3, Ratio_tower_3to10, Ratio_tower_justCOF, Ratio_tower_justOF2;
	int all_rec_cell_ene, all_Edif0_cell_ene, all_rec_cell_eMF, all_Edif0_cell_eMF, pile_ev;
	int pileup, pileup_event[37821], k;
	int just1cell_ene, just1cell_eMF, more1cell_ene, more1cell_eMF;
	int i , j, p, q, u, v, zp, g, n, zmu1_ene, zmu1_eMF, count, events_rem, corresponding_jet, zr1_dif, zr_dif, m;
	int pile01, pile02, pile03, pile12, pile13, pile23, pile43, pile53, pile63, pile54, pile64, pile65, total, only_COF, only_OF2;
	int n_ene, m_ene, k_ene, zi_ene, zii_ene, zr_ene, zmu_ene1, zmu_ene2, zmu_ene;
	int n_eMF, m_eMF, k_eMF, zi_eMF, zii_eMF, zr_eMF, zmu_eMF1, zmu_eMF2, zmu_eMF;
	Double_t jentry;
	int k_ene_bigger, k_eMF_bigger, k_ene_eMF_equal, muoncandidate_ene_bigger, muoncandidate_ene_eMF_equal, muoncandidate_eMF_bigger;
	int pileupcell_ene, pileupcell_eMF, cell_ene, cell_eMF;
	double Dpileupcell_ene, Dpileupcell_eMF, Dcell_ene, Dcell_eMF;
	int minus0_eMF, I0to100_eMF, I100to200_eMF, I200to300_eMF;
	int minus0_ene, I0to100_ene, I100to200_ene, I200to300_ene;
	int Eta_photon1, Eta_photon2, Phi_photon1, Phi_photon2;
	int njet_ene[1000], mjet_ene[1000], zjet_ene[1000], zr, cellpile, jetpile, jetmatch;
	int Zr_ene[1000], Zr_eMF[1000];
	int njet_eMF[1000], mjet_eMF[1000], zjet_eMF[1000];
	int pile_z1_ene = 0, jet_pile_ene = 0, mucand_pile_ene = 0, pile_z1_eMF = 0, jet_pile_eMF = 0, mucand_pile_eMF = 0, ismucand_ene[1000], ismucand_eMF[1000], all_pile_ene = 0, all_pile_eMF;
	int pile_z1_al_ene = 0, jet_pile_al_ene = 0, mucand_pile_al_ene = 0, pile_z1_al_eMF = 0, jet_pile_al_eMF = 0, mucand_pile_al_eMF = 0;
	float jet_ene[1000], etacm_ene[1000], phicm_ene[1000]; 
	float jet_eMF[1000], etacm_eMF[1000], phicm_eMF[1000];
	float Emu1_ene[1000], Emu2_ene[1000], Emu3_ene[1000];
	float Emu1_eMF[1000], Emu2_eMF[1000], Emu3_eMF[1000];
	float deltaR_ene, deltaR_eMF, phimiss, M_particle, ang_var;
	float theta_ene, px_ene, py_ene, pz_ene;
	float theta_eMF, px_eMF, py_eMF, pz_eMF;
	float theta_ene2, px_ene2, py_ene2, pz_ene2;
	float theta_eMF2, px_eMF2, py_eMF2, pz_eMF2;
	float Ejetsmu1_ene[100][4], Ejetsmu1_eMF[100][4];
	bool tfbin_ene[40][64], tfjet_ene[1000], zprint, done;
	bool tfbin_eMF[40][64], tfjet_eMF[1000], liquid;
	int OF2COFall, OF2COFpileup, OF2COFup1down2, OF2COFup2down3, OF2COFup3, OF2COFup0down1,  OF2COFup0down05, OF2COFup05down1, OF2COFup_1down0, OF2COFup_2down_1, OF2COFup_3down_2, OF2COFdown_3, OF2COFup05down15, n_channel;
	int muoncandidate_ene, muoncandidate_eMF, pile, pileevent, dif_muon, dif_jets, dif_z1, dif_zi, dif_zii, nopile;
	//noise cut
	float noisecut_eMF = 300, noisecut_ene = 300, b, r, t, pT_miss_ene1, pT_miss_ene2, eta_miss, theta_miss, sum, dif_jet, la, lbc, ld;
	float phi_miss_ene1, phi_miss_eMF1, pT_miss_eMF1, pT_miss_eMF2, phi_miss_ene2, phi_miss_eMF2, theta, ra, rbc, rd;
	float Ephoton1, Ephoton2, theta_particle, Ephoton;

	Ephotoncut=2000;
	intEphoton=2000;
	TString SEphoton;
	SEphoton = Form("%i",intEphoton);
	//SEphoton = "nocut";
	TString layerpoint;
	layerpoint = "end";

	/*
	//If we consider the middle of the layer, we will have for R/r:

		Double_t Rrlonglayer1 = 15.25;
		Double_t Rrlonglayer2 = 6.36;
		Double_t Rrlonglayer3 = 13;
		
		Double_t Rrmediumlayer1 = 15.25;
		Double_t Rrmediumlayer2 = 6.36;
		Double_t Rrmediumlayer3 = 7.61;

		Double_t Rrextendedlayer1 = 15.25;
		Double_t Rrextendedlayer2 = 9.44;
		Double_t Rrextendedlayer3 = 7.61;
	*/
	//If we consider the end of the layer, we will have for R/r:
	
		Double_t Rrlonglayer1 = 8.125;
		Double_t Rrlonglayer2 = 3.68;
		Double_t Rrlonglayer3 = 7;
		
		Double_t Rrmediumlayer1 = 8.125; 
		Double_t Rrmediumlayer2 = 3.68;
		Double_t Rrmediumlayer3 = 4.30;

		Double_t Rrextendedlayer1 = 8.125;
		Double_t Rrextendedlayer2 = 5.22;
		Double_t Rrextendedlayer3 = 4.30;
	
	char a[10];
	int intnoisecut_eMF = noisecut_eMF;
	TString Snoisecut_eMF;
	Snoisecut_eMF = Form("%i",intnoisecut_eMF);
	gStyle->SetTitleOffset(1.4,"X");	
	gStyle->SetTitleOffset(1.4,"Y");
	gStyle->SetTitleOffset(1.6,"Z");

	//Histograms of eta x phi distribuition of energy in the event 

	//All layers
	TH2F *Eta_Phi_E_ene = new TH2F("Eta_Phi_E_ene","Energy distribution in Eta x Phi OF2 QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E_eMF = new TH2F("Eta_Phi_E_eMF","Energy distribution in Eta x Phi COF QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	//Layer 1
	TH2F *Eta_Phi_E_ene1 = new TH2F("Eta_Phi_E_ene1","Energy distribution in Eta x Phi OF2 1st layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E_eMF1 = new TH2F("Eta_Phi_E_eMF1","Energy distribution in Eta x Phi COF 1st layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	//Layer 2
	TH2F *Eta_Phi_E_ene2 = new TH2F("Eta_Phi_E_ene2","Energy distribution in Eta x Phi OF2 2nd layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E_eMF2 = new TH2F("Eta_Phi_E_eMF2","Energy distribution in Eta x Phi COF 2nd layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	//Layer 3
	TH2F *Eta_Phi_E_ene3 = new TH2F("Eta_Phi_E_ene3","Energy distribution in Eta x Phi OF2 3rd layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E_eMF3 = new TH2F("Eta_Phi_E_eMF3","Energy distribution in Eta x Phi COF 3rd layer QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);

	//Pileup
	TH2F *Eta_Phi_E_pileup_ene = new TH2F("Eta_Phi_E_pileup_ene","Energy distribution in Eta x Phi OF2 QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E_pileup_eMF = new TH2F("Eta_Phi_E_pileup_eMF","Energy distribution in Eta x Phi COF QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);

	TH2F *Eta_Phi_E2_ene = new TH2F("Eta_Phi_E2_ene","Energy distribution in Eta x Phi OF2 QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);
	TH2F *Eta_Phi_E2_eMF = new TH2F("Eta_Phi_E2_eMF","Energy distribution in Eta x Phi COF QF < " + SQFcut, 40, -2, 2, 64, 0, 6.4);

	//Setting the X and Y axis' name
	Eta_Phi_E_ene->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_ene->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_ene->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_eMF->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_eMF->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_eMF->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_ene->GetXaxis()->CenterTitle();
	Eta_Phi_E_ene->GetYaxis()->CenterTitle();
	Eta_Phi_E_eMF->GetXaxis()->CenterTitle();
	Eta_Phi_E_eMF->GetYaxis()->CenterTitle();
	Eta_Phi_E_ene1->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_ene1->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_ene1->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_eMF1->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_eMF1->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_eMF1->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_ene1->GetXaxis()->CenterTitle();
	Eta_Phi_E_ene1->GetYaxis()->CenterTitle();
	Eta_Phi_E_eMF1->GetXaxis()->CenterTitle();
	Eta_Phi_E_eMF1->GetYaxis()->CenterTitle();
	Eta_Phi_E_ene2->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_ene2->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_ene2->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_eMF2->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_eMF2->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_eMF2->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_ene2->GetXaxis()->CenterTitle();
	Eta_Phi_E_ene2->GetYaxis()->CenterTitle();
	Eta_Phi_E_eMF2->GetXaxis()->CenterTitle();
	Eta_Phi_E_eMF2->GetYaxis()->CenterTitle();
	Eta_Phi_E_ene3->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_ene3->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_ene3->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_eMF3->GetXaxis()->SetTitle("Eta");
	Eta_Phi_E_eMF3->GetYaxis()->SetTitle("Phi");
	Eta_Phi_E_eMF3->GetZaxis()->SetTitle("E (MeV)");
	Eta_Phi_E_ene3->GetXaxis()->CenterTitle();
	Eta_Phi_E_ene3->GetYaxis()->CenterTitle();
	Eta_Phi_E_eMF3->GetXaxis()->CenterTitle();
	Eta_Phi_E_eMF3->GetYaxis()->CenterTitle();

	//Creating the histograms

	//Energy and number of jets
	//OF2
	TH1I *Njets_per_event_ene = new TH1I("Njets_per_event_ene","OF2 - Number of jets per event (QF<" + SQFcut + ")",40,0,40);
	TH1I *Njets_per_event_z1_ene = new TH1I("Njets_per_event_z1_ene","OF2 - Number of jets per event (Zr>1, QF<" + SQFcut + ")",40, 0, 40);
	TH1I *Njets_per_event_zequal1_ene = new TH1I("Njets_per_event_zequal1_ene","OF2 - Number of jets per event (Zr=1, QF<" + SQFcut + ")",40, 0, 40);
	TH1I *Njets_per_event_z2_ene = new TH1I("Njets_per_event_z2_ene","OF2 - Number of jets per event (Zr>2, QF<" + SQFcut + ")",40, 0, 40);

	TH1F *Ejets_ene = new TH1F("Ejets_ene","OF2 - Energy of the jets of all events QF<" + SQFcut + ")",100,0,300000);
	TH1F *Ejets_z1_ene = new TH1F("Ejets_z1_ene","OF2 - Energy of the jets (Zr>1, QF<" + SQFcut + ")",100,0,300000);
	TH1F *Ejets_zequal1_ene = new TH1F("Ejets_zequal1_ene","OF2 - Energy of the jets (Zr=1, QF<" + SQFcut + ")",100,0,100000);
	TH1F *Ejets_z2_ene = new TH1F("Ejets_z2_ene","OF2 - Energy of the jets (Zr>2, QF<" + SQFcut + ")",100,0,300000);

	//COF
	TH1I *Njets_per_event_eMF = new TH1I("Njets_per_event_eMF","COF - Number of jets per event (QF<" + SQFcut + ")",40,0,40);
	TH1I *Njets_per_event_z1_eMF = new TH1I("Njets_per_event_z1_eMF","COF - Number of jets per event (Zr>1 QF<" + SQFcut + ")",40, 0, 40);
	TH1I *Njets_per_event_zequal1_eMF = new TH1I("Njets_per_event_zequal1_eMF","COF - Number of jets per event (Zr=1 QF<" + SQFcut + ")",40, 0, 40);
	TH1I *Njets_per_event_z2_eMF = new TH1I("Njets_per_event_z2_eMF","COF - Number of jets per event (Zr>2 QF<" + SQFcut + ")",40, 0, 40);

	TH1F *Ejets_eMF = new TH1F("Ejets_eMF","COF - Energy of the jets of all events (QF<" + SQFcut + ")",100,0,300000);
	TH1F *Ejets_z1_eMF = new TH1F("Ejets_z1_eMF","COF - Energy of the jets (Zr>1, QF<" + SQFcut + ")",100,0,300000);
	TH1F *Ejets_zequal1_eMF = new TH1F("Ejets_zequal1_eMF","COF - Energy of the jets (Zr=1, QF<" + SQFcut + ")",100,0,100000);
	TH1F *Ejets_z2_eMF = new TH1F("Ejets_z2_eMF","COF - Energy of the jets (Zr>2, QF<" + SQFcut + ")",100,0,300000);

	//OF2xCOF (just number of jets)
	TH2F *Njets_per_event_ene_eMF = new TH2F("Njets_per_event_ene_eMF","Number of jets per event (QF<" + SQFcut + ")",40,0,40,40,0,40);
	TH2F *Njets_per_event_z1_ene_eMF = new TH2F("Njets_per_event_z1_ene_eMF","Number of jets per event (Zr>1, QF<" + SQFcut + ")",40,0,40,40,0,40);
	TH2F *Njets_per_event_zequal1_ene_eMF = new TH2F("Njets_per_event_zequal1_ene_eMF","Number of jets per event (Zr=1, QF<" + SQFcut + ")",40,0,40,40,0,40);
	TH2F *Njets_per_event_z2_ene_eMF = new TH2F("Njets_per_event_z2_ene_eMF","Number of jets per event (Zr>2, QF<" + SQFcut + ")",40,0,40,40,0,40);
	TH2F *Nmu_per_event_z1_ene_eMF1 = new TH2F("Nmu_per_event_z1_ene_eMF1","Number of muon candidates per event (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")",40,0,40,40,0,40);

	//Emiss histograms OF2 and COF
	TH1F *PT_miss_ene = new TH1F("PT_miss_ene","OF2 - pT miss QF<" + SQFcut,100,0,100000);
	TH1F *PT_miss_eMF = new TH1F("PT_miss_eMF","COF - pT miss QF<" + SQFcut,100,0,100000);
	TH1F *Phi_miss_ene = new TH1F("Phi_miss_ene","OF2 - phi miss QF<" + SQFcut,64,0,6.4);
	TH1F *Phi_miss_eMF = new TH1F("Phi_miss_eMF","COF - phi miss QF<" + SQFcut,64,0,6.4);
	TH1F *eta_miss_ene = new TH1F("eta_miss_ene","OF2 - eta miss QF<" + SQFcut,40,-2.,2.);
	TH1F *eta_miss_eMF = new TH1F("eta_miss_eMF","COF - eta miss QF<" + SQFcut,40,-2.,2.);

	//Emiss histograms OF2 and COF channel by channel
	TH1F *PT_miss_ene2 = new TH1F("PT_miss_ene2","OF2 - pT miss QF<" + SQFcut,100,0,100000);
	TH1F *PT_miss_eMF2 = new TH1F("PT_miss_eMF2","COF - pT miss QF<" + SQFcut,100,0,100000);
	TH1F *PT_miss_pileev_ene2 = new TH1F("PT_miss_pileev_ene2","OF2 - pT miss of pileup events QF<" + SQFcut,100,0,100000);
	TH1F *PT_miss_pileev_eMF2 = new TH1F("PT_miss_pileev_eMF2","COF - pT miss of pileup events QF<" + SQFcut,100,0,100000);
	TH1F *Phi_miss_ene2 = new TH1F("Phi_miss_ene2","OF2 - phi miss QF<" + SQFcut,64,0,6.4);
	TH1F *Phi_miss_eMF2 = new TH1F("Phi_miss_eMF2","COF - phi miss QF<" + SQFcut,64,0,6.4);
	TH1F *eta_miss_ene2 = new TH1F("eta_miss_ene2","OF2 - eta miss QF<" + SQFcut,40,-2.,2.);
	TH1F *eta_miss_eMF2 = new TH1F("eta_miss_eMF2","COF - eta miss QF<" + SQFcut,40,-2.,2.);

	//eta_cm and phi_cm histograms
	TH1F *phi_cm_ene = new TH1F("phi_cm_ene","OF2 - phi_jet QF<" + SQFcut,32,0,6.4);
	TH1F *eta_cm_ene = new TH1F("eta_cm_ene","OF2 - eta_jet QF<" + SQFcut,40,-2.,2.);
	TH1F *phi_cm_eMF = new TH1F("phi_cm_eMF","COF - phi_jet QF<" + SQFcut,32,0,6.4);
	TH1F *eta_cm_eMF = new TH1F("eta_cm_eMF","COF - eta_jet QF<" + SQFcut,40,-2.,2.);

	//TProfile
	EnexQF  = new TProfile("EnexQF","TProfile - OF2 - Energy versus QF", 75, 0, 30000, 0, 100);
	EMFxQF  = new TProfile("EMFxQF","TProfile - COF - Energy versus QF", 75, 0, 30000, 0, 100);

	TH2F *EnexQF2 = new TH2F("EnexQF2","OF2 - Energy versus QF", 75, 0, 30000, 100, 0, 50);
	TH2F *EMFxQF2 = new TH2F("EMFxQF2","COF - Energy versus QF", 75, 0, 30000, 100, 0, 50);

	EnexQF3  = new TProfile("EnexQF3","TProfile - OF2 - QF versus Energy",75,0,100,0,30000);
	EMFxQF3  = new TProfile("EMFxQF3","TProfile - COF - QF versus Energy",75,0,100,0,30000);

	EneXQF_100 = new TProfile("EnexQF_100","TProfile - OF2 - Energy versus 0 < E < 100 GeV and 0.01 < QF < 50",20,0,100000,0,50);
	EMFXQF_100 = new TProfile("EMFxQF_100","TProfile - COF - Energy versus 0 < E < 100 GeV and 0.01 < QF < 50",20,0,100000,0,50);

	EneXQF_50 = new TProfile("EnexQF_50","TProfile - OF2 - Energy versus 0 < E < 50 GeV and 0.01 < QF < 50",20,0,50000,0,100);
	EMFXQF_50 = new TProfile("EMFxQF_50","TProfile - COF - Energy versus 0 < E < 50 GeV and 0.01 < QF < 50",20,0,50000,0,100);

	EneXQF_30 = new TProfile("EnexQF_30","TProfile - OF2 - Energy versus 0 < E < 30 GeV and 0.01 < QF < 50",20,0,30000,0,100);
	EMFXQF_30 = new TProfile("EMFxQF_30","TProfile - COF - Energy versus 0 < E < 30 GeV and 0.01 < QF < 50",20,0,30000,0,100);

	//Check the layers
	TH1F *Ejets_all3layers_ene = new TH1F("Ejets_all3layers_ene","OF2 - Jets with energy in all layers (Emu1/Emu2 < 1.5 and Emu2/Emu3 < 1.5) (z=1 QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_all3layers_eMF = new TH1F("Ejets_all3layers_eMF","COF - Jets with energy in all layers (Emu1/Emu2 < 1.5 and Emu2/Emu3 < 1.5) (z=1 QF<"+SQFcut+")",100,0,15000);

	TH1F *Ejets_layers_ene = new TH1F("Ejets_layers_ene","OF2 - Single cell jet all TileCal layers (Zr=1 QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_layers_eMF = new TH1F("Ejets_layers_eMF","COF - Single cell jet all TileCal layers (Zr=1 QF<"+SQFcut+")",100,0,15000);

	TH1F *Emu_layer1_ene = new TH1F("Emu_layer1_ene","OF2 - Energy of muons in the first layer QF<"+SQFcut+")",100,0,10000);
	TH1F *Emu_layer2_ene = new TH1F("Emu_layer2_ene","OF2 - Energy of muons in the second layer QF<"+SQFcut+")",100,0,10000);
	TH1F *Emu_layer3_ene = new TH1F("Emu_layer3_ene","OF2 - Energy of muons in the third layer QF<"+SQFcut+")",100,0,10000);

	TH1F *Emu_layer1_eMF = new TH1F("Emu_layer1_eMF","COF - Energy of muons in the first layer QF<"+SQFcut+")",100,0,10000);
	TH1F *Emu_layer2_eMF = new TH1F("Emu_layer2_eMF","COF - Energy of muons in the second layer QF<"+SQFcut+")",100,0,10000);
	TH1F *Emu_layer3_eMF = new TH1F("Emu_layer3_eMF","COF - Energy of muons in the third layer QF<"+SQFcut+")",100,0,10000);

	//Zr distribution
	TH1I *Zr_dist_ene = new TH1I("Zr_dist_ene","OF2 - Zr distribuition QF<"+SQFcut,20,0,20);
	TH1I *Zr_dist_eMF = new TH1I("Zr_dist_eMF","COF - Zr distribuition QF<"+SQFcut,20,0,20);

	//DeltaR distribution
	TH1F *DeltaR_dist_ene = new TH1F("DeltaR_dist_ene","OF2 - Delta R distribution QF<"+SQFcut,18,0,2.6);
	TH1F *DeltaR_dist_eMF = new TH1F("DeltaR_dist_eMF","COF - Delta R distribution QF<"+SQFcut,18,0,2.6);

	//----------------AXION PARTICLE SEARCH
	TString Sang_var;
	ang_var=0.25;
	Sang_var = Form("%f",ang_var);

	//Energy of layers 
	TH1F *E_nomuon_l1_ene = new TH1F("E_nomuon_l1_ene","OF2 - Energy of no-muons in the first layer QF<"+SQFcut+")",100,0,3000);
	TH1F *E_nomuon_l1_eMF = new TH1F("E_nomuon_l1_eMF","COF - Energy of no-muons in the first layer QF<"+SQFcut+")",100,0,3000);
	TH1F *E_nomuon_l2_ene = new TH1F("E_nomuon_l2_ene","OF2 - Energy of no-muons in the second layer QF<"+SQFcut+")",100,0,3000);
	TH1F *E_nomuon_l2_eMF = new TH1F("E_nomuon_l2_eMF","COF - Energy of no-muons in the second layer QF<"+SQFcut+")",100,0,3000);
	TH1F *E_nomuon_l3_ene = new TH1F("E_nomuon_l3_ene","OF2 - Energy of no-muons in the third layer QF<"+SQFcut+")",100,0,3000);
	TH1F *E_nomuon_l3_eMF = new TH1F("E_nomuon_l3_eMF","COF - Energy of no-muons in the third layer QF<"+SQFcut+")",100,0,3000);

	//Energy/mm of each layer
	TH1F *Emm_nomuon_l1_ene = new TH1F("Emm_nomuon_l1_ene","OF2 - Energy per cm of no-muons in the first layer QF<"+SQFcut+")",100,0,50);
	TH1F *Emm_nomuon_l1_eMF = new TH1F("Emm_nomuon_l1_eMF","COF - Energy per cm of no-muons in the first layer QF<"+SQFcut+")",100,0,50);
	TH1F *Emm_nomuon_l2_ene = new TH1F("Emm_nomuon_l2_ene","OF2 - Energy per cm of no-muons in the second layer QF<"+SQFcut+")",100,0,50);
	TH1F *Emm_nomuon_l2_eMF = new TH1F("Emm_nomuon_l2_eMF","COF - Energy per cm of no-muons in the second layer QF<"+SQFcut+")",100,0,50);
	TH1F *Emm_nomuon_l3_ene = new TH1F("Emm_nomuon_l3_ene","OF2 - Energy per cm of no-muons in the third layer QF<"+SQFcut+")",100,0,50);
	TH1F *Emm_nomuon_l3_eMF = new TH1F("Emm_nomuon_l3_eMF","COF - Energy per cm of no-muons in the third layer QF<"+SQFcut+")",100,0,50);

	TH2F *Emm_nomuon_l12_ene = new TH2F("Emm_nomuon_l12_ene","OF2 - Energy per cm of no-muons in layer 1 versus layer 2 QF<"+SQFcut+")",100,0,50,100,0,50);
	TH2F *Emm_nomuon_l12_eMF = new TH2F("Emm_nomuon_l12_eMF","COF - Energy per cm of no-muons in layer 1 versus layer 2 QF<"+SQFcut+")",100,0,50,100,0,50);
	TH2F *Emm_nomuon_l13_ene = new TH2F("Emm_nomuon_l13_ene","OF2 - Energy per cm of no-muons in layer 1 versus layer 3 QF<"+SQFcut+")",100,0,50,100,0,50);
	TH2F *Emm_nomuon_l13_eMF = new TH2F("Emm_nomuon_l13_eMF","COF - Energy per cm of no-muons in layer 1 versus layer 3 QF<"+SQFcut+")",100,0,50,100,0,50);
	TH2F *Emm_nomuon_l23_ene = new TH2F("Emm_nomuon_l23_ene","OF2 - Energy per cm of no-muons in layer 2 versus layer 3 QF<"+SQFcut+")",100,0,50,100,0,50);
	TH2F *Emm_nomuon_l23_eMF = new TH2F("Emm_nomuon_l23_eMF","COF - Energy per cm of no-muons in layer 2 versus layer 3 QF<"+SQFcut+")",100,0,50,100,0,50);

	TH1F *Emm_all_nocut_ene = new TH1F("Emm_all_nocut_ene","OF2 - Energy per cm QF<"+SQFcut+")",100,0,25);
	TH1F *Emm_all_nocut_eMF = new TH1F("Emm_all_nocut_eMF","COF - Energy per cm QF<"+SQFcut+")",100,0,25);

	//Energy just on layer 1, 2, 3, 1 and 2 or 2 and 3
	TH1F *Ejets_2layers_ene = new TH1F("Ejets_2layers_ene","OF2 - Jets with energy just in the second layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_2layers_eMF = new TH1F("Ejets_2layers_eMF","COF - Jets with energy just in the second layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);

	TH1F *Ejets_3layers_ene = new TH1F("Ejets_3layers_ene","OF2 - Jets with energy just in the third layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_3layers_eMF = new TH1F("Ejets_3layers_eMF","COF - Jets with energy just in the third layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);

	TH1F *Ejets_1layers_ene = new TH1F("Ejets_1layers_ene","OF2 - Jets with energy just in the first layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_1layers_eMF = new TH1F("Ejets_1layers_eMF","COF - Jets with energy just in the first layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);

	TH1F *Ejets_12layers_ene = new TH1F("Ejets_12layers_ene","OF2 - Jets with energy just in the first and second layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_12layers_eMF = new TH1F("Ejets_12layers_eMF","COF - Jets with energy just in the first and second layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);

	TH1F *Ejets_23layers_ene = new TH1F("Ejets_23layers_ene","OF2 - Jets with energy just in the second and third layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);
	TH1F *Ejets_23layers_eMF = new TH1F("Ejets_23layers_eMF","COF - Jets with energy just in the second and third layer (Zr=1 or Zr=1+1cell QF<"+SQFcut+")",100,0,15000);

	TH1F *N_channel_per_jet_ene = new TH1F("N_channel_per_jet_ene", "OF2 - Number of channels with energy per jet (E just layer 1) (Zr=1 or Zr=1+1cell) QF<"+SQFcut+")",10,0,10);
	TH1F *N_channel_per_jet_eMF = new TH1F("N_channel_per_jet_eMF", "COF - Number of channels with energy per jet (E just layer 1) (Zr=1 or Zr=1+1cell) QF<"+SQFcut+")",10,0,10);
	TH2F *N_channel_per_jet_Energy_ene = new TH2F("N_channel_per_jet_Energy_ene", "OF2 - Number of channels with energy per jet versus energy (E just layer 1) (Zr=1 or Zr=1+1cell) QF<"+SQFcut+")",10,0,10,100,0,15000);
	TH2F *N_channel_per_jet_Energy_eMF = new TH2F("N_channel_per_jet_Energy_eMF", "COF - Number of channels with energy per jet versus energy (E just layer 1) (Zr=1 or Zr=1+1cell) QF<"+SQFcut+")",10,0,10,100,0,15000);

	//Theta particle chute
	TH1F *THETA_Particle_ene = new TH1F ("THETA_Particle_ene", "OF2 - Theta particle (E just layer 1 or 3, 2 channels) QF<"+SQFcut+")",32,0,3.2);
	TH1F *THETA_Particle_eMF = new TH1F ("THETA_Particle_eMF", "COF - Theta particle (E just layer 1 or 3, 2 channels) QF<"+SQFcut+")",32,0,3.2);

	//Mass particle chute
	TH1F *MASS_Particle_ene = new TH1F ("MASS_Particle_ene", "OF2 - Invariant mass gamma+gamma of Zr = 1+1 QF<"+SQFcut+")",100, 0, 3000); // PEDRO!
	TH1F *MASS_Particle_eMF = new TH1F ("MASS_Particle_eMF", "COF - Invariant mass gamma+gamma of Zr = 1+1 QF<"+SQFcut+")",100, 0, 3000); // PEDRO!
	TH1F *MASS_Particle_l1_ene = new TH1F ("MASS_Particle_l1_ene", "OF2 - Invariant mass gamma+gamma of Zr = 1+1 L1 QF<"+SQFcut+")",100, 0, 3000); 
	TH1F *MASS_Particle_l1_eMF = new TH1F ("MASS_Particle_l1_eMF", "COF - Invariant mass gamma+gamma of Zr = 1+1 L1 QF<"+SQFcut+")",100, 0, 3000); 
	TH1F *MASS_Particle_l2_ene = new TH1F ("MASS_Particle_l2_ene", "OF2 - Invariant mass gamma+gamma of Zr = 1+1 L2 QF<"+SQFcut+")",100, 0, 3000); 
	TH1F *MASS_Particle_l2_eMF = new TH1F ("MASS_Particle_l2_eMF", "COF - Invariant mass gamma+gamma of Zr = 1+1 L2 QF<"+SQFcut+")",100, 0, 3000); 
	TH1F *MASS_Particle_l3_ene = new TH1F ("MASS_Particle_l3_ene", "OF2 - Invariant mass gamma+gamma of Zr = 1+1 L3 QF<"+SQFcut+")",100, 0, 3000); 
	TH1F *MASS_Particle_l3_eMF = new TH1F ("MASS_Particle_l3_eMF", "COF - Invariant mass gamma+gamma of Zr = 1+1 L3 QF<"+SQFcut+")",100, 0, 3000); 

	//----------------------------------------

	//Unused---

	//Energy in the layers
	TH2F *Elayers12_ene = new TH2F("Elayers12_ene","OF2 - Jets with energy only in first and second layers (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers12_eMF = new TH2F("Elayers12_eMF","COF - Jets with energy only in first and second layers (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH2F *Elayers13_ene = new TH2F("Elayers13_ene","OF2 - Jets with energy only in first and third layers (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers13_eMF = new TH2F("Elayers13_eMF","COF - Jets with energy only in first and third layers (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH2F *Elayers23_ene = new TH2F("Elayers23_ene","OF2 - Jets with energy only in second and third layers (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers23_eMF = new TH2F("Elayers23_eMF","COF - Jets with energy only in second and third layers (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH2F *Elayers1_ene = new TH2F("Elayers1_ene","OF2 - Jets with energy only in first layer (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers1_eMF = new TH2F("Elayers1_eMF","COF - Jets with energy only in first layer (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH2F *Elayers2_ene = new TH2F("Elayers2_ene","OF2 - Jets with energy only in second layer (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers2_eMF = new TH2F("Elayers2_eMF","COF - Jets with energy only in second layer (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH2F *Elayers3_ene = new TH2F("Elayers3_ene","OF2 - Jets with energy only in third layer (QF<"+SQFcut+")",100,0,15000,16,0,16);
	TH2F *Elayers3_eMF = new TH2F("Elayers3_eMF","COF - Jets with energy only in third layer (QF<"+SQFcut+")",100,0,15000,16,0,16);

	TH1F *Elayersall_ene = new TH1F("Elayersall_ene","OF2 - Jets with energy in all layers (Zr = 1 QF<"+SQFcut+")", 100, 0, 15000);
	TH1F *Elayersall_eMF = new TH1F("Elayersall_eMF","OF2 - Jets with energy in all layers (Zr = 1 QF<"+SQFcut+")", 100, 0, 15000);

	//-----------


	TH1F *Ratio_E_l_ene = new TH1F("Ratio_E_l_ene","OF2 - Energy/Length in each layer (Zr = 1 QF<"+SQFcut+")", 100, 0, 100);
	TH1F *Ratio_E_l_eMF = new TH1F("Ratio_E_l_eMF","COF - Energy/Length in each layer (Zr = 1 QF<"+SQFcut+")", 100, 0, 100);

	//phi_miss comparison
	TH2F *Phi_miss_ene_eMF = new TH2F("Phi_miss_ene_eMF", "phi miss OF2 versus COF", 64,0,6.4, 64,0,6.4);
	TH2F *PT_miss_ene_eMF = new TH2F("PT_miss_ene_eMF", "pT miss OF2 versus COF", 100,0,50000, 100,0,50000);
	TH2F *Phi_miss_comparison_ene = new TH2F("Phi_miss_comparison_ene", "OF2 - phi miss two methods comparison", 64,0,6.4, 64,0,6.4);
	TH2F *Phi_miss_comparison_eMF = new TH2F("Phi_miss_comparison_eMF", "COF - phi miss two methods comparison", 64,0,6.4, 64,0,6.4);
	TH2F *PT_miss_comparison_ene = new TH2F("PT_miss_comparison_ene", "OF2 - pT miss two methods comparison", 100,0,50000, 100,0,50000);
	TH2F *PT_miss_comparison_eMF = new TH2F("PT_miss_comparison_eMF", "COF - pT miss two methods comparison", 100,0,50000, 100,0,50000);

	//Ejet comparison
	TH2F *E_jets_disp_ene_eMF = new TH2F("E_jets_disp_ene_eMF", "Energy of jets ",100,500,2500,100,500,2500);

	//Unused---
	TH2F *E123_R1_ene = new TH2F("E123_R1_ene","Energy in each layer, Region 1, OF2 QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	TH2F *E123_R2_ene = new TH2F("E123_R2_ene","Energy in each layer, Region 2, OF2 QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	TH2F *E123_R3_ene = new TH2F("E123_R3_ene","Energy in each layer, Region 3, OF2 QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	TH2F *E123_R1_eMF = new TH2F("E123_R1_eMF","Energy in each layer, Region 1, COF QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	TH2F *E123_R2_eMF = new TH2F("E123_R2_eMF","Energy in each layer, Region 2, COF QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	TH2F *E123_R3_eMF = new TH2F("E123_R3_eMF","Energy in each layer, Region 3, COF QF < " + SQFcut, 100, 0, 2000, 100, 0, 2000);
	//---
	//E and nr of muons
	TH1F *Ejets_mu_cand1_ene = new TH1F("Ejets_mu_cand1_ene","OF2 - Energy of muons candidates (Zr = 1, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Ejets_mu_cand1_eMF = new TH1F("Ejets_mu_cand1_eMF","COF - Energy of muons candidates (Zr = 1, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Nmu_per_event1_ene = new TH1F("Nmu_per_event1_ene","OF2 - Number of muon candidates per event (Zr = 1, QF<" + SQFcut + ")", 10, 0, 10);
	TH1F *Nmu_per_event1_eMF = new TH1F("Nmu_per_event1_eMF","COF - Number of muon candidates per event (Zr = 1, QF<" + SQFcut + ")", 10, 0, 10);
	TH1F *Ejets_mu_cand2_ene = new TH1F("Ejets_mu_cand2_ene","OF2 - Energy of muons candidates (Zr = 1 + 1 cell, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Ejets_mu_cand2_eMF = new TH1F("Ejets_mu_cand2_eMF","COF - Energy of muons candidates (Zr = 1 + 1 cell, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Nmu_per_event2_ene = new TH1F("Nmu_per_event2_ene","OF2 - Number of muon candidates per event (Zr = 1 + 1 cell, QF<" + SQFcut + ")", 10, 0, 10);
	TH1F *Nmu_per_event2_eMF = new TH1F("Nmu_per_event2_eMF","COF - Number of muon candidates per event (Zr = 1 + 1 cell, QF<" + SQFcut + ")", 10, 0, 10);
	TH1F *Ejets_mu_cand_ene = new TH1F("Ejets_mu_cand_ene","OF2 - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Ejets_mu_cand_eMF = new TH1F("Ejets_mu_cand_eMF","COF - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")", 100, 0, 3000);
	TH1F *Nmu_per_event_ene = new TH1F("Nmu_per_event_ene","OF2 - Number of muon candidates per event (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")", 10, 0, 10);
	TH1F *Nmu_per_event_eMF = new TH1F("Nmu_per_event_eMF","COF - Number of muon candidates per event (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")", 10, 0, 10);
	TH2F *Nmu_per_event_ene_eMF = new TH2F("Nmu_per_event_ene_eMF","Number of muon candidates per event (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")", 10, 0, 10, 10, 0, 10);

	TH2F *Energy_mu_ene_eMF = new TH2F("Energy_mu_ene_eMF","Muon candidates energy OF2/COF QF<" + SQFcut + ")",100,0,100,100,-3,3);

	//Pileup analysis

	TH1F *Ejets_pileup_ene = new TH1F("Ejets_pileup_ene","OF2 - Energy of jets - pileup events (QF<" + SQFcut + ")", 100, 0, 100000);
	TH1F *Ejets_pileup_eMF = new TH1F("Ejets_pileup_eMF","COF - Energy of jets - pileup events (QF<" + SQFcut + ")", 100, 0, 100000);
	TH2F *Ejets_pileup_ene_eMF = new TH2F("Ejets_pileup_ene_eMF","Energy of jets - pileup events (QF<" + SQFcut + ")", 100, 0, 100000, 100, 0, 100000);
	TH1F *Ejets_pileup_z1_ene = new TH1F("Ejets_pileup_z1_ene","OF2 - Energy of jets Zr = 1 - pileup events (Zr = 1, QF<" + SQFcut + ")", 100, 0, 1000);
	TH1F *Ejets_pileup_z1_eMF = new TH1F("Ejets_pileup_z1_eMF","COF - Energy of jets Zr = 1 - pileup events (Zr = 1, QF<" + SQFcut + ")", 100, 0, 1000);
	TH1F *Ejets_pileup_muons_ene = new TH1F("Ejets_pileup_muons_ene","OF2 - Energy of muon candidates - pileup events (Zr = 1, QF<" + SQFcut + ")", 100, 0, 3500);
	TH1F *Ejets_pileup_muons_eMF = new TH1F("Ejets_pileup_muons_eMF","COF - Energy of muon candidates - pileup events (Zr = 1, QF<" + SQFcut + ")", 100, 0, 3500);
	TH2F *Ejets_pileup_muons = new TH2F("Ejets_pileup_muons","Energy of muons candidates - pileup - QF < " + SQFcut, 100, 0, 3000, 100, 0, 3000);
	TH1F *Ejets_pileup_muons_enexeMF = new TH1F("Ejets_pileup_muons_enexeMF","Energy of muon candidates OF2/COF - pileup events (Zr = 1, QF<" + SQFcut + ")", 100, -2, 20);

	TH1F *Nmuon_COFminusOF2_pileup = new TH1F("Nmuon_COFminusOF2_pileup","Number of muon candidates COF minus OF2 (QF<" + SQFcut + ")",20,-10,10);
	TH1F *Njets_COFminusOF2_pileup = new TH1F("Njets_COFminusOF2_pileup","Number of jets COF minus OF2 (QF<" + SQFcut + ")",20,-10,10);
	TH1F *Njets_z1_COFminusOF2_pileup = new TH1F("Njets_z1_COFminusOF2_pileup","Number of jets Zr=1 COF minus OF2 (QF<" + SQFcut + ")",20,-10,10);

	//Ratio of pileup cell, tower and jet

	TH1F *Ratio_cell_pile_energy = new TH1F("Ratio_cell_pile_energy","Ratio of cells energy with pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);
	TH1F *Ratio_tower_pile_energy = new TH1F("Ratio_tower_pile_energy","Ratio of tower energy with pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);
	TH1F *Ratio_jet_pile_energy = new TH1F("Ratio_jet_pile_energy","Ratio of jet energy with pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);

	//Ratio of no pileup cell, tower and jet

	TH1F *Ratio_cell_nopile_energy = new TH1F("Ratio_cell_nopile_energy","Ratio of cells energy without pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);
	TH1F *Ratio_tower_nopile_energy = new TH1F("Ratio_tower_nopile_energy","Ratio of tower energy without pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);
	TH1F *Ratio_jet_nopile_energy = new TH1F("Ratio_jet_nopile_energy","Ratio of jet energy without pileup (COF/OF2) (noisecut_COF="+Snoisecut_eMF+" QF<" + SQFcut + ")",115,-1.5,10);

	//Cell, tower and jet pile up and no pileup energy

	TH2F *Ecell_pileup_ene_eMF = new TH2F("Ecell_pileup_ene_eMF","Energy of pileup cells QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Ecell_nopileup_ene_eMF = new TH2F("Ecell_nopileup_ene_eMF","Energy of no pileup cells QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Etower_pileup_ene_eMF = new TH2F("Etower_pileup_ene_eMF","Energy of pileup towers QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Etower_nopileup_ene_eMF = new TH2F("Etower_nopileup_ene_eMF","Energy of no pileup towers QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Ejet_pileup_ene_eMF = new TH2F("Ejet_pileup_ene_eMF","Energy of pileup jets QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Ejet_nopileup_ene_eMF = new TH2F("Ejet_nopileup_ene_eMF","Energy of no pileup jets QF<" + SQFcut + ")",100,0,5000,100,0,5000);

	TH2F *QFOF2xQFCOF_1 = new TH2F("QFOF2xQFCOF_1","chi2_OF2 x chi2_COF E_OF2>300 E_COF>300",100,0,300,100,0,300);
	TH2F *QFOF2xQFCOF_2 = new TH2F("QFOF2xQFCOF_2","chi2_OF2 x chi2_COF E_OF2>250 E_COF>300",100,0,300,100,0,300);	
	TH2F *QFOF2xQFCOF_3 = new TH2F("QFOF2xQFCOF_3","chi2_OF2 x chi2_COF E_OF2>200 E_COF>300",100,0,300,100,0,300);
	TH2F *QFOF2xQFCOF_4 = new TH2F("QFOF2xQFCOF_4","chi2_OF2 x chi2_COF E_OF2>150 E_COF>300",100,0,300,100,0,300);
	TH2F *QFOF2xQFCOF_5 = new TH2F("QFOF2xQFCOF_5","chi2_OF2 x chi2_COF E_OF2>300 E_COF>0",100,0,300,100,0,300);
	TH2F *QFOF2xQFCOF_6 = new TH2F("QFOF2xQFCOF_6","chi2_OF2 x chi2_COF E_OF2>250 E_COF>0",100,0,300,100,0,300);	
	TH2F *QFOF2xQFCOF_7 = new TH2F("QFOF2xQFCOF_7","chi2_OF2 x chi2_COF E_OF2>200 E_COF>0",100,0,300,100,0,300);
	TH2F *QFOF2xQFCOF_8 = new TH2F("QFOF2xQFCOF_8","chi2_OF2 x chi2_COF E_OF2>150 E_COF>0",100,0,300,100,0,300);

	TH2F *E_OF2_versus_E_COF_chi2 = new TH2F("E_OF2_versus_E_COF_chi2","Energy of channel QF_OF2<" + SQFcut_OF2 + ")",100,0,5000,100,0,5000);
	TH1F *E_OF2_chi2 = new TH1F("E_OF2_chi2","Energy of channel (QF_OF2<" + SQFcut_OF2 + ")",100,0,5000);
	TH1F *E_COF_chi2 = new TH1F("E_COF_chi2","Energy of channel (QF_OF2<" + SQFcut_OF2 + ")",100,0,5000);

	TH2F *Ecell_pileup_OF2alone_ene_eMF = new TH2F("Ecell_pileup_OF2alone_ene_eMF","Energy of pileup cells when just OF2 reconstruct QF<" + SQFcut + ")",100,noisecut_ene,1600,100,-500,noisecut_eMF);
	TH2F *Ecell_pileup_COFalone_eMF_ene = new TH2F("Ecell_pileup_COFalone_eMF_ene","Energy of pileup cells when just COF reconstruct QF<" + SQFcut + ")",100,noisecut_eMF,1600,100,-500,noisecut_ene);

	TH1F *Zr_jetdif_pileup = new TH1F("Zr_jetdif_pileup","Jets: Zr_COF minus Zr_OF2, pileup events (QF<" + SQFcut + ")",20,-10,10);
	TH1F *Energy_alone_jet_pileup_ene = new TH1F("Energy_alone_jet_pileup_ene","Energy of OF2 jets, where COF doesn't reconstruct, pileup events (QF<" + SQFcut + ")",100,300,10000);
	TH1F *Energy_alone_jet_pileup_eMF = new TH1F("Energy_alone_jet_pileup_eMF","Energy of COF jets, where OF2 doesn't reconstruct, pileup events (QF<" + SQFcut + ")",100,300,10000);

	TH1F *Cell_pileup_percent_ene = new TH1F("Cell_pileup_percent_ene", "OF2 - Percent of cells with pileup per event (QF<" + SQFcut + ")",100,0,100);
	TH1F *Cell_pileup_percent_eMF = new TH1F("Cell_pileup_percent_eMF", "COF - Percent of cells with pileup per event (QF<" + SQFcut + ")",100,0,10);
	TH1F *Cell_pileup_eMF_minus_ene = new TH1F("Cell_pileup_eMF_minus_ene", "Number of pileup cells per event COF minus OF2 (QF<" + SQFcut + ")",100,-100,100);

	TH1F *Energy_pile_z1_ene = new TH1F("Energy_pile_z1_ene","OF2 - Energy of jets (Zr=1) with pile up (QF<" + SQFcut + ")",100,0,10000);
	TH1F *Energy_pile_z1_eMF = new TH1F("Energy_pile_z1_eMF","COF - Energy of jets (Zr=1) with pile up (QF<" + SQFcut + ")",100,0,10000);
	TH2F *Energy_pile_z1_ene_eMF = new TH2F("Energy_pile_z1_ene_eMF","Energy of jets (Zr=1) with pile up (noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,10000,100,0,10000);

	TH1F *Energy_jet_pile_ene = new TH1F("Energy_jet_pile_ene","OF2 - Energy of jets (Zr>1) with pile up (QF<" + SQFcut + ")",100,0,50000);
	TH1F *Energy_jet_pile_eMF = new TH1F("Energy_jet_pile_eMF","COF - Energy of jets (Zr>1) with pile up (QF<" + SQFcut + ")",100,0,50000);
	TH2F *Energy_jet_pile_ene_eMF = new TH2F("Energy_jet_pile_ene_eMF","Energy of jets (Zr>1) with pile up (noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,100000,100,0,100000);

	TH1F *Energy_mucand_pile_ene = new TH1F("Energy_mucand_pile_ene","OF2 - Energy of muon candidates with pile up (QF<" + SQFcut + ")",100,0,10000);
	TH1F *Energy_mucand_pile_eMF = new TH1F("Energy_mucand_pile_eMF","COF - Energy of muon candidates with pile up (QF<" + SQFcut + ")",100,0,10000);
	TH2F *Energy_mucand_pile_ene_eMF = new TH2F("Energy_mucand_pile_ene_eMF","Energy of muon candidates with pile up (noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,5000,100,0,5000);

	TH1F *Number_pile_z1_ene = new TH1F("Number_pile_z1_ene","OF2 - Number of pileup cells per jets (Zr=1) (QF<" + SQFcut + ")",16,0,16);
	TH1F *Number_pile_z1_eMF = new TH1F("Number_pile_z1_eMF","COF - Number of pileup cells per jets (Zr=1) (QF<" + SQFcut + ")",16,0,16);
	TH2F *Number_pile_z1_ene_eMF = new TH2F("Number_pile_z1_ene_eMF","Number of pileup cells per jets (Zr=1) (noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",16,0,16,16,0,16);

	TH1F *Number_jet_pile_ene = new TH1F("Number_jet_pile_ene","OF2 - Number of pileup cells per jets (Zr>1) (QF<" + SQFcut + ")",50,0,50);
	TH1F *Number_jet_pile_eMF = new TH1F("Number_jet_pile_eMF","COF - Number of pileup cells per jets (Zr>1) (QF<" + SQFcut + ")",50,0,50);
	TH2F *Number_jet_pile_ene_eMF = new TH2F("Number_jet_pile_ene_eMF","Number of pileup cells per jets (Zr>1) (noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",50,0,50,50,0,50);

	TH1F *Number_mucand_pile_ene = new TH1F("Number_mucand_pile_ene","OF2 - Number of pileup cells per jets (muon candidates QF<" + SQFcut + ")",16,0,16);
	TH1F *Number_mucand_pile_eMF = new TH1F("Number_mucand_pile_eMF","COF - Number of pileup cells per jets (muon candidates QF<" + SQFcut + ")",16,0,16);
	TH2F *Number_mucand_pile_ene_eMF = new TH2F("Number_mucand_pile_ene_eMF","Number of pileup cells per jets (muon candidates noisecutCOF="+Snoisecut_eMF+" QF<" + SQFcut + ")",16,0,16,16,0,16);

	TH1F *Energy_alone_pile_z1_ene = new TH1F("Energy_alone_pile_z1_ene","OF2 - Energy of alone jets (Zr=1) with pile up (QF<" + SQFcut + ")",100,0,5000);
	TH1F *Energy_alone_pile_z1_eMF = new TH1F("Energy_alone_pile_z1_eMF","COF - Energy of alone jets (Zr=1) with pile up (QF<" + SQFcut + ")",100,0,5000);
	TH1F *Energy_alone_jet_pile_ene = new TH1F("Energy_alone_jet_pile_ene","OF2 - Energy of alone jets (Zr>1) with pile up (QF<" + SQFcut + ")",100,0,50000);
	TH1F *Energy_alone_jet_pile_eMF = new TH1F("Energy_alone_jet_pile_eMF","COF - Energy of alone jets (Zr>1) with pile up (QF<" + SQFcut + ")",100,0,50000);
	TH1F *Energy_alone_mucand_pile_ene = new TH1F("Energy_alone_mucand_pile_ene","OF2 - Energy of alone muon candidates with pile up (QF<" + SQFcut + ")",100,0,2000);
	TH1F *Energy_alone_mucand_pile_eMF = new TH1F("Energy_alone_mucand_pile_eMF","COF - Energy of alone muon candidates with pile up (QF<" + SQFcut + ")",100,0,2000);
	TH1F *Ejets_alone_pileup_ene = new TH1F("Ejets_alone_pileup_ene","OF2 - Energy of alone jets - pileup events (QF<" + SQFcut + ")", 100, 0, 100000);
	TH1F *Ejets_alone_pileup_eMF = new TH1F("Ejets_alone_pileup_eMF","COF - Energy of alone jets - pileup events (QF<" + SQFcut + ")", 100, 0, 100000);

	TH1F *Number_alone_pile_z1_ene = new TH1F("Number_alone_pile_z1_ene","OF2 - Number of pileup cells per alone jets (Zr=1) (QF<" + SQFcut + ")",20,0,20);
	TH1F *Number_alone_pile_z1_eMF = new TH1F("Number_alone_pile_z1_eMF","COF - Number of pileup cells per alone jets (Zr=1) (QF<" + SQFcut + ")",20,0,20);

	TH1F *Number_alone_jet_pile_ene = new TH1F("Number_alone_jet_pile_ene","OF2 - Number of pileup cells per alone jets (Zr>1) (QF<" + SQFcut + ")",20,0,20);
	TH1F *Number_alone_jet_pile_eMF = new TH1F("Number_alone_jet_pile_eMF","COF - Number of pileup cells per alone jets (Zr>1) (QF<" + SQFcut + ")",20,0,20);

	TH1F *Number_alone_mucand_pile_ene = new TH1F("Number_alone_mucand_pile_ene","OF2 - Number of pileup cells per alone jets (muon candidates QF<" + SQFcut + ")",10,0,10);
	TH1F *Number_alone_mucand_pile_eMF = new TH1F("Number_alone_mucand_pile_eMF","COF - Number of pileup cells per alone jets (muon candidates QF<" + SQFcut + ")",10,0,10);

	TH2F *Energy_alone_cell_pile_ene_eMF = new TH2F("Energy_alone_cell_pile_ene_eMF","Energy of cells where OF2 reconstruct and COF don't versus COF with E_cell>0 (QF<" + SQFcut + ")",100,0,6000,100,0,300);
	TH2F *Energy_alone_muon_pile_ene_eMF = new TH2F("Energy_alone_muon_pile_ene_eMF","Energy of muon where OF2 reconstruct and COF don't versus COF with E_cell>0 (QF<" + SQFcut + ")",100,0,6000,100,0,300);
	TH2F *Energy_alone_tower_pile_ene_eMF = new TH2F("Energy_tower_cell_pile_ene_eMF","Energy of tower where OF2 reconstruct a cell and COF don't versus COF with E_cell>0 (QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Energy_alone_tower2_pile_ene_eMF = new TH2F("Energy_tower2_cell_pile_ene_eMF","Tower Energy (at least one cell with Energy < 300 MeV ????? noisecut_eMF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Energy_alone_tower2_pileCOF_ene_eMF = new TH2F("Energy_tower2_cell_pileCOF_ene_eMF","Tower Energy (at least one cell with E_OF2 < 300 MeV noisecut_eMF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,5000,100,0,5000);
	TH2F *Energy_alone_tower2_pileOF2_ene_eMF = new TH2F("Energy_tower2_cell_pileOF2_ene_eMF","Tower Energy (at least one cell with E_COF < 300 MeV noisecut_eMF="+Snoisecut_eMF+" QF<" + SQFcut + ")",100,0,5000,100,0,5000);

	TH2F *Energy_alone_muon_pile_ene = new TH2F("Energy_alone_muon_pile_ene","OF2 - Energy of muon where OF2 reconstruct and COF don't versus COF with E_cell>0 (QF<" + SQFcut + ")",100,0,6000,100,0,300);
	TH2F *Energy_alone_muon_pile_eMF = new TH2F("Energy_alone_muon_pile_eMF","COF - Energy of muon with pileup and noisecut_eMF>"+Snoisecut_eMF+" E_cell>0 (QF<" + SQFcut + ")",100,0,6000,100,0,300);

	TH1F *White = new TH1F("White", "White",100,-100,100);
	White->SetLineColor(0);

	//Setting the line color
	Njets_per_event_ene->SetLineColor(2);
	Njets_per_event_z1_ene->SetLineColor(2);
	Njets_per_event_zequal1_ene->SetLineColor(2);
	Njets_per_event_z2_ene->SetLineColor(2);
	Ejets_ene->SetLineColor(2);
	Ejets_z1_ene->SetLineColor(2);
	Ejets_zequal1_ene->SetLineColor(2);
	Ejets_z2_ene->SetLineColor(2);
	phi_cm_ene->SetLineColor(2);
	eta_cm_ene->SetLineColor(2);
	Ratio_E_l_ene->SetLineColor(2);
	Ejets_mu_cand1_ene->SetLineColor(2);
	Nmu_per_event1_ene->SetLineColor(2);
	Ejets_mu_cand2_ene->SetLineColor(2);
	Nmu_per_event2_ene->SetLineColor(2);
	Ejets_mu_cand_ene->SetLineColor(2);
	Nmu_per_event_ene->SetLineColor(2);
	Ejets_pileup_ene->SetLineColor(2);
	Ejets_pileup_muons_ene->SetLineColor(2);
	Ejets_pileup_z1_ene->SetLineColor(2);
	Energy_alone_jet_pileup_ene->SetLineColor(2);
	Cell_pileup_percent_ene->SetLineColor(2);
	Energy_pile_z1_ene->SetLineColor(2);
	Energy_jet_pile_ene->SetLineColor(2);
	Energy_mucand_pile_ene->SetLineColor(2);
	Energy_alone_pile_z1_ene->SetLineColor(2);
	Energy_alone_jet_pile_ene->SetLineColor(2);
	Energy_alone_mucand_pile_ene->SetLineColor(2);	
	Number_pile_z1_ene->SetLineColor(2);
	Number_jet_pile_ene->SetLineColor(2);
	Number_mucand_pile_ene->SetLineColor(2);

	Njets_per_event_eMF->SetLineColor(4);
	Njets_per_event_z1_eMF->SetLineColor(4);
	Njets_per_event_zequal1_eMF->SetLineColor(4);
	Njets_per_event_z2_eMF->SetLineColor(4);
	Ejets_eMF->SetLineColor(4);
	Ejets_z1_eMF->SetLineColor(4);
	Ejets_zequal1_eMF->SetLineColor(4);
	Ejets_z2_eMF->SetLineColor(4);
	phi_cm_eMF->SetLineColor(4);
	eta_cm_eMF->SetLineColor(4);
	PT_miss_eMF->SetLineColor(4);
	PT_miss_pileev_ene2->SetLineColor(2);
	PT_miss_pileev_eMF2->SetLineColor(4);
	Phi_miss_eMF->SetLineColor(4);
	eta_miss_eMF->SetLineColor(4);
	PT_miss_eMF2->SetLineColor(4);
	Phi_miss_eMF2->SetLineColor(4);
	eta_miss_eMF2->SetLineColor(4);
	Ratio_E_l_eMF->SetLineColor(4);
	Ejets_mu_cand1_eMF->SetLineColor(4);
	Nmu_per_event1_eMF->SetLineColor(4);
	Ejets_mu_cand2_eMF->SetLineColor(4);
	Nmu_per_event2_eMF->SetLineColor(4);
	Ejets_mu_cand_eMF->SetLineColor(4);
	Nmu_per_event_eMF->SetLineColor(4);
	Ejets_pileup_eMF->SetLineColor(4);
	Ejets_pileup_muons_eMF->SetLineColor(4);
	Ejets_pileup_z1_eMF->SetLineColor(4);
	Energy_alone_jet_pileup_eMF->SetLineColor(4);
	Cell_pileup_percent_eMF->SetLineColor(4);
	Energy_pile_z1_eMF->SetLineColor(4);
	Energy_jet_pile_eMF->SetLineColor(4);
	Energy_mucand_pile_eMF->SetLineColor(4);
	Energy_alone_pile_z1_eMF->SetLineColor(4);
	Energy_alone_jet_pile_eMF->SetLineColor(4);
	Energy_alone_mucand_pile_eMF->SetLineColor(4);
	Number_pile_z1_eMF->SetLineColor(4);
	Number_jet_pile_eMF->SetLineColor(4);
	Number_mucand_pile_eMF->SetLineColor(4);

	Ejets_all3layers_ene->SetLineColor(2);
	Ejets_all3layers_eMF->SetLineColor(4);
	Ejets_layers_ene->SetLineColor(2);
	Ejets_layers_eMF->SetLineColor(4);
	PT_miss_ene->SetLineColor(2);
	Phi_miss_ene->SetLineColor(2);
	eta_miss_ene->SetLineColor(2);
	PT_miss_ene2->SetLineColor(2);
	Phi_miss_ene2->SetLineColor(2);
	eta_miss_ene2->SetLineColor(2);

	Zr_dist_ene->SetLineColor(2);
	Zr_dist_eMF->SetLineColor(4);

	DeltaR_dist_ene->SetLineColor(2);
	DeltaR_dist_eMF->SetLineColor(4);

	Elayers12_ene->SetLineColor(2);
	Elayers13_ene->SetLineColor(2);
	Elayers23_ene->SetLineColor(2);
	Elayersall_ene->SetLineColor(2);
	Elayers12_eMF->SetLineColor(4);
	Elayers13_eMF->SetLineColor(4);
	Elayers23_eMF->SetLineColor(4);
	Elayersall_eMF->SetLineColor(4);

	Emu_layer1_eMF->SetLineColor(4);
	Emu_layer2_eMF->SetLineColor(4);
	Emu_layer3_eMF->SetLineColor(4);
	Emu_layer1_ene->SetLineColor(2);
	Emu_layer2_ene->SetLineColor(2);
	Emu_layer3_ene->SetLineColor(2);

	Elayers1_ene->SetLineColor(2);
	Elayers2_ene->SetLineColor(2);
	Elayers3_ene->SetLineColor(2);
	Elayers1_eMF->SetLineColor(4);
	Elayers2_eMF->SetLineColor(4);
	Elayers3_eMF->SetLineColor(4);


	E_OF2_chi2->SetLineColor(2);
	E_COF_chi2->SetLineColor(4);
	E_OF2_versus_E_COF_chi2->GetXaxis()->SetTitle("OF2");
	E_OF2_versus_E_COF_chi2->GetYaxis()->SetTitle("COF");

	//AXION PARTICLE SEARCH------------------
	E_nomuon_l1_ene->SetLineColor(2);
	E_nomuon_l2_ene->SetLineColor(2);
	E_nomuon_l3_ene->SetLineColor(2);
	Emm_nomuon_l1_ene->SetLineColor(2);
	Emm_nomuon_l2_ene->SetLineColor(2);
	Emm_nomuon_l3_ene->SetLineColor(2);
	E_nomuon_l1_eMF->SetLineColor(4);
	E_nomuon_l2_eMF->SetLineColor(4);
	E_nomuon_l3_eMF->SetLineColor(4);
	Emm_nomuon_l1_eMF->SetLineColor(4);
	Emm_nomuon_l2_eMF->SetLineColor(4);
	Emm_nomuon_l3_eMF->SetLineColor(4);
	E_nomuon_l1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	E_nomuon_l2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	E_nomuon_l3_ene->GetXaxis()->SetTitle("Energy (MeV)");
	E_nomuon_l1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	E_nomuon_l2_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	E_nomuon_l3_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Emm_nomuon_l1_ene->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	Emm_nomuon_l2_ene->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	Emm_nomuon_l3_ene->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	Emm_nomuon_l1_eMF->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	Emm_nomuon_l2_eMF->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	Emm_nomuon_l3_eMF->GetXaxis()->SetTitle("Energy/cm (MeV/cm)");
	QFOF2xQFCOF_1->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_1->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_2->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_2->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_3->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_3->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_4->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_4->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_5->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_5->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_6->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_6->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_7->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_7->GetYaxis()->SetTitle("chi2_COF");
	QFOF2xQFCOF_8->GetXaxis()->SetTitle("chi2_OF2");
	QFOF2xQFCOF_8->GetYaxis()->SetTitle("chi2_COF");

	Emm_nomuon_l12_ene->GetXaxis()->SetTitle("MeV/cm layer 1");
	Emm_nomuon_l12_ene->GetYaxis()->SetTitle("MeVcm layer 2");
	Emm_nomuon_l13_ene->GetXaxis()->SetTitle("MeV/cm layer 1");
	Emm_nomuon_l13_ene->GetYaxis()->SetTitle("MeVcm layer 3");
	Emm_nomuon_l23_ene->GetXaxis()->SetTitle("MeV/cm layer 2");
	Emm_nomuon_l23_ene->GetYaxis()->SetTitle("MeV/cm layer 3");
	Emm_nomuon_l12_eMF->GetXaxis()->SetTitle("MeV/cm layer 1");
	Emm_nomuon_l12_eMF->GetYaxis()->SetTitle("MeV/cm layer 2");
	Emm_nomuon_l13_eMF->GetXaxis()->SetTitle("MeV/cm layer 1");
	Emm_nomuon_l13_eMF->GetYaxis()->SetTitle("MeV/cm layer 3");
	Emm_nomuon_l23_eMF->GetXaxis()->SetTitle("MeV/cm layer 2");
	Emm_nomuon_l23_eMF->GetYaxis()->SetTitle("MeV/cm layer 3");

	Emm_all_nocut_ene->SetLineColor(2);
	Emm_all_nocut_eMF->SetLineColor(4);
	Emm_all_nocut_ene->GetXaxis()->SetTitle("MeV/cm");
	Emm_all_nocut_eMF->GetXaxis()->SetTitle("MeV/cm");

	Ejets_1layers_ene->SetLineColor(2);
	Ejets_1layers_eMF->SetLineColor(4);
	Ejets_2layers_ene->SetLineColor(2);
	Ejets_2layers_eMF->SetLineColor(4);
	Ejets_12layers_ene->SetLineColor(2);
	Ejets_12layers_eMF->SetLineColor(4);
	Ejets_23layers_ene->SetLineColor(2);
	Ejets_23layers_eMF->SetLineColor(4);
	Ejets_3layers_ene->SetLineColor(2);
	Ejets_3layers_eMF->SetLineColor(4);
	Ejets_1layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_1layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_2layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_2layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_12layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_12layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_23layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_23layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_3layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_3layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	N_channel_per_jet_ene->GetXaxis()->SetTitle("Number of channels");
	N_channel_per_jet_eMF->GetXaxis()->SetTitle("Number of channels");
	N_channel_per_jet_ene->SetLineColor(2);
	N_channel_per_jet_eMF->SetLineColor(4);
	N_channel_per_jet_Energy_ene->GetXaxis()->SetTitle("Number of channels");
	N_channel_per_jet_Energy_ene->GetYaxis()->SetTitle("Energy (MeV)");
	N_channel_per_jet_Energy_eMF->GetXaxis()->SetTitle("Number of channels");
	N_channel_per_jet_Energy_eMF->GetYaxis()->SetTitle("Energy (MeV)");

	THETA_Particle_ene->SetLineColor(2);
	THETA_Particle_eMF->SetLineColor(4);
	THETA_Particle_ene->GetXaxis()->SetTitle("theta (rad)");
	THETA_Particle_eMF->GetXaxis()->SetTitle("theta (rad)");
	MASS_Particle_ene->SetLineColor(2);
	MASS_Particle_eMF->SetLineColor(4);
	MASS_Particle_ene->GetXaxis()->SetTitle("MeV");
	MASS_Particle_eMF->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l1_ene->SetLineColor(2);
	MASS_Particle_l1_eMF->SetLineColor(4);
	MASS_Particle_l1_ene->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l1_eMF->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l2_ene->SetLineColor(2);
	MASS_Particle_l2_eMF->SetLineColor(4);
	MASS_Particle_l2_ene->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l2_eMF->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l3_ene->SetLineColor(2);
	MASS_Particle_l3_eMF->SetLineColor(4);
	MASS_Particle_l3_ene->GetXaxis()->SetTitle("MeV");
	MASS_Particle_l3_eMF->GetXaxis()->SetTitle("MeV");
	//---------------------------------------

	//Setting the X and Y axis

	Cell_pileup_percent_ene->GetXaxis()->SetTitle("%");
	Cell_pileup_percent_eMF->GetXaxis()->SetTitle("%");

	Njets_per_event_ene->GetXaxis()->SetTitle("Number of jets");
	Ejets_pileup_muons->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ejets_pileup_muons->GetYaxis()->SetTitle("COF Energy (MeV)");
	Njets_per_event_zequal1_ene->GetXaxis()->SetTitle("Number of jets");
	Njets_per_event_z2_ene->GetXaxis()->SetTitle("Number of jets");
	Njets_per_event_eMF->GetXaxis()->SetTitle("Number of jets");
	Njets_per_event_z1_eMF->GetXaxis()->SetTitle("Number of jets");
	Njets_per_event_zequal1_eMF->GetXaxis()->SetTitle("Number of jets");
	Njets_per_event_z2_eMF->GetXaxis()->SetTitle("Number of jets");

	Nmu_per_event_ene_eMF->GetXaxis()->SetTitle("OF2");
	Nmu_per_event_ene_eMF->GetYaxis()->SetTitle("COF");

	Ejets_mu_cand_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_pileup_z1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_pileup_z1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_mu_cand_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Nmu_per_event_ene->GetXaxis()->SetTitle("Number of muons candidates");
	Nmu_per_event_eMF->GetXaxis()->SetTitle("Number of muons candidates");
	Ejets_mu_cand1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_mu_cand1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Nmu_per_event1_ene->GetXaxis()->SetTitle("Number of muons candidates");
	Nmu_per_event1_eMF->GetXaxis()->SetTitle("Number of muons candidates");
	Ejets_mu_cand2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_mu_cand2_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Nmu_per_event2_ene->GetXaxis()->SetTitle("Number of muons candidates");
	Nmu_per_event2_eMF->GetXaxis()->SetTitle("Number of muons candidates");
	Energy_pile_z1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_jet_pile_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_mucand_pile_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_pile_z1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_jet_pile_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_mucand_pile_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Number_pile_z1_ene->GetXaxis()->SetTitle("Number of pile up cells");
	Number_jet_pile_ene->GetXaxis()->SetTitle("Number of pile up cells");
	Number_mucand_pile_ene->GetXaxis()->SetTitle("Number of pile up cells");
	Number_pile_z1_eMF->GetXaxis()->SetTitle("Number of pile up cells");
	Number_jet_pile_eMF->GetXaxis()->SetTitle("Number of pile up cells");
	Number_mucand_pile_eMF->GetXaxis()->SetTitle("Number of pile up cells");
	Energy_alone_pile_z1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_alone_jet_pile_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_alone_mucand_pile_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_alone_pile_z1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_alone_jet_pile_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Energy_alone_mucand_pile_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_pileup_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_pileup_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_pileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ejets_pileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");

	Energy_pile_z1_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Energy_jet_pile_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Energy_mucand_pile_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Energy_pile_z1_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Energy_jet_pile_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Energy_mucand_pile_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Number_pile_z1_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Number_jet_pile_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Number_mucand_pile_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Number_pile_z1_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Number_jet_pile_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Number_mucand_pile_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Energy_alone_tower2_pile_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Energy_alone_tower2_pile_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Energy_alone_tower2_pileCOF_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Energy_alone_tower2_pileOF2_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Ecell_pileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ecell_pileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Ecell_nopileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ecell_nopileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Etower_pileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Etower_pileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Etower_nopileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Etower_nopileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Ejet_pileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ejet_pileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Ejet_nopileup_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ejet_nopileup_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");

	Ejets_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_z1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_zequal1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_z2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_z1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_zequal1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_z2_eMF->GetXaxis()->SetTitle("Energy (MeV)");

	Ejets_all3layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_all3layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_layers_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Ejets_layers_eMF->GetXaxis()->SetTitle("Energy (MeV)");

	Emu_layer1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Emu_layer2_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Emu_layer3_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Emu_layer1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Emu_layer2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Emu_layer3_ene->GetXaxis()->SetTitle("Energy (MeV)");

	Ecell_pileup_OF2alone_ene_eMF->GetXaxis()->SetTitle("OF2 Energy (MeV)");
	Ecell_pileup_OF2alone_ene_eMF->GetYaxis()->SetTitle("COF Energy (MeV)");
	Ecell_pileup_COFalone_eMF_ene->GetXaxis()->SetTitle("COF Energy (MeV)");
	Ecell_pileup_COFalone_eMF_ene->GetYaxis()->SetTitle("OF2 Energy (MeV)");

	Phi_miss_ene->GetXaxis()->SetTitle("Phi");
	Phi_miss_eMF->GetXaxis()->SetTitle("Phi");
	PT_miss_eMF->GetXaxis()->SetTitle("pTmiss (MeV)");
	PT_miss_ene->GetXaxis()->SetTitle("pTmiss (MeV)");
	PT_miss_pileev_ene2->GetXaxis()->SetTitle("pTmiss (MeV)");
	PT_miss_pileev_eMF2->GetXaxis()->SetTitle("pTmiss (MeV)");

	Phi_miss_ene2->GetXaxis()->SetTitle("Phi");
	Phi_miss_eMF2->GetXaxis()->SetTitle("Phi");
	PT_miss_eMF2->GetXaxis()->SetTitle("pTmiss (MeV)");
	PT_miss_ene2->GetXaxis()->SetTitle("pTmiss (MeV)");
	Njets_per_event_ene_eMF->GetXaxis()->SetTitle("OF2");
	Njets_per_event_ene_eMF->GetYaxis()->SetTitle("COF");
	Njets_per_event_z1_ene_eMF->GetXaxis()->SetTitle("OF2");
	Njets_per_event_z1_ene_eMF->GetYaxis()->SetTitle("COF");
	Njets_per_event_zequal1_ene_eMF->GetXaxis()->SetTitle("OF2");
	Njets_per_event_zequal1_ene_eMF->GetYaxis()->SetTitle("COF");
	Njets_per_event_z2_ene_eMF->GetXaxis()->SetTitle("OF2");
	Njets_per_event_z2_ene_eMF->GetYaxis()->SetTitle("COF");
	Nmu_per_event_z1_ene_eMF1->GetXaxis()->SetTitle("OF2");
	Nmu_per_event_z1_ene_eMF1->GetYaxis()->SetTitle("COF");

	phi_cm_ene->GetXaxis()->SetTitle("Phi");
	phi_cm_eMF->GetXaxis()->SetTitle("Phi");
	eta_cm_ene->GetXaxis()->SetTitle("Eta");
	eta_cm_eMF->GetXaxis()->SetTitle("Eta");
	
	EMFxQF->GetXaxis()->SetTitle("Energy (MeV)");
	EnexQF->GetXaxis()->SetTitle("Energy (MeV)");
	EMFxQF2->GetXaxis()->SetTitle("Energy (MeV)");
	EnexQF2->GetXaxis()->SetTitle("Energy (MeV)");
	EMFxQF2->GetYaxis()->SetTitle("Quality Factor");
	EnexQF2->GetYaxis()->SetTitle("Quality Factor");

	EneXQF_100->GetXaxis()->SetTitle("Energy (MeV)");
	EMFXQF_100->GetXaxis()->SetTitle("Energy (MeV)");
	EneXQF_50->GetXaxis()->SetTitle("Energy (MeV)");
	EMFXQF_50->GetXaxis()->SetTitle("Energy (MeV)");
	EneXQF_30->GetXaxis()->SetTitle("Energy (MeV)");
	EMFXQF_30->GetXaxis()->SetTitle("Energy (MeV)");

	Zr_dist_ene->GetXaxis()->SetTitle("Zr");
	Zr_dist_eMF->GetXaxis()->SetTitle("Zr");

	DeltaR_dist_ene->GetXaxis()->SetTitle("Delta R");
	DeltaR_dist_eMF->GetXaxis()->SetTitle("Delta R");

	Elayers12_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers12_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers23_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers23_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers13_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers13_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers1_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers1_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers2_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers3_eMF->GetXaxis()->SetTitle("Energy (MeV)");
	Elayers3_ene->GetXaxis()->SetTitle("Energy (MeV)");

	Elayersall_ene->GetXaxis()->SetTitle("Energy (MeV)");
	Elayersall_eMF->GetXaxis()->SetTitle("Energy (MeV)");

	Elayers12_ene->GetYaxis()->SetTitle("Zr");
	Elayers23_ene->GetYaxis()->SetTitle("Zr");
	Elayers13_ene->GetYaxis()->SetTitle("Zr");
	Elayers1_ene->GetYaxis()->SetTitle("Zr");
	Elayers2_ene->GetYaxis()->SetTitle("Zr");
	Elayers3_ene->GetYaxis()->SetTitle("Zr");
	Elayers12_eMF->GetYaxis()->SetTitle("Zr");
	Elayers23_eMF->GetYaxis()->SetTitle("Zr");
	Elayers13_eMF->GetYaxis()->SetTitle("Zr");
	Elayers1_eMF->GetYaxis()->SetTitle("Zr");
	Elayers2_eMF->GetYaxis()->SetTitle("Zr");
	Elayers3_eMF->GetYaxis()->SetTitle("Zr");

	Elayers12_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers23_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers13_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers1_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers2_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers3_ene->GetYaxis()->SetTitleOffset(2.0);
	Elayers12_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers23_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers13_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers1_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers2_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers3_eMF->GetYaxis()->SetTitleOffset(2.0);
	Elayers12_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers23_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers13_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers1_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers2_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers3_ene->GetXaxis()->SetTitleOffset(2.0);
	Elayers12_eMF->GetXaxis()->SetTitleOffset(2.0);
	Elayers23_eMF->GetXaxis()->SetTitleOffset(2.0);
	Elayers13_eMF->GetXaxis()->SetTitleOffset(2.0);
	Elayers1_eMF->GetXaxis()->SetTitleOffset(2.0);
	Elayers2_eMF->GetXaxis()->SetTitleOffset(2.0);
	Elayers3_eMF->GetXaxis()->SetTitleOffset(2.0);

	Elayers12_ene->GetYaxis()->CenterTitle();
	Elayers23_ene->GetYaxis()->CenterTitle();
	Elayers13_ene->GetYaxis()->CenterTitle();
	Elayers1_ene->GetYaxis()->CenterTitle();
	Elayers2_ene->GetYaxis()->CenterTitle();
	Elayers3_ene->GetYaxis()->CenterTitle();
	Elayers12_eMF->GetYaxis()->CenterTitle();
	Elayers23_eMF->GetYaxis()->CenterTitle();
	Elayers13_eMF->GetYaxis()->CenterTitle();
	Elayers1_eMF->GetYaxis()->CenterTitle();
	Elayers2_eMF->GetYaxis()->CenterTitle();
	Elayers3_eMF->GetYaxis()->CenterTitle();
	Elayers12_ene->GetXaxis()->CenterTitle();
	Elayers23_ene->GetXaxis()->CenterTitle();
	Elayers13_ene->GetXaxis()->CenterTitle();
	Elayers1_ene->GetXaxis()->CenterTitle();
	Elayers2_ene->GetXaxis()->CenterTitle();
	Elayers3_ene->GetXaxis()->CenterTitle();
	Elayers12_eMF->GetXaxis()->CenterTitle();
	Elayers23_eMF->GetXaxis()->CenterTitle();
	Elayers13_eMF->GetXaxis()->CenterTitle();
	Elayers1_eMF->GetXaxis()->CenterTitle();
	Elayers2_eMF->GetXaxis()->CenterTitle();
	Elayers3_eMF->GetXaxis()->CenterTitle();

	Ratio_E_l_ene->GetXaxis()->SetTitle("Energy/Length (MeV/cm)");
	Ratio_E_l_eMF->GetXaxis()->SetTitle("Energy/Length (MeV/cm)");

	E123_R1_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R1_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R1_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R2_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R2_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R2_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R3_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R3_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R3_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R1_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R1_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R1_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R2_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R2_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R2_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R3_ene->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R3_ene->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R3_ene->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R1_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R1_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R1_eMF->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R2_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R2_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R2_eMF->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R3_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R3_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R3_eMF->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R1_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R1_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R1_eMF->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R2_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R2_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R2_eMF->GetZaxis()->SetTitle("E3 (MeV)");
	E123_R3_eMF->GetXaxis()->SetTitle("E1 (MeV)");
	E123_R3_eMF->GetYaxis()->SetTitle("E2 (MeV)");
	E123_R3_eMF->GetZaxis()->SetTitle("E3 (MeV)");


	Phi_miss_ene_eMF->GetXaxis()->SetTitle("phi miss OF2");
	Phi_miss_ene_eMF->GetYaxis()->SetTitle("phi miss COF");
	PT_miss_ene_eMF->GetXaxis()->SetTitle("pT miss OF2 (MeV)");
	PT_miss_ene_eMF->GetYaxis()->SetTitle("pT miss COF (MeV)");

	Phi_miss_comparison_ene->GetXaxis()->SetTitle("cell by cell");
	Phi_miss_comparison_ene->GetYaxis()->SetTitle("center of energy method");
	Phi_miss_comparison_eMF->GetXaxis()->SetTitle("cell by cell");
	Phi_miss_comparison_eMF->GetYaxis()->SetTitle("center of energy method");
	PT_miss_comparison_ene->GetXaxis()->SetTitle("cell by cell");
	PT_miss_comparison_ene->GetYaxis()->SetTitle("center of energy method");
	PT_miss_comparison_eMF->GetXaxis()->SetTitle("cell by cell");
	PT_miss_comparison_eMF->GetYaxis()->SetTitle("center of energy method");

	TH1F *E_isolated = new TH1F("E_isolated","Isolated Energy", 100, 0, 3000);

	if (noisecut_eMF == 100) E_isolated->SetTitle("Isolated Energy with Noise Cut = 100");
	if (noisecut_eMF == 200) E_isolated->SetTitle("Isolated Energy with Noise Cut = 200");
	if (noisecut_eMF == 300) E_isolated->SetTitle("Isolated Energy with Noise Cut = 300");

	//End made by Joao Pedro!

	cout<<"histo files created"<< endl;

	cout<<"output txt files created"<< endl;

	float const PI=3.141596;

	// need to fill the eta (rapidity) information for each channel
	// according to the figure each channel has a width of 0.1 in eta therefore I am seting the middle value for it
	// for example, if LBA D0(channel=0) eta varies between 0.0 and 0.1 therefore I set ETA=0.05

	//partitions LBA and LBC

	for(int ipart=0; ipart<2; ipart++){

		partner[ipart][0]=0;	
		partner[ipart][1]=4;
		partner[ipart][2]=3;
		partner[ipart][3]=2;
		partner[ipart][4]=1;
		partner[ipart][5]=8;
		partner[ipart][6]=7;
		partner[ipart][7]=6;
		partner[ipart][8]=5;
		partner[ipart][9]=10;
		partner[ipart][10]=9;
		partner[ipart][11]=12;
		partner[ipart][12]=11;
		partner[ipart][13]=14;
		partner[ipart][14]=13;
		partner[ipart][15]=18;
		partner[ipart][16]=17;
		partner[ipart][17]=16;
		partner[ipart][18]=15;
		partner[ipart][19]=20;
		partner[ipart][20]=19;
		partner[ipart][21]=22;
		partner[ipart][22]=21;
		partner[ipart][23]=26;
		partner[ipart][24]=25;
		partner[ipart][25]=24;
		partner[ipart][26]=23;
		partner[ipart][27]=28;
		partner[ipart][28]=27;
		partner[ipart][29]=32;
		partner[ipart][30]=0;
		partner[ipart][31]=0;
		partner[ipart][32]=29;
		partner[ipart][33]=34;
		partner[ipart][34]=33;
		partner[ipart][35]=38;
		partner[ipart][36]=37;
		partner[ipart][37]=36;
		partner[ipart][38]=35;
		partner[ipart][39]=40;
		partner[ipart][40]=39;
		partner[ipart][41]=44;
		partner[ipart][42]=47;
		partner[ipart][43]=0;
		partner[ipart][44]=41;
		partner[ipart][45]=46;
		partner[ipart][46]=45;
		partner[ipart][47]=42;

    	for(int imodule=0; imodule<64; imodule++){
   		 	ETA[ipart][imodule][0]=0.05;
    		ETA[ipart][imodule][1]=0.05;
    		ETA[ipart][imodule][2]=0.05;
    		ETA[ipart][imodule][3]=0.05;
    		ETA[ipart][imodule][4]=0.05;
    		ETA[ipart][imodule][5]=0.15;
    		ETA[ipart][imodule][6]=0.15;
    		ETA[ipart][imodule][7]=0.15;
    		ETA[ipart][imodule][8]=0.15;
   		 	ETA[ipart][imodule][9]=0.25;
    		ETA[ipart][imodule][10]=0.25;
    		ETA[ipart][imodule][11]=0.25;
    		ETA[ipart][imodule][12]=0.25;
   			ETA[ipart][imodule][13]=0.20;
    		ETA[ipart][imodule][14]=0.20;
    		ETA[ipart][imodule][15]=0.35;
    		ETA[ipart][imodule][16]=0.35;
    		ETA[ipart][imodule][17]=0.35;
    		ETA[ipart][imodule][18]=0.35;
    		ETA[ipart][imodule][19]=0.45;
    		ETA[ipart][imodule][20]=0.45;
    		ETA[ipart][imodule][21]=0.45;
    		ETA[ipart][imodule][22]=0.45;
    		ETA[ipart][imodule][23]=0.55;
    		ETA[ipart][imodule][24]=0.40;
    		ETA[ipart][imodule][25]=0.40;
    		ETA[ipart][imodule][26]=0.55;
    		ETA[ipart][imodule][27]=0.55;
    		ETA[ipart][imodule][28]=0.55;
   			ETA[ipart][imodule][29]=0.65;
   		 	ETA[ipart][imodule][30]=0.;
    		ETA[ipart][imodule][31]=0.;
    		ETA[ipart][imodule][32]=0.65;
    		ETA[ipart][imodule][33]=0.65;
    		ETA[ipart][imodule][34]=0.65;
    		ETA[ipart][imodule][35]=0.75;
    		ETA[ipart][imodule][36]=0.85;
    		ETA[ipart][imodule][37]=0.85;
    		ETA[ipart][imodule][38]=0.75;
    		ETA[ipart][imodule][39]=0.75;
    		ETA[ipart][imodule][40]=0.75;
    		ETA[ipart][imodule][41]=0.60;
    		ETA[ipart][imodule][42]=0.85;
    		ETA[ipart][imodule][43]=0.;
    		ETA[ipart][imodule][44]=0.60;
    		ETA[ipart][imodule][45]=0.95;
    		ETA[ipart][imodule][46]=0.95;
    		ETA[ipart][imodule][47]=0.85;

			// There are 3 layers in the calorimeter and to be able to recognize them I am setting the value for the radius R, named Rmeters here
			// Rmeters=2.5 layer A Rmeters=3.0 layer BC Rmeters=3.5 layer D
	   		Rmeters[ipart][imodule][0]=3.5;
	    	Rmeters[ipart][imodule][1]=2.5;
	    	Rmeters[ipart][imodule][2]=3.0;
	    	Rmeters[ipart][imodule][3]=3.0;
	    	Rmeters[ipart][imodule][4]=2.5;
	    	Rmeters[ipart][imodule][5]=2.5;
	    	Rmeters[ipart][imodule][6]=3.0;
	    	Rmeters[ipart][imodule][7]=3.0;
	    	Rmeters[ipart][imodule][8]=2.5;
	    	Rmeters[ipart][imodule][9]=2.5;
	    	Rmeters[ipart][imodule][10]=2.5;
	    	Rmeters[ipart][imodule][11]=3.0;
	    	Rmeters[ipart][imodule][12]=3.0;
    		Rmeters[ipart][imodule][13]=3.5;
    		Rmeters[ipart][imodule][14]=3.5;
    		Rmeters[ipart][imodule][15]=2.5;
    		Rmeters[ipart][imodule][16]=3.0;
    		Rmeters[ipart][imodule][17]=3.0;
    		Rmeters[ipart][imodule][18]=2.5;
    		Rmeters[ipart][imodule][19]=2.5;
    		Rmeters[ipart][imodule][20]=2.5;
    		Rmeters[ipart][imodule][21]=3.0;
    		Rmeters[ipart][imodule][22]=3.0;
    		Rmeters[ipart][imodule][23]=2.5;
    		Rmeters[ipart][imodule][24]=3.5;
    		Rmeters[ipart][imodule][25]=3.5;
    		Rmeters[ipart][imodule][26]=2.5;
    		Rmeters[ipart][imodule][27]=3.0;
    		Rmeters[ipart][imodule][28]=3.0;
    		Rmeters[ipart][imodule][29]=2.5;
    		Rmeters[ipart][imodule][30]=0.;
    		Rmeters[ipart][imodule][31]=0.;
    		Rmeters[ipart][imodule][32]=2.5;
    		Rmeters[ipart][imodule][33]=3.0;
    		Rmeters[ipart][imodule][34]=3.0;
    		Rmeters[ipart][imodule][35]=2.5;
    		Rmeters[ipart][imodule][36]=2.5;
    		Rmeters[ipart][imodule][37]=2.5;
    		Rmeters[ipart][imodule][38]=2.5;
    		Rmeters[ipart][imodule][39]=3.0;
    		Rmeters[ipart][imodule][40]=3.0;
    		Rmeters[ipart][imodule][41]=3.5;
    		Rmeters[ipart][imodule][42]=3.0;
    		Rmeters[ipart][imodule][43]=0.;
    		Rmeters[ipart][imodule][44]=3.5;
    		Rmeters[ipart][imodule][45]=2.5;
    		Rmeters[ipart][imodule][46]=2.5;
    		Rmeters[ipart][imodule][47]=3.0;	
    	}
	}

	//partitions LBA and LBC

	for(int ipart=2; ipart<4; ipart++){
		partner[ipart][0]=0;
		partner[ipart][1]=0;
		partner[ipart][2]=3;
		partner[ipart][3]=2;
		partner[ipart][4]=5;
		partner[ipart][5]=4;
		partner[ipart][6]=7;
		partner[ipart][7]=6;
		partner[ipart][8]=9;
		partner[ipart][9]=8;
		partner[ipart][10]=11;
		partner[ipart][11]=10;
		partner[ipart][12]=0;
		partner[ipart][13]=0;
		partner[ipart][14]=15;
		partner[ipart][15]=14;
		partner[ipart][16]=17;
		partner[ipart][17]=16;
		partner[ipart][18]=0;
		partner[ipart][19]=0;
		partner[ipart][20]=21;
		partner[ipart][21]=20;
		partner[ipart][22]=23;
		partner[ipart][23]=22;
		partner[ipart][24]=0;
		partner[ipart][25]=0;
		partner[ipart][26]=0;
		partner[ipart][27]=0;
		partner[ipart][28]=0;
		partner[ipart][29]=0;
		partner[ipart][30]=35;
		partner[ipart][31]=32;
		partner[ipart][32]=31;
		partner[ipart][33]=0;
		partner[ipart][34]=0;
		partner[ipart][35]=30;
		partner[ipart][36]=39;
		partner[ipart][37]=38;
		partner[ipart][38]=37;
		partner[ipart][39]=36;
		partner[ipart][40]=41;
		partner[ipart][41]=40;
		partner[ipart][42]=0;
		partner[ipart][43]=0;
		partner[ipart][44]=0;
		partner[ipart][45]=0;
		partner[ipart][46]=0;
		partner[ipart][47]=0;
	    for(int imodule=0; imodule<64; imodule++){
	    	ETA[ipart][imodule][0]=1.30;
   		 	ETA[ipart][imodule][1]=1.50;
   			ETA[ipart][imodule][2]=0.85;
   			ETA[ipart][imodule][3]=0.85;
   			ETA[ipart][imodule][4]=0.95;
   			ETA[ipart][imodule][5]=0.95;
   		 	ETA[ipart][imodule][6]=1.20;
   		 	ETA[ipart][imodule][7]=1.20;
    		ETA[ipart][imodule][8]=1.10;
    		ETA[ipart][imodule][9]=1.10;
    		ETA[ipart][imodule][10]=1.25;
    		ETA[ipart][imodule][11]=1.25;
    		ETA[ipart][imodule][12]=1.05;
    		ETA[ipart][imodule][13]=1.15;
    		ETA[ipart][imodule][14]=1.15;
    		ETA[ipart][imodule][15]=1.15;
    		ETA[ipart][imodule][16]=1.00;
   		 	ETA[ipart][imodule][17]=1.00;
    		ETA[ipart][imodule][18]=0.;
   	 		ETA[ipart][imodule][19]=0.;
   			ETA[ipart][imodule][20]=1.35;
   			ETA[ipart][imodule][21]=1.35;
   			ETA[ipart][imodule][22]=1.25;
   			ETA[ipart][imodule][23]=1.25;
   			ETA[ipart][imodule][24]=0.;
   			ETA[ipart][imodule][25]=0.;
   			ETA[ipart][imodule][26]=0.;
   			ETA[ipart][imodule][27]=0.;
   			ETA[ipart][imodule][28]=0.;
   			ETA[ipart][imodule][29]=0.;
   			ETA[ipart][imodule][30]=1.35;
   			ETA[ipart][imodule][31]=1.45;
   			ETA[ipart][imodule][32]=1.45;
   			ETA[ipart][imodule][33]=0.;
   			ETA[ipart][imodule][34]=0.;
   			ETA[ipart][imodule][35]=1.35;
   			ETA[ipart][imodule][36]=1.45;
    		ETA[ipart][imodule][37]=1.25;
    		ETA[ipart][imodule][38]=1.25;
    		ETA[ipart][imodule][39]=1.45;
    		ETA[ipart][imodule][40]=1.55;
    		ETA[ipart][imodule][41]=1.55;
    		ETA[ipart][imodule][42]=0.;
    		ETA[ipart][imodule][43]=0.;
    		ETA[ipart][imodule][44]=0.;
    		ETA[ipart][imodule][45]=0.;
    		ETA[ipart][imodule][46]=0.;
    		ETA[ipart][imodule][47]=0.;

			// There are 3 layers in the calorimeter and to be able to recognize them I am setting the value for the radius R, named Rmeters here
			// Rmeters=2.5 layer A Rmeters=3.0 layer B Rmeters=3.5 layer D. For the scintillators I am using Rmeters= 1.5(E4) 2.0(E3) 2.5(E2) 3.0(E1)

		    Rmeters[ipart][imodule][0]=2.0;
		    Rmeters[ipart][imodule][1]=1.5;
		    Rmeters[ipart][imodule][2]=3.5;
		    Rmeters[ipart][imodule][3]=3.5;
		    Rmeters[ipart][imodule][4]=3.0;
		    Rmeters[ipart][imodule][5]=3.0;
		    Rmeters[ipart][imodule][6]=2.5;
	   		Rmeters[ipart][imodule][7]=2.5;
		    Rmeters[ipart][imodule][8]=3.0;
		    Rmeters[ipart][imodule][9]=3.0;
		    Rmeters[ipart][imodule][10]=2.5;
		    Rmeters[ipart][imodule][11]=2.5;
		    Rmeters[ipart][imodule][12]=3.0;
		    Rmeters[ipart][imodule][13]=2.5;
		    Rmeters[ipart][imodule][14]=3.0;
		    Rmeters[ipart][imodule][15]=3.0;
		    Rmeters[ipart][imodule][16]=3.5;
		    Rmeters[ipart][imodule][17]=3.5;
		    Rmeters[ipart][imodule][18]=0.;
		    Rmeters[ipart][imodule][19]=0.;
		    Rmeters[ipart][imodule][20]=2.5;
		    Rmeters[ipart][imodule][21]=2.5;
		    Rmeters[ipart][imodule][22]=3.0;
		    Rmeters[ipart][imodule][23]=3.0;
		    Rmeters[ipart][imodule][24]=0.;
		    Rmeters[ipart][imodule][25]=0.;
		    Rmeters[ipart][imodule][26]=0.;
		    Rmeters[ipart][imodule][27]=0.;
	   		Rmeters[ipart][imodule][28]=0.;
	    	Rmeters[ipart][imodule][29]=0.;
	    	Rmeters[ipart][imodule][30]=3.0;
	    	Rmeters[ipart][imodule][31]=2.5;
	    	Rmeters[ipart][imodule][32]=2.5;
	    	Rmeters[ipart][imodule][33]=0.;
	    	Rmeters[ipart][imodule][34]=0.;
	    	Rmeters[ipart][imodule][35]=3.0;
	    	Rmeters[ipart][imodule][36]=3.0;
	    	Rmeters[ipart][imodule][37]=3.5;
	    	Rmeters[ipart][imodule][38]=3.5;
	    	Rmeters[ipart][imodule][39]=3.0;
	    	Rmeters[ipart][imodule][40]=2.5;
    		Rmeters[ipart][imodule][41]=2.5;
    		Rmeters[ipart][imodule][42]=0.;
    		Rmeters[ipart][imodule][43]=0.;
    		Rmeters[ipart][imodule][44]=0.;
    		Rmeters[ipart][imodule][45]=0.;
    		Rmeters[ipart][imodule][46]=0.;
    		Rmeters[ipart][imodule][47]=0.;
			//	
  		}
	}

	for(int jpart=0;jpart<4;jpart++){
	 for(int jmodule=0;jmodule<64;jmodule++){
	    for (int jchannel=0;jchannel<48;jchannel++){
		    // calculate the phi coordinate for each module
	        PHI[jpart][jmodule][jchannel]=jmodule*2*PI/64.;
		    // calculate the eta negative coordinate on partition 1 and 3
	        if(jpart==1)ETA[jpart][jmodule][jchannel]=ETA[jpart][jmodule][jchannel]*(-1.);
	        if(jpart==3)ETA[jpart][jmodule][jchannel]=ETA[jpart][jmodule][jchannel]*(-1.);
	    }
	  }
	}

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		pileupcell_ene = 0, pileupcell_eMF = 0;
		cell_ene = 0, cell_eMF = 0;
		//jet_pile_ene = 0, jet_pile_eMF = 0;
		//pile_z1_ene = 0, pile_z1_eMF = 0;
		//mucand_pile_ene = 0, mucand_pile_eMF = 0;
	    Long64_t ientry = LoadTree(jentry);
	    if (ientry < 0) break;
	    nb = fChain->GetEntry(jentry);   nbytes += nb;
		//cout<< "nb= "<<nb<<"   "<<"jentry= "<<jentry<<endl;

		//Made by Joao Pedro!

		for (int i=0;i<1000;i++){
			njet_ene[i] = 0;
			mjet_ene[i] = 0;
			zjet_ene[i] = 0;
			Zr_ene[i] = 0; 
			Zr_eMF[i] = 0;
			njet_eMF[i]= 0;
			mjet_eMF[i] = 0;
			zjet_eMF[i] = 0;
			ismucand_ene[i] = 0; 
			ismucand_eMF[i] = 0;
			jet_ene[i] = 0; 
			etacm_ene[i] = 0;
			phicm_ene[i] = 0; 
			jet_eMF[i] = 0; 
			etacm_eMF[i] = 0;
			phicm_eMF[i] = 0;
			Emu1_ene[i] = 0; 
			Emu2_ene[i] = 0; 
			Emu3_ene[i] = 0;
			Emu1_eMF[i] = 0; 
			Emu2_eMF[i] = 0; 
			Emu3_eMF[i] = 0;
			tfjet_ene[i] = 0;
			tfjet_eMF[i] = 0;
		}	

		// **************************************************************************
		//---> For Muons: 
		//1.8 GeV at third layer (?)
		//< 6 GeV E total
		//zr = 1
		//
		//barrel:
		//< 1.0 GeV 1st layer
		//< 3.0 GeV 2nd layer
		//< 1.5 GeV 3rd layers
		//
		//eta between 0.8 and 1.0:
		//< 1.0 GeV 1st layer
		//< 3.0 GeV 2nd layer
		//< 2.5 GeV 3rd layer
		//
		//extended:
		//< 1.0 GeV 1st layer
		//< 2.0 GeV 2nd layer
		//< 2.5 GeV 3rd layer
		//
		// ****************************************************************************

		muoncandidate_ene = 0, muoncandidate_eMF = 0;

		if (pile == 0) {
			//cout << jentry << endl;
			nopile++;
		}	

		if (pile==1) {
			pileevent++;
			//cout << jentry << endl;
		}
		cout << "jentry = " << jentry << endl;
		px_ene2 = 0, py_ene2 = 0; pz_ene2 = 0;
		px_eMF2 = 0, py_eMF2 = 0; pz_eMF2 = 0;

    	for(int part=0; part<4;part++){
    	  for(int module=0; module<64; module++) {
    	    for(int channel=0; channel<48; channel++) {
				if (tMF[part][module][channel][0] > -25 && tMF[part][module][channel][0] < 25 && pedMF[part][module][channel] < 100) {

					if (ene[part][module][channel]>300 && eMF[part][module][channel][0]>300) QFOF2xQFCOF_1->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>250 && eMF[part][module][channel][0]>300) QFOF2xQFCOF_2->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>200 && eMF[part][module][channel][0]>300) QFOF2xQFCOF_3->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>150 && eMF[part][module][channel][0]>300) QFOF2xQFCOF_4->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>300 && eMF[part][module][channel][0]>0) QFOF2xQFCOF_5->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>250 && eMF[part][module][channel][0]>0) QFOF2xQFCOF_6->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>200 && eMF[part][module][channel][0]>0) QFOF2xQFCOF_7->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);
					if (ene[part][module][channel]>150 && eMF[part][module][channel][0]>0) QFOF2xQFCOF_8->Fill(chi2[part][module][channel],chi2MF[part][module][channel]);

					if (ene[part][module][channel]>300 && eMF[part][module][channel][0]>300 && chi2MF[part][module][channel]<100) {
						if (chi2[part][module][channel]<QFcut_OF2) E_OF2_versus_E_COF_chi2->Fill(ene[part][module][channel],eMF[part][module][channel][0]);
						if (chi2[part][module][channel]<QFcut_OF2) E_OF2_chi2->Fill(ene[part][module][channel]);
						if (chi2[part][module][channel]<QFcut_OF2) E_COF_chi2->Fill(eMF[part][module][channel][0]);
					}

					if (chi2MF[part][module][channel] > 0.01 && chi2MF[part][module][channel] < QFcut){
						int scinti = 0;
						//no scintilators, scinti must to be 0
						if (channel == 0) scinti = 1;
						if (channel == 1) scinti = 1;
						if (channel == 12) scinti = 1;
						if (channel == 13) scinti = 1;
						theta_ene2 = 2*atan(exp(-1*ETA[part][module][channel]));
						theta_eMF2 = 2*atan(exp(-1*ETA[part][module][channel]));

						//Filling TProfile  
						if (ene[part][module][channel] > 0 && scinti == 0 && chi2[part][module][channel] < QFcut) {
							EnexQF2->Fill(ene[part][module][channel],chi2[part][module][channel]);
							EnexQF3->Fill(chi2[part][module][channel],ene[part][module][channel],1);
							EnexQF->Fill(ene[part][module][channel],chi2[part][module][channel],1);
							if (ene[part][module][channel] < 100000) EneXQF_100->Fill(ene[part][module][channel],chi2[part][module][channel],1);
								if (ene[part][module][channel] < 50000) EneXQF_50->Fill(ene[part][module][channel],chi2[part][module][channel],1);
							if (ene[part][module][channel] < 30000) EneXQF_30->Fill(ene[part][module][channel],chi2[part][module][channel],1);
						}
						if (eMF[part][module][channel][0] > 0 && scinti == 0 && chi2MF[part][module][channel] > 0.01 && chi2MF[part][module][channel] < QFcut) {
							EMFxQF2->Fill(eMF[part][module][channel][0],chi2MF[part][module][channel]);
							EMFxQF3->Fill(chi2MF[part][module][channel],eMF[part][module][channel][0],1);
							EMFxQF->Fill(eMF[part][module][channel][0],chi2MF[part][module][channel],1);
							if (eMF[part][module][channel][0] < 100000) EMFXQF_100->Fill(eMF[part][module][channel][0],chi2MF[part][module][channel],1);
							if (eMF[part][module][channel][0] < 50000) EMFXQF_50->Fill(eMF[part][module][channel][0],chi2MF[part][module][channel],1);
							if (eMF[part][module][channel][0] < 30000) EMFXQF_30->Fill(eMF[part][module][channel][0],chi2MF[part][module][channel],1);
						}

						//Filling Eta_Phi_E histograms, for each layer and also for all. Sum px, py and pz of all channels
						if (eMF[part][module][channel][0] > noisecut_eMF && chi2MF[part][module][channel] > 0.01 && chi2MF[part][module][channel] < QFcut){
						if (scinti == 0) {
							Eta_Phi_E_eMF->Fill(ETA[part][module][channel],PHI[part][module][channel],eMF[part][module][channel][0]);
							if (Rmeters[part][module][channel] == 2.5) Eta_Phi_E_eMF1->Fill(ETA[part][module][channel],PHI[part][module][channel],eMF[part][module][channel][0]); 
							if (Rmeters[part][module][channel] == 3.0) Eta_Phi_E_eMF2->Fill(ETA[part][module][channel],PHI[part][module][channel],eMF[part][module][channel][0]); 
							if (Rmeters[part][module][channel] == 3.5) Eta_Phi_E_eMF3->Fill(ETA[part][module][channel],PHI[part][module][channel],eMF[part][module][channel][0]); 	
							px_eMF2 += (eMF[part][module][channel][0])*(sin(theta_eMF2))*(cos(PHI[part][module][channel]));
							py_eMF2 += (eMF[part][module][channel][0])*(sin(theta_eMF2))*(sin(PHI[part][module][channel]));
							pz_eMF2 += (eMF[part][module][channel][0])*(cos(theta_eMF2));
						}
						}
     				    if (ene[part][module][channel] > noisecut_ene && chi2[part][module][channel] < QFcut) {
							if (scinti == 0) {
								Eta_Phi_E_ene->Fill(ETA[part][module][channel],PHI[part][module][channel],ene[part][module][channel]);
								if (Rmeters[part][module][channel] == 2.5) Eta_Phi_E_ene1->Fill(ETA[part][module][channel],PHI[part][module][channel],ene[part][module][channel]); 
								if (Rmeters[part][module][channel] == 3.0) Eta_Phi_E_ene2->Fill(ETA[part][module][channel],PHI[part][module][channel],ene[part][module][channel]); 
								if (Rmeters[part][module][channel] == 3.5) Eta_Phi_E_ene3->Fill(ETA[part][module][channel],PHI[part][module][channel],ene[part][module][channel]); 
								px_ene2 += (ene[part][module][channel])*(sin(theta_ene2))*(cos(PHI[part][module][channel]));
								py_ene2 += (ene[part][module][channel])*(sin(theta_ene2))*(sin(PHI[part][module][channel]));
								pz_ene2 += (ene[part][module][channel])*(cos(theta_ene2));
							}
						}
					}
				}
			}
		  }
		}
		//Pile-up Analysis
		pileup = 0;
		//cout << "Event " << jentry << " part mod ch   s0 s1 s2 s3 s4 s5 s6 " << endl;
		pileup_event[jentry] = pileup;
		//if (pileup==0) cout << "Event number " << jentry << " " << pileup_event[jentry] << endl;
		//px_ene2 = 0, py_ene2 = 0; pz_ene2 = 0;
		//px_eMF2 = 0, py_eMF2 = 0; pz_eMF2 = 0;

		pile = 0;
		pile_ev=0;
    	for(int jpart=0; jpart<4;jpart++){
    	  for(int jmodule=0; jmodule<64; jmodule++){
    	    for(int jchannel=0; jchannel<48; jchannel++){
				if (tMF[jpart][jmodule][jchannel][0] > -25 && tMF[jpart][jmodule][jchannel][0] < 25 && pedMF[jpart][jmodule][jchannel] < 100 && chi2MF[jpart][jmodule][jchannel] > 0.01 && chi2MF[jpart][jmodule][jchannel] < QFcut && chi2[jpart][jmodule][jchannel] < QFcut_OF2 && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) {
					OF2COFall++;
					if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
						OF2COFpileup++;
						pile = 1;

						if (ene[jpart][jmodule][jchannel] > noisecut_ene) {
							if (eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel] < 0.5 || eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel] > 1.5) pile_ev=1;
						}
						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 0.5 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 1.5) OF2COFup05down15++;
						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 1 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 2) OF2COFup1down2++; 
						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 2 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 3) OF2COFup2down3++; 
						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 3) OF2COFup3++;	
		 				if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 0 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 1) OF2COFup0down1++;
 						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 0 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 0.5) OF2COFup0down05++;
	 					if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 0.5 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 1) OF2COFup05down1++;
 						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > -1 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 0) OF2COFup_1down0++;
 						if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > -2 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < -1) OF2COFup_2down_1++;
					 	if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > -3 && ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < -2) OF2COFup_3down_2++;
					 	if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < -3) OF2COFdown_3++;
					}
	   		    }
	   		}
	  	  }
		}

		//Calculation pT_miss and phi_miss from channel by channel momentum
		//px_miss = -px_sum
		px_ene2 = -1*px_ene2;
		py_ene2 = -1*py_ene2;

		pT_miss_ene2 = sqrt(px_ene2*px_ene2 + py_ene2*py_ene2);		
		if (px_ene2 == 0 && py_ene2 > 0) phi_miss_ene2 = PI/2;
		if (px_ene2 == 0 && py_ene2 < 0) phi_miss_ene2 = 3*PI/2;
		if (px_ene2 > 0 && py_ene2 == 0) phi_miss_ene2 = 0;
		if (px_ene2 < 0 && py_ene2 == 0) phi_miss_ene2 = PI;		
		if (px_ene2 > 0 && py_ene2 > 0) phi_miss_ene2 = atan(py_ene2/px_ene2);
		if (px_ene2 < 0 && py_ene2 < 0) phi_miss_ene2 = atan(py_ene2/px_ene2) + PI;
		if (px_ene2 > 0 && py_ene2 < 0) phi_miss_ene2 = atan(py_ene2/px_ene2) + 2*PI;
		if (px_ene2 < 0 && py_ene2 > 0) phi_miss_ene2 = atan(py_ene2/px_ene2) + PI;
		Phi_miss_ene2->Fill(phi_miss_ene2);
		PT_miss_ene2->Fill(pT_miss_ene2);
		if (pile_ev==1) PT_miss_pileev_ene2->Fill(pT_miss_ene2);
		px_eMF2 = -1*px_eMF2;
		py_eMF2 = -1*py_eMF2;
		pT_miss_eMF2 = sqrt(px_eMF2*px_eMF2 + py_eMF2*py_eMF2);
		if (px_eMF2 == 0 && py_eMF2 > 0) phi_miss_eMF2 = PI/2;
		if (px_eMF2 == 0 && py_eMF2 < 0) phi_miss_eMF2 = 3*PI/2;
		if (px_eMF2 > 0 && py_eMF2 == 0) phi_miss_eMF2 = 0;
		if (px_eMF2 < 0 && py_eMF2 == 0) phi_miss_eMF2 = PI;
		if (px_eMF2 > 0 && py_eMF2 > 0) phi_miss_eMF2 = atan(py_eMF2/px_eMF2);
		if (px_eMF2 < 0 && py_eMF2 < 0) phi_miss_eMF2 = atan(py_eMF2/px_eMF2) + PI;
		if (px_eMF2 > 0 && py_eMF2 < 0) phi_miss_eMF2 = atan(py_eMF2/px_eMF2) + 2*PI;
		if (px_eMF2 < 0 && py_eMF2 > 0) phi_miss_eMF2 = atan(py_eMF2/px_eMF2) + PI;
		Phi_miss_eMF2->Fill(phi_miss_eMF2);
		PT_miss_eMF2->Fill(pT_miss_eMF2);
		if (pile_ev==1) PT_miss_pileev_eMF2->Fill(pT_miss_eMF2);
		Phi_miss_ene_eMF->Fill(phi_miss_ene2,phi_miss_eMF2);
		PT_miss_ene_eMF->Fill(pT_miss_ene2,pT_miss_eMF2);
		//-------------------------------------------

		//Creating a root file for each Event, with Eta x Phi distribution
		TString s;
		s = Form("%lli",jentry);

		//Writting and drawing the remarkable events
		if (jentry == 100000000000000000){

		
   	    	//for(int part=0; part<4;part++) {
   	          //for(int module=0; module<64; module++){
       		  	//for(int channel=0; channel<48; channel++){
					//if (eMF[part][module][channel][0]>noisecut_eMF && tMF[part][module][channel][0] > -25 && tMF[part][module][channel][0] < 25 && pedMF[part][module][channel] < 100 && eMF[part][module][channel][0] > 0 && chi2MF[part][module][channel] > 0.01 && chi2MF[part][module][channel] < 50){
						//cout << "part module channel s0 s1 s2 s3 s4 s5 s6 ene   eMF" << endl;
 						//cout << part << " " << module << " " << channel << " " << " " << sample[part][module][channel][0] << " " << sample[part][module][channel][1] << " " << sample[part][module][channel][2] << " " << sample[part][module][channel][3] << " " << sample[part][module][channel][4] << " " << sample[part][module][channel][5] << " " << sample[part][module][channel][6] << " " << eMF[part][module][channel][0] << " " << ene[part][module][channel] << endl;
					//}
				//}
			  //}
			//}


	      	TCanvas *eventOF21 = new TCanvas("c1","c1",900,700);
	      	Eta_Phi_E_ene->SetTitle("OF2 - Energy distribution in Eta x Phi Event " + s + " QF<" + SQFcut);	
			Eta_Phi_E_ene->Draw("LEGO");
			eventOF21->SaveAs("Remarkable_Events/Event" + s + "OF2QF" + SQFcut + "(1).jpeg","jpeg");
			eventOF21->Close();
			eventOF21->Delete();

	      	TCanvas *eventOF22 = new TCanvas("c2","c2",900,700);
	      	Eta_Phi_E_ene->SetTitle("OF2 - Energy distribution in Eta x Phi Event " + s + " QF<" + SQFcut);	
			Eta_Phi_E_ene->Draw("COLZ");
			eventOF22->SaveAs("Remarkable_Events/Event" + s + "OF2QF" + SQFcut + "(2).jpeg","jpeg");
			eventOF22->Close();
			eventOF22->Delete();

    	  	TCanvas *eventCOF1 = new TCanvas("c3","c3",900,700);
    	  	Eta_Phi_E_eMF->SetTitle("COF - Energy distribution in Eta x Phi Event " + s + " QF<" + SQFcut);	
			Eta_Phi_E_eMF->Draw("LEGO");
			eventCOF1->SaveAs("Remarkable_Events/Event" + s + "COFQF" + SQFcut + "(1).jpeg","jpeg");;
			eventCOF1->Close();
			eventCOF1->Delete();

      		TCanvas *eventCOF2 = new TCanvas("c4","c4",900,700);
      		Eta_Phi_E_eMF->SetTitle("COF - Energy distribution in Eta x Phi Event " + s + " QF<" + SQFcut);	
			Eta_Phi_E_eMF->Draw("COLZ");
			eventCOF2->SaveAs("Remarkable_Events/Event" + s + "COFQF" + SQFcut + "(2).jpeg","jpeg");;
			eventCOF2->Close();
			eventCOF2->Delete();
		}


		//Filling the array tfbin with true for bins with energy
		for (i = 0; i < 40; i++) {
			for (j = 0; j < 64; j++) {
				if (Eta_Phi_E_ene->GetBinContent(i+1,j+1) > 0) tfbin_ene[i][j] = true;
				if (Eta_Phi_E_eMF->GetBinContent(i+1,j+1) > 0) tfbin_eMF[i][j] = true;
			}
		}

		//OF2 jets

		zr_ene = 0, zi_ene = 0, zii_ene = 0, zmu_ene1 = 0, zmu_ene2 = 0, zmu_ene = 0, k_ene = 0, px_ene = 0, py_ene = 0, pz_ene = 0, eta_miss = 0, theta_ene = 0, theta_miss = 0, n_ene = 0, m_ene = 0, i = 0, j = 0, u = 0, v = 0, g = 0, zmu1_ene = 0, Emu1_ene[k_ene] = 0, Emu2_ene[k_ene] = 0, Emu3_ene[k_ene] = 0, r = 0;
		b=1;
		while (b != 0){
			Emu1_ene[k_ene] = 0, Emu2_ene[k_ene] = 0, Emu3_ene[k_ene] = 0;
			b=0;
			//Finding the bin with highest energy and call it b "n x m".
			for (i = 0; i < 40; i++) {
				for (j = 0; j < 64; j++) {
					if (Eta_Phi_E_ene->GetBinContent(i+1,j+1) > b && tfbin_ene[i][j] == true){
						b = Eta_Phi_E_ene->GetBinContent(i+1,j+1);
						n_ene = i;
						m_ene = j;
					}
				}
			}
			//Reconstructing the jets
			if (b != 0) {

				njet_ene[k_ene] = n_ene;
				mjet_ene[k_ene] = m_ene;

				jet_ene[k_ene] = Eta_Phi_E_ene->GetBinContent(n_ene+1,m_ene+1);

				Emu1_ene[k_ene] = Eta_Phi_E_ene1->GetBinContent(n_ene+1,m_ene+1);
				Emu2_ene[k_ene] = Eta_Phi_E_ene2->GetBinContent(n_ene+1,m_ene+1);
				Emu3_ene[k_ene] = Eta_Phi_E_ene3->GetBinContent(n_ene+1,m_ene+1);
				tfbin_ene[n_ene][m_ene] = false;
				//zr = "radius" of the jet, it will grow up until the energy of external ring be 0
				zr_ene = 0;
				while ((jet_ene[k_ene] - r) != 0){
					zr_ene++;
					r = jet_ene[k_ene];
					for (u = -zr_ene; u <= zr_ene; u++){
						for (v = -zr_ene; v <= zr_ene; v++){
							if(tfbin_ene[n_ene+u][m_ene+v] == true && (n_ene+u) < 40 && (m_ene+v) < 64 && (n_ene+u) > -1 && (m_ene+v) > -1){ 
								tfbin_ene[n_ene+u][m_ene+v] = false;
						 		jet_ene[k_ene] += Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
								Emu1_ene[k_ene] += Eta_Phi_E_ene1->GetBinContent(n_ene+u+1,m_ene+v+1);
								Emu2_ene[k_ene] += Eta_Phi_E_ene2->GetBinContent(n_ene+u+1,m_ene+v+1);
								Emu3_ene[k_ene] += Eta_Phi_E_ene3->GetBinContent(n_ene+u+1,m_ene+v+1);
							}					
						}
					}
				}//from while ((jet[k] - r) != 0) end the cells of the jets count

				//ENERGY IN EACH LAYER
				if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0) Elayers1_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);
				if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Elayers2_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);
				if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0) Elayers3_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);
		      	if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Elayers12_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);
				if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0) Elayers13_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);
				if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] != 0) Elayers23_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene],zr_ene);

				Zr_dist_ene->Fill(zr_ene);
				deltaR_ene = sqrt(2)*(zr_ene-1)*0.1;
				DeltaR_dist_ene->Fill(deltaR_ene);

				//Calculating the center of energy of OF2 jets		
				sum = 0;
				for (n = -zr_ene; n<=zr_ene; n++){
					sum = 0;
					for (i = -zr_ene; i <= zr_ene; i++){
						sum += Eta_Phi_E_ene->GetBinContent(njet_ene[k_ene]+n+1,mjet_ene[k_ene]+i+1);
					}
					etacm_ene[k_ene] += (njet_ene[k_ene] + n)*sum;
				}
				etacm_ene[k_ene] = etacm_ene[k_ene]/jet_ene[k_ene];
				sum = 0;
				for (n = -zr_ene; n<=zr_ene; n++){
				sum = 0;
					for (i = -zr_ene; i <= zr_ene; i++){
						sum += Eta_Phi_E_ene->GetBinContent(njet_ene[k_ene]+i+1,mjet_ene[k_ene]+n+1);
					}
					phicm_ene[k_ene] += (mjet_ene[k_ene] + n)*sum;
				}
				phicm_ene[k_ene] = phicm_ene[k_ene]/jet_ene[k_ene];
				//	

				//Calculating the sum of pT OF2 events
	
				etacm_ene[k_ene] = (etacm_ene[k_ene])*0.1-2.0;
				phicm_ene[k_ene] = (phicm_ene[k_ene])*0.1;
				eta_cm_ene->Fill(etacm_ene[k_ene]);	
				phi_cm_ene->Fill(phicm_ene[k_ene]);
				theta_ene = 2*atan(exp(-1*etacm_ene[k_ene]));
				px_ene += jet_ene[k_ene]*sin(theta_ene)*cos(phicm_ene[k_ene]);
				py_ene += jet_ene[k_ene]*sin(theta_ene)*sin(phicm_ene[k_ene]);
				pz_ene += jet_ene[k_ene]*cos(theta_ene);


				//Looking for the energy in diferent layers

				if (zr_ene == 1){
					Ejets_all3layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	
					liquid = false;
		        	if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0) liquid = true;
					if (liquid == false && (Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene]) !=0 ) Ejets_layers_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene]);
					if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] != 0) Elayersall_ene->Fill(Emu1_ene[k_ene]+Emu2_ene[k_ene]+Emu3_ene[k_ene]);		
					zmu1_ene++;
		       	}	

				//Filling the histogram of energy and counting the number of jets
				zjet_ene[k_ene] = zr_ene;
				tfjet_ene[k_ene] = true;

				//LOOKING FOR MUONS 

				if (zr_ene == 1){
					//AXION PARTICLE SEARCH-----------------------------
					//Energy just on layer 1, 2 or 3:

					//Just Energy on layer 1
	    	    	if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0) {
						Ejets_1layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	
						//Check the number of channels with energy per jet
						n_channel=0;
				    	for(int jpart=0;jpart<4;jpart++){
    					  for(int jmodule=0;jmodule<64;jmodule++){
        				    for (int jchannel=0;jchannel<48;jchannel++){
								if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.15) {
									if (ene[jpart][jmodule][jchannel]>noisecut_ene) n_channel++;
									if (n_channel == 1) Ephoton1 = ene[jpart][jmodule][jchannel];
									if (n_channel == 2) Ephoton2 = ene[jpart][jmodule][jchannel];
								}
						    }
						  }		
						}
						if (n_channel == 2) {
							//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_ene->Fill(M_particle);
							//THETA_Particle_ene->Fill(theta_particle);
						}
						N_channel_per_jet_ene->Fill(n_channel);
						N_channel_per_jet_Energy_ene->Fill(n_channel,Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
					}	

					//Just Energy on layer 2
					if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Ejets_2layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	

					//Just Energy on layer 3
	   				if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0) {
						Ejets_3layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
						//Check the number of channels with energy per jet
						n_channel=0;
						for(int jpart=0;jpart<4;jpart++){
    					  for(int jmodule=0;jmodule<64;jmodule++){
        			 	   for (int jchannel=0;jchannel<48;jchannel++){
								if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.15) {
									if (ene[jpart][jmodule][jchannel]>noisecut_ene) n_channel++;
									if (n_channel == 1) Ephoton1 = ene[jpart][jmodule][jchannel];
									if (n_channel == 2) Ephoton2 = ene[jpart][jmodule][jchannel];
								}
					 		}
						  }		
						}

						if (n_channel == 2){
			 				//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_ene->Fill(M_particle);
							//THETA_Particle_ene->Fill(theta_particle);
						}
						N_channel_per_jet_ene->Fill(n_channel);
						N_channel_per_jet_Energy_ene->Fill(n_channel,Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
					}

					//Energy just on layer 1 and 2 or 2 and 3:
	    	    	if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Ejets_12layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	
	    	    	if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] != 0) Ejets_23layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
					//----------------------------------------------------------------------------------------------------------
				}

				if (zr_ene == 2) {
					g=0;
					for (u = -1; u <= 1; u++){
						for (v = -1; v <= 1; v++){
							if (Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1) > 0) {
								g++;
								//FIND MASS 135.
								if (g==1) {
									Eta_photon1=n_ene+u+1;
									Phi_photon1=m_ene+u+1;
									Ephoton1 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
								}
								if (g==2) {
									Ephoton2 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
									Eta_photon2=n_ene+u+1;
									Phi_photon2=m_ene+u+1;
								}
							}	
						}
					}
				    if (g == 2){
				    	//cout << "ZR=1+1CELL : event " << jentry << endl;
				    	//AXION PARTICLE SEARCH-----------------------------
				    	//Energy just on layer 1, 2 or 3:
				    	//Just Enegy on layer 1:
						if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
							if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(8.125*tan(0.141421356/2));
							if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(8.125*tan(0.1/2));
							M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_ene->Fill(M_particle);
							Ejets_1layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	
							//Check the number of channels with energy per jet
							n_channel=0;
					    	for(int jpart=0;jpart<4;jpart++){
    					  	  for(int jmodule=0;jmodule<64;jmodule++){
    	    			   		for (int jchannel=0;jchannel<48;jchannel++){
									if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.15) {
										if (ene[jpart][jmodule][jchannel]>noisecut_ene) n_channel++;
										if (n_channel == 1) Ephoton1 = ene[jpart][jmodule][jchannel];
										if (n_channel == 2) Ephoton2 = ene[jpart][jmodule][jchannel];
									}
								}
							  }		
							}
							if (n_channel == 2) {
								//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
								//MASS_Particle_ene->Fill(M_particle);
								//THETA_Particle_ene->Fill(theta_particle);
							}
					
							N_channel_per_jet_ene->Fill(n_channel);
							N_channel_per_jet_Energy_ene->Fill(n_channel,Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
						}	

						//Just Energy on Layer 2
						if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Ejets_2layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);		
		
						//Just Energy on Layer 3
		    	 	    if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
				 			Ejets_3layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
							//Check the number of channels with energy per jet
							n_channel=0;
						 	for(int jpart=0;jpart<4;jpart++){
    						  for(int jmodule=0;jmodule<64;jmodule++){
    		    			    for (int jchannel=0;jchannel<48;jchannel++){
									if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.15) {
										if (ene[jpart][jmodule][jchannel]>noisecut_ene) n_channel++;
										if (n_channel == 1) Ephoton1 = ene[jpart][jmodule][jchannel];
										if (n_channel == 2) Ephoton2 = ene[jpart][jmodule][jchannel];
									}
								}
							  }		
							}

							if (n_channel == 2) {
								//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
								//MASS_Particle_ene->Fill(M_particle);
								//THETA_Particle_ene->Fill(theta_particle);
							}
							N_channel_per_jet_ene->Fill(n_channel);
							N_channel_per_jet_Energy_ene->Fill(n_channel,Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
						}
						//Energy just on layer 1 and 2 or 2 and 3:
		        		if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0) Ejets_12layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);	
		        		if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] != 0) Ejets_23layers_ene->Fill(Emu1_ene[k_ene] + Emu2_ene[k_ene] + Emu3_ene[k_ene]);
						//--------------------------------------------------
					}
				}

				//CHECK MUONS LAYER BY LAYER
				ismucand_ene[k_ene] = 0;

				la = 0, lbc = 0, ld = 0, theta = 0, ra = 0, rbc = 0, rd = 0;

				//barrel region
				if ((abs((njet_ene[k_ene])*0.1-2.0)) <= 0.8){ //condition to LB Region

					if (zr_ene == 1){

						theta = 2*atan(exp(-1*(abs((njet_ene[k_ene])*0.1-2.0)))); // try with etacm_ene[k_ene]
						la = 31.05/sin(theta);
						lbc = 84.87/sin(theta);
						ld = 37.26/sin(theta);
						ra = Emu1_ene[k_ene]/la;
						rbc = Emu2_ene[k_ene]/lbc;
						rd = Emu3_ene[k_ene]/ld;	

						if (ra>0) Emm_all_nocut_ene->Fill(ra);
						if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
						if (rd>0) Emm_all_nocut_ene->Fill(rd);	

						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0){
								zmu_ene++;
								zmu_ene1++;
								ismucand_ene[k_ene] = 1;
								Ejets_mu_cand1_ene->Fill(jet_ene[k_ene]);
								Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_ene[k_ene]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l2_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
								}
								//---------------------------------------------------------

								if (pile == 1) muoncandidate_ene++;

								for(int jpart=0;jpart<4;jpart++){
    						  	  for(int jmodule=0;jmodule<64;jmodule++){
        					   		for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.1 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.1 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.1 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.1) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_ene->Fill(jet_ene[k_ene]);
											}
										}
					 				}
					 			  }
					   			}
							}
						}
						else{	

							//AXION PARTICLE SEARCH------------------------------------
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
							if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
							if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
							if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
							//---------------------------------------------------------
						}

						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0) E123_R1_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);

						if (ra > 0) Ratio_E_l_ene->Fill(ra);
						if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
						if (rd > 0) Ratio_E_l_ene->Fill(rd);	

						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0){
							Emu_layer1_ene->Fill(Emu1_ene[k_ene]);	
							Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
							Emu_layer3_ene->Fill(Emu3_ene[k_ene]);
						}
					}

					if (zr_ene == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1) > 0) {
									g++;
									if (g==1) {
										Eta_photon1=n_ene+u+1;
										Phi_photon1=m_ene+u+1;
										Ephoton1 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
										Eta_photon2=n_ene+u+1;
										Phi_photon2=m_ene+u+1;
									}
								}	
							}
						}
		   				if (g == 2){

			  				if (Emu1_ene[k_ene] > 0 || Emu2_ene[k_ene] > 0 || Emu3_ene[k_ene] > 0){

								Emu_layer1_ene->Fill(Emu1_ene[k_ene]);
								Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
								Emu_layer3_ene->Fill(Emu3_ene[k_ene]);

								theta = 2*atan(exp(-1*(etacm_ene[k_ene]))); // try with etacm_ene[k_ene]
								la = 31.05/sin(theta);
								lbc = 53.82/sin(theta);
								ld = 68.31/sin(theta);
								ra = Emu1_ene[k_ene]/la;
								rbc = Emu2_ene[k_ene]/lbc;
								rd = Emu3_ene[k_ene]/ld;
	
								if (ra>0) Emm_all_nocut_ene->Fill(ra);
								if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
								if (rd>0) Emm_all_nocut_ene->Fill(rd);

								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_ene++;
										zmu_ene2++;
										ismucand_ene[k_ene] = 1;
										Ejets_mu_cand2_ene->Fill(jet_ene[k_ene]);
										Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);
										if (pile == 1) muoncandidate_ene++;

										//AXION PARTICLE SEARCH------------------------------------
										if (jet_ene[k_ene]>600){

											//AXION PARTICLE SEARCH-----------------------------
				    						//Energy just on layer 1, 2 or 3:
				  					  		//Just Enegy on layer 1:
											if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l1_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}
											//Just Enegy on layer 2:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
											if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer2*tan(0.141421356/2));
											if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer2*tan(0.1/2));
											M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
											MASS_Particle_l2_ene->Fill(M_particle);
											MASS_Particle_ene->Fill(M_particle);

											}
											}
								    		//Just Enegy on layer 3:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}
											if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
											if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
											if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
											if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
											if (rd >0) Emm_nomuon_l3_ene->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{

									//AXION PARTICLE SEARCH-----------------------------
						    		//Energy just on layer 1, 2 or 3:
						    		//Just Enegy on layer 1:
									if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l1_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}
									//Just Enegy on layer 2:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);

									}
									}
				    				//Just Enegy on layer 3:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}


									//AXION PARTICLE SEARCH------------------------------------
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
									//---------------------------------------------------------
								}

								if (ra > 0) Ratio_E_l_ene->Fill(ra);
								if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
								if (rd > 0) Ratio_E_l_ene->Fill(rd);
								if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0) E123_R1_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);
		  					}	
		   				}	
						if (g > 2){	
							Ejets_z1_ene->Fill(jet_ene[k_ene]);
							zii_ene++;
						}
					}	
				}

				//0.8 < eta < 1.0 region
				if ((abs((njet_ene[k_ene])*0.1-2.0)) > 0.8 && (abs((njet_ene[k_ene])*0.1-2.0)) < 1.0){ //long a and bc, extended d

					if (zr_ene == 1){

						theta = 2*atan(exp(-1*(abs((njet_ene[k_ene])*0.1-2.0)))); // try with etacm_ene[k_ene]
						la = 31.05/sin(theta);
						lbc = 84.87/sin(theta);
						ld = 68.31/sin(theta);

						ra = Emu1_ene[k_ene]/la;
						rbc = Emu2_ene[k_ene]/lbc;
						rd = Emu3_ene[k_ene]/ld;

						if (ra>0) Emm_all_nocut_ene->Fill(ra);
						if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
						if (rd>0) Emm_all_nocut_ene->Fill(rd);

						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0) {
								zmu_ene++;
								zmu_ene1++;
								ismucand_ene[k_ene] = 1;
								Ejets_mu_cand1_ene->Fill(jet_ene[k_ene]);
								Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);
								if (pile == 1) muoncandidate_ene++;

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_ene[k_ene]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);	
								}
								//---------------------------------------------------------

								for(int jpart=0;jpart<4;jpart++){ 
    						  	  for(int jmodule=0;jmodule<64;jmodule++){
        					  		for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.1 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.1 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.1 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.1) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_ene->Fill(jet_ene[k_ene]);
											}
										}
			 						}
			  					  }
			   					}
							}
						}
						else{
							//AXION PARTICLE SEARCH------------------------------------
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
							if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
							if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
							if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
							//---------------------------------------------------------
						}

						if (ra > 0) Ratio_E_l_ene->Fill(ra);
						if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
						if (rd > 0) Ratio_E_l_ene->Fill(rd);

						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0)  E123_R2_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);

						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0){
							Emu_layer1_ene->Fill(Emu1_ene[k_ene]);	
							Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
							Emu_layer3_ene->Fill(Emu3_ene[k_ene]);
						}
					}
					if (zr_ene == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1) > 0){
									g++;
									if (g==1) {
										Eta_photon1=n_ene+u+1;
										Phi_photon1=m_ene+u+1;
										Ephoton1 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
										Eta_photon2=n_ene+u+1;
										Phi_photon2=m_ene+u+1;
									}	
								}
							}
						}
						if (g == 2){	
			    		
			  				if (Emu1_ene[k_ene] > 0 || Emu2_ene[k_ene] > 0 || Emu3_ene[k_ene] > 0){

								Emu_layer1_ene->Fill(Emu1_ene[k_ene]);
								Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
								Emu_layer3_ene->Fill(Emu3_ene[k_ene]);		

								theta = 2*atan(exp(-1*(etacm_ene[k_ene]))); // try with etacm_ene[k_ene]
								la = 31.05/sin(theta);
								lbc = 53.82/sin(theta);
								ld = 68.31/sin(theta);
								ra = Emu1_ene[k_ene]/la;
								rbc = Emu2_ene[k_ene]/lbc;
								rd = Emu3_ene[k_ene]/ld;	

								if (ra>0) Emm_all_nocut_ene->Fill(ra);
								if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
								if (rd>0) Emm_all_nocut_ene->Fill(rd);
				
								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_ene++;
										zmu_ene2++;
										ismucand_ene[k_ene] = 1;
										Ejets_mu_cand2_ene->Fill(jet_ene[k_ene]);
										Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);
										if (pile == 1) muoncandidate_ene++;

										//AXION PARTICLE SEARCH-----------------------------------
										if (jet_ene[k_ene]>600){

			 					   			//Energy just on layer 1, 2 or 3:
							    			//Just Enegy on layer 1:
											if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l1_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}
											//Just Enegy on layer 2:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer2*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer2*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l2_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);								
											}
											}	
							    			//Just Enegy on layer 3:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}

											if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
											if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
											if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
											if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
											if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{

			 				   		//Energy just on layer 1, 2 or 3:
						    		//Just Enegy on layer 1:
									if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l1_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}
									//Just Enegy on layer 2:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);								
									}
									}	
					    			//Just Enegy on layer 3:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}

									//AXION PARTICLE SEARCH------------------------------------
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
									//---------------------------------------------------------
								}
								if (ra > 0) Ratio_E_l_ene->Fill(ra);
								if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
								if (rd > 0) Ratio_E_l_ene->Fill(rd);
								if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0)  E123_R2_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);
		   					}	
						}	
						if (g > 2){	
							Ejets_z1_ene->Fill(jet_ene[k_ene]);
							zii_ene++;
						}
					}	
				}

				//extended region
				if ((abs((njet_ene[k_ene])*0.1-2.0)) >= 1.0){ //condition to EB region

					if (zr_ene == 1){

						theta = 2*atan(exp(-1*(etacm_ene[k_ene]))); // try with etacm_ene[k_ene]
						la = 31.05/sin(theta);
						lbc = 53.82/sin(theta);
						ld = 68.31/sin(theta);
						ra = Emu1_ene[k_ene]/la;
						rbc = Emu2_ene[k_ene]/lbc;
						rd = Emu3_ene[k_ene]/ld;

						if (ra>0) Emm_all_nocut_ene->Fill(ra);
						if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
						if (rd>0) Emm_all_nocut_ene->Fill(rd);

						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0){
								zmu_ene++;
								zmu_ene1++;
								ismucand_ene[k_ene] = 1;
								Ejets_mu_cand1_ene->Fill(jet_ene[k_ene]);
								Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);
								if (pile == 1) muoncandidate_ene++;

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_ene[k_ene]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
								}
								//---------------------------------------------------------
								for(int jpart=0;jpart<4;jpart++){
    					 		  for(int jmodule=0;jmodule<64;jmodule++){
        				   			for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_ene[k_ene] <= ETA[jpart][jmodule][jchannel]+0.1 && etacm_ene[k_ene] >= ETA[jpart][jmodule][jchannel]-0.1 && phicm_ene[k_ene] <= PHI[jpart][jmodule][jchannel]+0.1 && phicm_ene[k_ene] >= PHI[jpart][jmodule][jchannel]-0.1) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_ene->Fill(jet_ene[k_ene]);
											}
										}
				 					}
				  			  	  }
			   					}
							}
						}
						else{
							//AXION PARTICLE SEARCH
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
							if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
							if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
							if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
						}
						//----------------------------------------------------------
						if (ra > 0) Ratio_E_l_ene->Fill(ra);
						if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
						if (rd > 0) Ratio_E_l_ene->Fill(rd);
						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0)  E123_R3_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);
						if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0){
							Emu_layer1_ene->Fill(Emu1_ene[k_ene]);	
							Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
							Emu_layer3_ene->Fill(Emu3_ene[k_ene]);
						}
					}
					if (zr_ene == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1) > 0) {
									g++;	
									if (g==1) {
										Eta_photon1=n_ene+u+1;
										Phi_photon1=m_ene+u+1;
										Ephoton1 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_ene->GetBinContent(n_ene+u+1,m_ene+v+1);
										Eta_photon2=n_ene+u+1;
										Phi_photon2=m_ene+u+1;
									}					
								}
							}
						}
		  		 		if (g == 2){
		  					if (Emu1_ene[k_ene] > 0 || Emu2_ene[k_ene] > 0 || Emu3_ene[k_ene] > 0){	

								Emu_layer1_ene->Fill(Emu1_ene[k_ene]);
								Emu_layer2_ene->Fill(Emu2_ene[k_ene]);
								Emu_layer3_ene->Fill(Emu3_ene[k_ene]);

								theta = 2*atan(exp(-1*(etacm_ene[k_ene]))); // try with etacm_ene[k_ene]
								la = 31.05/sin(theta);
								lbc = 53.82/sin(theta);
								ld = 68.31/sin(theta);
								ra = Emu1_ene[k_ene]/la;
								rbc = Emu2_ene[k_ene]/lbc;
								rd = Emu3_ene[k_ene]/ld;

								if (ra>0) Emm_all_nocut_ene->Fill(ra);
								if (rbc>0) Emm_all_nocut_ene->Fill(rbc);
								if (rd>0) Emm_all_nocut_ene->Fill(rd);
								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_ene++;
										zmu_ene2++;		
										ismucand_ene[k_ene] = 1;
										Ejets_mu_cand2_ene->Fill(jet_ene[k_ene]);
										Ejets_mu_cand_ene->Fill(jet_ene[k_ene]);
										if (pile == 1) muoncandidate_ene++;	

										//AXION PARTICLE SEARCH------------------------------------
										if (jet_ene[k_ene]>600){

				 					   		//AXION PARTICLE SEARCH-----------------------------
				 				   			//Energy just on layer 1, 2 or 3:
				    						//Just Enegy on layer 1:
											if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l1_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}
											//Just Enegy on layer 2:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer2*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer2*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l2_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}
				    						//Just Enegy on layer 3:
											if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_ene->Fill(M_particle);
												MASS_Particle_ene->Fill(M_particle);
											}
											}

											if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);	
											if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
											if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
											if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
											if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{

						    		//AXION PARTICLE SEARCH-----------------------------
				    				//Energy just on layer 1, 2 or 3:
				    				//Just Enegy on layer 1:
									if (Emu1_ene[k_ene] != 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l1_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}
									//Just Enegy on layer 2:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] != 0 && Emu3_ene[k_ene] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}
						    		//Just Enegy on layer 3:
									if (Emu1_ene[k_ene] == 0 && Emu2_ene[k_ene] == 0 && Emu3_ene[k_ene] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_ene->Fill(M_particle);
										MASS_Particle_ene->Fill(M_particle);
									}
									}


									if (ra != 0 || rbc != 0) Emm_nomuon_l12_ene->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_ene->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_ene->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_ene->Fill(Emu1_ene[k_ene]);
									if (rbc>0) E_nomuon_l2_ene->Fill(Emu2_ene[k_ene]);
									if (rd>0)  E_nomuon_l3_ene->Fill(Emu3_ene[k_ene]);
									if (ra>0)  Emm_nomuon_l1_ene->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_ene->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_ene->Fill(rd);
									//---------------------------------------------------------
								}

								if (ra > 0) Ratio_E_l_ene->Fill(ra);
								if (rbc > 0) Ratio_E_l_ene->Fill(rbc);
								if (rd > 0) Ratio_E_l_ene->Fill(rd);

								if (Emu1_ene[k_ene] > 0 && Emu2_ene[k_ene] > 0 && Emu3_ene[k_ene] > 0)  E123_R3_ene->Fill(Emu1_ene[k_ene],Emu2_ene[k_ene],Emu3_ene[k_ene]);
						    }	
		   				}	
						if (g > 2){	
							Ejets_z1_ene->Fill(jet_ene[k_ene]);
							zii_ene++;
						}
					}	
				}

				//end muons 

				Zr_ene[k_ene]=zr_ene;

				if ( zr_ene == 1) Ejets_zequal1_ene->Fill(jet_ene[k_ene]);
				if ( zr_ene > 1 ) {
					Ejets_z1_ene->Fill(jet_ene[k_ene]);
					zi_ene++;
				}
				if ( zr_ene > 2) {
					Ejets_z2_ene->Fill(jet_ene[k_ene]);
					zii_ene++;
				}

				Ejets_ene->Fill(jet_ene[k_ene]);
				k_ene++;
			}//from if (b != 0)
		}//from while (b != 0) end the jets count

		px_ene = -1*px_ene;
		py_ene = -1*py_ene;

		pT_miss_ene1 = sqrt(px_ene*px_ene + py_ene*py_ene);

		if (px_ene == 0 && py_ene > 0) phi_miss_ene1 = PI/2;
		if (px_ene == 0 && py_ene < 0) phi_miss_ene1 = 3*PI/2;
		if (px_ene > 0 && py_ene == 0) phi_miss_ene1 = 0;
		if (px_ene < 0 && py_ene == 0) phi_miss_ene1 = PI;		
		if (px_ene > 0 && py_ene > 0) phi_miss_ene1 = atan(py_ene/px_ene);
		if (px_ene < 0 && py_ene < 0) phi_miss_ene1 = atan(py_ene/px_ene) + PI;
		if (px_ene > 0 && py_ene < 0) phi_miss_ene1 = atan(py_ene/px_ene) + 2*PI;
		if (px_ene < 0 && py_ene > 0) phi_miss_ene1 = atan(py_ene/px_ene) + PI;

		Phi_miss_ene->Fill(phi_miss_ene1);
		PT_miss_ene->Fill(pT_miss_ene1);
		eta_miss_ene->Fill(eta_miss);
		Njets_per_event_ene->Fill(k_ene);
		Njets_per_event_z1_ene->Fill(zi_ene);
		Njets_per_event_zequal1_ene->Fill(zmu1_ene);
		Njets_per_event_z2_ene->Fill(zii_ene);
		Nmu_per_event_ene->Fill(zmu_ene);
		Nmu_per_event1_ene->Fill(zmu_ene1);
		Nmu_per_event2_ene->Fill(zmu_ene2);

		//---------------------------------------------------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------------------------------------------------

		//COF
		zr_eMF=0, zi_eMF=0, zii_eMF=0, zmu_eMF1 = 0, zmu_eMF2 = 0, zmu_eMF = 0, k_eMF=0, px_eMF=0, py_eMF=0, pz_eMF=0, theta_eMF=0, eta_miss=0, theta_miss=0, zmu1_eMF=0, Emu1_eMF[k_eMF] = 0, Emu2_eMF[k_eMF] = 0, Emu3_eMF[k_eMF] = 0; 
		b=1;

		while (b != 0){
			Emu1_eMF[k_eMF] = 0, Emu2_eMF[k_eMF] = 0, Emu3_eMF[k_eMF] = 0;
			b=0;
			//Finding the bin with highest energy and call it b "n x m".
			for (i = 0; i < 40; i++) {
				for (j = 0; j < 64; j++) {
					if (Eta_Phi_E_eMF->GetBinContent(i+1,j+1) > b && tfbin_eMF[i][j] == true){
						b = Eta_Phi_E_eMF->GetBinContent(i+1,j+1);
						n_eMF = i;
						m_eMF = j;
					}
				}
			}

			if (b != 0){

				njet_eMF[k_eMF] = n_eMF;
				mjet_eMF[k_eMF] = m_eMF;

				jet_eMF[k_eMF] = Eta_Phi_E_eMF->GetBinContent(n_eMF+1,m_eMF+1);
				Emu1_eMF[k_eMF] = Eta_Phi_E_eMF1->GetBinContent(n_eMF+1,m_eMF+1);
				Emu2_eMF[k_eMF] = Eta_Phi_E_eMF2->GetBinContent(n_eMF+1,m_eMF+1);
				Emu3_eMF[k_eMF] = Eta_Phi_E_eMF3->GetBinContent(n_eMF+1,m_eMF+1);

			    tfbin_eMF[n_eMF][m_eMF] = false;
				//zr = "radius" of the jet, it will grow up until the energy of external ring be 0		
				zr_eMF = 0;
				while ((jet_eMF[k_eMF] - r) != 0){
					zr_eMF++;
					r = jet_eMF[k_eMF];
					for (u = -zr_eMF; u <= zr_eMF; u++){
						for (v = -zr_eMF; v <= zr_eMF; v++){
							if(tfbin_eMF[n_eMF+u][m_eMF+v] == true && (n_eMF+u) < 40 && (m_eMF+v) < 64 && (n_eMF+u) > -1 && (m_eMF+v) > -1){ 
								tfbin_eMF[n_eMF+u][m_eMF+v] = false;
							 	jet_eMF[k_eMF] += Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
							 	Emu1_eMF[k_eMF] += Eta_Phi_E_eMF1->GetBinContent(n_eMF+u+1,m_eMF+v+1);
							 	Emu2_eMF[k_eMF] += Eta_Phi_E_eMF2->GetBinContent(n_eMF+u+1,m_eMF+v+1);
							 	Emu3_eMF[k_eMF] += Eta_Phi_E_eMF3->GetBinContent(n_eMF+u+1,m_eMF+v+1);
							}					
						}
					}
				}//end a jet count
			
				Zr_dist_eMF->Fill(zr_eMF);
				deltaR_eMF = (sqrt(2))*(zr_eMF-1)*0.1;
				DeltaR_dist_eMF->Fill(deltaR_eMF);

				//ENERGY IN EACH LAYER

				if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0) Elayers1_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);
				if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Elayers2_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);
				if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0) Elayers3_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);
	    	  	if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Elayers12_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);
				if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0) Elayers13_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);
				if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] != 0) Elayers23_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF],zr_eMF);

				//Calculating the center of energy of COF jets		
				sum = 0;
				for (n = -zr_eMF; n<=zr_eMF; n++){
					sum = 0;
					for (i = -zr_eMF; i <= zr_eMF; i++){
						sum += Eta_Phi_E_eMF->GetBinContent(njet_eMF[k_eMF]+1+n,mjet_eMF[k_eMF]+1+i);
					}
	
					etacm_eMF[k_eMF] += (njet_eMF[k_eMF] + n)*sum;
				}
				etacm_eMF[k_eMF] = etacm_eMF[k_eMF]/jet_eMF[k_eMF];
	
				sum = 0;
				for (n = -zr_eMF; n<=zr_eMF; n++){
					sum = 0;
					for (i = -zr_eMF; i <= zr_eMF; i++){
						sum += Eta_Phi_E_eMF->GetBinContent(njet_eMF[k_eMF]+1+i,mjet_eMF[k_eMF]+1+n);
	    	    	}
					phicm_eMF[k_eMF] += (mjet_eMF[k_eMF] + n)*sum;
				}
				phicm_eMF[k_eMF] = phicm_eMF[k_eMF]/jet_eMF[k_eMF];
				//

				//Calculating the sum of pT 

				etacm_eMF[k_eMF] = (etacm_eMF[k_eMF])*0.1-2.0;
				phicm_eMF[k_eMF] = (phicm_eMF[k_eMF])*0.1;

				eta_cm_eMF->Fill(etacm_eMF[k_eMF]);	
				phi_cm_eMF->Fill(phicm_eMF[k_eMF]);

				theta_eMF = 2*atan(exp(-1*etacm_eMF[k_eMF]));

				px_eMF += jet_eMF[k_eMF]*sin(theta_eMF)*cos(phicm_eMF[k_eMF]);
				py_eMF += jet_eMF[k_eMF]*sin(theta_eMF)*sin(phicm_eMF[k_eMF]);
				pz_eMF += jet_eMF[k_eMF]*cos(theta_eMF);

				//Looking to the energy in diferent layers
			
				if (zr_eMF == 1){
					Ejets_all3layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	
					liquid = false;
		        	if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0) liquid = true;
					if (liquid == false && (Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF]) !=0 ) Ejets_layers_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF]);
					if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] != 0) Elayersall_eMF->Fill(Emu1_eMF[k_eMF]+Emu2_eMF[k_eMF]+Emu3_eMF[k_eMF]);		
					zmu1_eMF++;
		       	}	

				//Filling the histogram of energy and counting the number of jets
				zjet_eMF[k_eMF] = zr_eMF;
				tfjet_eMF[k_eMF] = true;

				//LOOKING FOR MUONS 

				la = 0, lbc = 0, ld = 0, theta = 0, ra = 0, rbc = 0, rd = 0;

				ismucand_eMF[k_eMF] = 0;

				if (zr_eMF == 1){
					//AXION PARTICLE SEARCH-----------------------------
					//Energy just on layer 1, 2 or 3:
					//JUST LAYER 1
		    	    if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0) {
						Ejets_1layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	
						//Check the number of channels with energy per jet
						n_channel=0;
						for(int jpart=0;jpart<4;jpart++){
    					  for(int jmodule=0;jmodule<64;jmodule++){
    	    			   for (int jchannel=0;jchannel<48;jchannel++){
								if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.15) {
									if (eMF[jpart][jmodule][jchannel][0]>noisecut_eMF) n_channel++;
									if (n_channel == 1) Ephoton1 = eMF[jpart][jmodule][jchannel][0];
									if (n_channel == 2) Ephoton2 = eMF[jpart][jmodule][jchannel][0];
								}
							}
						  }		
						}
						if (n_channel == 2) {
							//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_eMF->Fill(M_particle);
							//THETA_Particle_eMF->Fill(theta_particle);
						}

						N_channel_per_jet_eMF->Fill(n_channel);
						N_channel_per_jet_Energy_eMF->Fill(n_channel,Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
					}

					//JUST LAYER 2
	    	  	    if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Ejets_2layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);			

	    	    	//JUST LAYER 3
		        	if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
			 			Ejets_3layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
						//Check the number of channels with energy per jet
						n_channel=0;
						for(int jpart=0;jpart<4;jpart++){
    					  for(int jmodule=0;jmodule<64;jmodule++){
        				   for (int jchannel=0;jchannel<48;jchannel++){        			
								if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.15) {
									if (eMF[jpart][jmodule][jchannel][0]>noisecut_eMF) n_channel++;
									if (n_channel == 1) Ephoton1 = eMF[jpart][jmodule][jchannel][0];
									if (n_channel == 2) Ephoton2 = eMF[jpart][jmodule][jchannel][0];
								}
							}
						  }		
						}

						if (n_channel == 2) {
							//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_eMF->Fill(M_particle);
							//THETA_Particle_eMF->Fill(theta_particle);
						}

						N_channel_per_jet_eMF->Fill(n_channel);
						N_channel_per_jet_Energy_eMF->Fill(n_channel,Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
					}

					//Energy just on layer 1 and 2 or 2 and 3:
		    	    if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Ejets_12layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	
	   			   	if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] != 0) Ejets_23layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
					//-------------------------------------------------
				}
				if (zr_eMF == 2) {
					g=0;
					for (u = -1; u <= 1; u++){
						for (v = -1; v <= 1; v++){
							if (Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1) > 0) {
								g++;
								//FIND MASS 135.
								if (g==1) {
									Eta_photon1=n_eMF+u+1;
									Phi_photon1=m_eMF+u+1;
									Ephoton1 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
								}
								if (g==2) {
									Ephoton2 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
									Eta_photon2=n_eMF+u+1;
									Phi_photon2=m_eMF+u+1;
								}
							}
						}	
					}
		  			if (g == 2){
						if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0) {
							if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(8.125*tan(0.141421356/2));
							if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(8.125*tan(0.1/2));
							M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
							//MASS_Particle_eMF->Fill(M_particle);
						}

						//AXION PARTICLE SEARCH-----------------------------
						//Energy just on layer 1:
				        if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0) {
							Ejets_1layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	
							//Check the number of channels with energy per jet
							n_channel=0;

							for(int jpart=0;jpart<4;jpart++){
    				 		  for(int jmodule=0;jmodule<64;jmodule++){
        			 	  		for (int jchannel=0;jchannel<48;jchannel++){
									if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.15) {
										if (eMF[jpart][jmodule][jchannel][0]>noisecut_eMF) n_channel++;
										if (n_channel == 1) Ephoton1 = eMF[jpart][jmodule][jchannel][0];
										if (n_channel == 2) Ephoton2 = eMF[jpart][jmodule][jchannel][0];
									}
								}
							  }		
							}	

							if (n_channel == 2) {
								//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
								//MASS_Particle_eMF->Fill(M_particle);
								//THETA_Particle_eMF->Fill(theta_particle);
							}	

							N_channel_per_jet_eMF->Fill(n_channel);
							N_channel_per_jet_Energy_eMF->Fill(n_channel,Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
						}	

						//Energy just on layer 2
			   		    if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Ejets_2layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);				

				        //Energy just on layer 3
				        if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
							Ejets_3layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	

							//Check the number of channels with energy per jet
							n_channel=0;
							for(int jpart=0;jpart<4;jpart++){
    						  for(int jmodule=0;jmodule<64;jmodule++){
       							for (int jchannel=0;jchannel<48;jchannel++){
									if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.15) {
										if (eMF[jpart][jmodule][jchannel][0]>noisecut_eMF) n_channel++;
										if (n_channel == 1) Ephoton1 = eMF[jpart][jmodule][jchannel][0];
										if (n_channel == 2) Ephoton2 = eMF[jpart][jmodule][jchannel][0];
									}
								}
				  			  }		
							}

							if (n_channel == 2) {
								//M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
								//MASS_Particle_eMF->Fill(M_particle);
								//THETA_Particle_eMF->Fill(theta_particle);
							}

							N_channel_per_jet_eMF->Fill(n_channel);
							N_channel_per_jet_Energy_eMF->Fill(n_channel,Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
						}
						//Energy just on layer 1 and 2 or 2 and 3:
					    if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0) Ejets_12layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);	
				    	if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] != 0) Ejets_23layers_eMF->Fill(Emu1_eMF[k_eMF] + Emu2_eMF[k_eMF] + Emu3_eMF[k_eMF]);
						//--------------------------------------------------
					}
				}

				//barrel region
				if ((abs((njet_eMF[k_eMF])*0.1-2.0)) <= 0.8){ //condition to LB region

					if (zr_eMF == 1){

						theta = 2*atan(exp(-1*(abs((njet_eMF[k_eMF])*0.1-2.0)))); // try with etacm_eMF[k_eMF]
						la = 31.05/sin(theta);
						lbc = 84.87/sin(theta);
						ld = 37.26/sin(theta);

						ra = Emu1_eMF[k_eMF]/la;
						rbc = Emu2_eMF[k_eMF]/lbc;
						rd = Emu3_eMF[k_eMF]/ld;

						if (ra>0) Emm_all_nocut_eMF->Fill(ra);
						if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
						if (rd>0) Emm_all_nocut_eMF->Fill(rd);

						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0){
								zmu_eMF++;
								zmu_eMF1++;
								ismucand_eMF[k_eMF] = 1;
								Ejets_mu_cand1_eMF->Fill(jet_eMF[k_eMF]);
								Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
								if (pile == 1) muoncandidate_eMF++;

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_eMF[k_eMF]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
								}
								//---------------------------------------------------------

								for(int jpart=0;jpart<4;jpart++){
    							  for(int jmodule=0;jmodule<64;jmodule++){
        							for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.15 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.15 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.15 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.15) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_eMF->Fill(jet_eMF[k_eMF]);
											}
										}
									}
			 					  }
								}
							}
						}
						else{

							//AXION PARTICLE SEARCH------------------------------------
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
							if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
							if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
							if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
							//---------------------------------------------------------
						}

						if (ra > 0) Ratio_E_l_eMF->Fill(ra);
						if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
						if (rd > 0) Ratio_E_l_eMF->Fill(rd);

						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0) E123_R1_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0){
							Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);	
							Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
							Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);
						}
					}

					if (zr_eMF == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1) > 0){ 
									g++;
									if (g==1) {
										Eta_photon1=n_eMF+u+1;
										Phi_photon1=m_eMF+u+1;
										Ephoton1 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
										Eta_photon2=n_eMF+u+1;
										Phi_photon2=m_eMF+u+1;
									}											
								}
							}
						}
		 			    if (g == 2){
	 		    	
						    if (Emu1_eMF[k_eMF] > 0 || Emu2_eMF[k_eMF] > 0 || Emu3_eMF[k_eMF] > 0){

								Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);
								Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
								Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);

								theta = 2*atan(exp(-1*(etacm_eMF[k_eMF]))); // try with etacm_eMF[k_eMF]
								la = 31.05/sin(theta);
								lbc = 84.87/sin(theta);
								ld = 37.26/sin(theta);
								ra = Emu1_eMF[k_eMF]/la;
								rbc = Emu2_eMF[k_eMF]/lbc;
								rd = Emu3_eMF[k_eMF]/ld;

								if (ra>0) Emm_all_nocut_eMF->Fill(ra);
								if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
								if (rd>0) Emm_all_nocut_eMF->Fill(rd);

	
								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_eMF++;
										zmu_eMF2++;
										ismucand_eMF[k_eMF] = 1;
										Ejets_mu_cand2_eMF->Fill(jet_eMF[k_eMF]);
										Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
										if (pile == 1) muoncandidate_eMF++;

										//AXION PARTICLE SEARCH------------------------------------
										if (jet_eMF[k_eMF]>600){

						    				//AXION PARTICLE SEARCH-----------------------------
								    		//Energy just on layer 1, 2 or 3:
								    		//Just Enegy on layer 1:
											if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l1_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);
											}
											}
											//Just Enegy on layer 2:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer2*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer2*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l2_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);								
											}
											}
								    		//Just Enegy on layer 3:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);
											}
											}	



											if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
											if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
											if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
											if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
											if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{


		    						//AXION PARTICLE SEARCH-----------------------------
						    		//Energy just on layer 1, 2 or 3:
						    		//Just Enegy on layer 1:
									if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l1_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);
									}
									}
									//Just Enegy on layer 2:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);								
									}
									}
						    		//Just Enegy on layer 3:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrlonglayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrlonglayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);
									}
									}	


									//AXION PARTICLE SEARCH------------------------------------
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
									//---------------------------------------------------------
								}

								if (ra > 0) Ratio_E_l_eMF->Fill(ra);
								if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
								if (rd > 0) Ratio_E_l_eMF->Fill(rd);
								if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0) E123_R1_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
						   }	
		   				}	
						if (g > 2){	
							Ejets_z1_eMF->Fill(jet_eMF[k_eMF]);
							zii_eMF++;
						}
					}	
				}


				//0.8 < eta < 1.0 region
				if ((abs((njet_eMF[k_eMF])*0.1-2.0)) > 0.8 && (abs((njet_eMF[k_eMF])*0.1-2.0)) < 1.0){ //long a and bc, extended d

					if (zr_eMF == 1){

						theta = 2*atan(exp(-1*(abs((njet_eMF[k_eMF])*0.1-2.0)))); // try with etacm_eMF[k_eMF]
						la = 31.05/sin(theta);
						lbc = 84.87/sin(theta);
						ld = 68.31/sin(theta);
						ra = Emu1_eMF[k_eMF]/la;
						rbc = Emu2_eMF[k_eMF]/lbc;
						rd = Emu3_eMF[k_eMF]/ld;

						if (ra>0) Emm_all_nocut_eMF->Fill(ra);
						if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
						if (rd>0) Emm_all_nocut_eMF->Fill(rd);


						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0){
								zmu_eMF++;
								zmu_eMF1++;
								ismucand_eMF[k_eMF] = 1;
								Ejets_mu_cand1_eMF->Fill(jet_eMF[k_eMF]);
								Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
								if (pile == 1) muoncandidate_eMF++;

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_eMF[k_eMF]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
								}
								//---------------------------------------------------------

			 					for(int jpart=0;jpart<4;jpart++){
    							  for(int jmodule=0;jmodule<64;jmodule++){
        		    				for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.1 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.1 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.1 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.1) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_eMF->Fill(jet_eMF[k_eMF]);
											}
										}
									}
			  					  }
			    				}
							}
						}
						else{
							//AXION PARTICLE SEARCH------------------------------------
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
							if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
							if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
							if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
							//---------------------------------------------------------
						}
						if (ra > 0) Ratio_E_l_eMF->Fill(ra);
						if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
						if (rd > 0) Ratio_E_l_eMF->Fill(rd);

						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0) E123_R2_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0){
							Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);	
							Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
							Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);
						}
					}
					if (zr_eMF == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1) > 0){
									g++;
									if (g==1) {
										Eta_photon1=n_eMF+u+1;
										Phi_photon1=m_eMF+u+1;
										Ephoton1 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
										Eta_photon2=n_eMF+u+1;
										Phi_photon2=m_eMF+u+1;
									}																		
								}
							}
						}
				   		if (g == 2){
					    	
				 		    if (Emu1_eMF[k_eMF] > 0 || Emu2_eMF[k_eMF] > 0 || Emu3_eMF[k_eMF] > 0){
								Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);
								Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
								Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);	

								theta = 2*atan(exp(-1*(etacm_eMF[k_eMF]))); // try with etacm_eMF[k_eMF]
								la = 31.05/sin(theta);
								lbc = 84.87/sin(theta);
								ld = 68.31/sin(theta);
								ra = Emu1_eMF[k_eMF]/la;
								rbc = Emu2_eMF[k_eMF]/lbc;
								rd = Emu3_eMF[k_eMF]/ld;	

								if (ra>0) Emm_all_nocut_eMF->Fill(ra);
								if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
								if (rd>0) Emm_all_nocut_eMF->Fill(rd);

								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_eMF++;
										zmu_eMF2++;
										Ejets_mu_cand2_eMF->Fill(jet_eMF[k_eMF]);
										Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
										ismucand_eMF[k_eMF] = 1;
										if (pile == 1) muoncandidate_eMF++;

										//AXION PARTICLE SEARCH------------------------------------
										if (jet_eMF[k_eMF]>600){

		    								//AXION PARTICLE SEARCH-----------------------------
			    							//Energy just on layer 1, 2 or 3:
			    							//Just Enegy on layer 1:
											if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										 		MASS_Particle_l1_eMF->Fill(M_particle);
										 		MASS_Particle_eMF->Fill(M_particle);
										 	}	
											}
											//Just Enegy on layer 2:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer2*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer2*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l2_eMF->Fill(M_particle);
										 		MASS_Particle_eMF->Fill(M_particle);								
											}
											}
			    							//Just Enegy on layer 3:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_eMF->Fill(M_particle);
										 		MASS_Particle_eMF->Fill(M_particle);
											}
											}




											if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
											if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
											if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
											if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
											if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{

		    						//AXION PARTICLE SEARCH-----------------------------
					    			//Energy just on layer 1, 2 or 3:
					    			//Just Enegy on layer 1:
									if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
								 		MASS_Particle_l1_eMF->Fill(M_particle);
								 		MASS_Particle_eMF->Fill(M_particle);
								 	}	
									}
									//Just Enegy on layer 2:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_eMF->Fill(M_particle);
								 		MASS_Particle_eMF->Fill(M_particle);								
									}
									}
			    					//Just Enegy on layer 3:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrmediumlayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrmediumlayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_eMF->Fill(M_particle);
								 		MASS_Particle_eMF->Fill(M_particle);
									}
									}




									//AXION PARTICLE SEARCH------------------------------------
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
									//---------------------------------------------------------
								}
	
								if (ra > 0) Ratio_E_l_eMF->Fill(ra);
								if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
								if (rd > 0) Ratio_E_l_eMF->Fill(rd);
								if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0)E123_R2_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
				    		}	
				    	}	
						if (g > 2){	
							Ejets_z1_eMF->Fill(jet_eMF[k_eMF]);
							zii_eMF++;
						}
					}	
				}

				//extended region
				if ((abs((njet_eMF[k_eMF])*0.1-2.0)) >= 1.0){ //condition to EB region

					if (zr_eMF == 1){

						theta = 2*atan(exp(-1*(etacm_eMF[k_eMF]))); // try with etacm_eMF[k_eMF]
						la = 31.05/sin(theta);
						lbc = 53.82/sin(theta);
						ld = 68.31/sin(theta);
						ra = Emu1_eMF[k_eMF]/la;
						rbc = Emu2_eMF[k_eMF]/lbc;
						rd = Emu3_eMF[k_eMF]/ld;

						if (ra>0) Emm_all_nocut_eMF->Fill(ra);
						if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
						if (rd>0) Emm_all_nocut_eMF->Fill(rd);


						if (ra < 15. && rbc < 15. && rd < 15.){
							if (ra > 0 || rbc > 0 || rd > 0){
								zmu_eMF++;
								zmu_eMF1++;
								ismucand_eMF[k_eMF] = 1;
								Ejets_mu_cand1_eMF->Fill(jet_eMF[k_eMF]);
								Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
								if (pile == 1) muoncandidate_eMF++;

								//AXION PARTICLE SEARCH------------------------------------
								if (jet_eMF[k_eMF]>600){
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
								}
								//---------------------------------------------------------
								for(int jpart=0;jpart<4;jpart++){
    							  for(int jmodule=0;jmodule<64;jmodule++){
        			    			for (int jchannel=0;jchannel<48;jchannel++){
										if (etacm_eMF[k_eMF] <= ETA[jpart][jmodule][jchannel]+0.1 && etacm_eMF[k_eMF] >= ETA[jpart][jmodule][jchannel]-0.1 && phicm_eMF[k_eMF] <= PHI[jpart][jmodule][jchannel]+0.1 && phicm_eMF[k_eMF] >= PHI[jpart][jmodule][jchannel]-0.1) {
											if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
												Ejets_pileup_muons_eMF->Fill(jet_eMF[k_eMF]);
											}
										}
			 						}
								  }
			   					}	
							}
						}
						else{
							//AXION PARTICLE SEARCH------------------------------------
							if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
							if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
							if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
							if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
							if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
							if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
							if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
							if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
							if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
							//---------------------------------------------------------
						}
						if (ra > 0) Ratio_E_l_eMF->Fill(ra);
						if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
						if (rd > 0) Ratio_E_l_eMF->Fill(rd);

						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0)E123_R3_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
						if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0){
							Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);	
							Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
							Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);
						}
					}
					if (zr_eMF == 2) {
						g=0;
						for (u = -1; u <= 1; u++){
							for (v = -1; v <= 1; v++){
								if (Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1) > 0){
									g++;	
									if (g==1) {
										Eta_photon1=n_eMF+u+1;
										Phi_photon1=m_eMF+u+1;
										Ephoton1 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
									}
									if (g==2) {
										Ephoton2 = Eta_Phi_E_eMF->GetBinContent(n_eMF+u+1,m_eMF+v+1);
										Eta_photon2=n_eMF+u+1;
										Phi_photon2=m_eMF+u+1;
									}						
								}
							}
						}
			   			if (g == 2){

		   					if (Emu1_eMF[k_eMF] > 0 || Emu2_eMF[k_eMF] > 0 || Emu3_eMF[k_eMF] > 0){
								Emu_layer1_eMF->Fill(Emu1_eMF[k_eMF]);
								Emu_layer2_eMF->Fill(Emu2_eMF[k_eMF]);
								Emu_layer3_eMF->Fill(Emu3_eMF[k_eMF]);

								theta = 2*atan(exp(-1*(etacm_eMF[k_eMF]))); // try with etacm_eMF[k_eMF]
								la = 31.05/sin(theta);
								lbc = 53.82/sin(theta);
								ld = 68.31/sin(theta);
								ra = Emu1_eMF[k_eMF]/la;
								rbc = Emu2_eMF[k_eMF]/lbc;
								rd = Emu3_eMF[k_eMF]/ld;

								if (ra>0) Emm_all_nocut_eMF->Fill(ra);
								if (rbc>0) Emm_all_nocut_eMF->Fill(rbc);
								if (rd>0) Emm_all_nocut_eMF->Fill(rd);

								if (ra < 15. && rbc < 15. && rd < 15.){
									if (ra > 0 || rbc > 0 || rd > 0){
										zmu_eMF++;
										zmu_eMF2++;
										ismucand_eMF[k_eMF] = 1;
										Ejets_mu_cand2_eMF->Fill(jet_eMF[k_eMF]);
										Ejets_mu_cand_eMF->Fill(jet_eMF[k_eMF]);
										if (pile == 1) muoncandidate_eMF++;

										//AXION PARTICLE SEARCH------------------------------------
										if (jet_eMF[k_eMF]>600){
		    
		    								//AXION PARTICLE SEARCH-----------------------------
								    		//Energy just on layer 1, 2 or 3:
								    		//Just Enegy on layer 1:
											if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																								
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer1*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer1*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l1_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);
											}
											}	
											//Just Enegy on layer 2:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																	
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer2*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer2*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l2_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);								
											}
				    						}
				    						//Just Enegy on layer 3:
											if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
											if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
												if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer3*tan(0.141421356/2));
												if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer3*tan(0.1/2));
												M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
												MASS_Particle_l3_eMF->Fill(M_particle);
												MASS_Particle_eMF->Fill(M_particle);								
											}
											}





											if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
											if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
											if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
											if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
											if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
											if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
											if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
											if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
											if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
										}
										//---------------------------------------------------------
									}
								}
								else{


		    						//AXION PARTICLE SEARCH-----------------------------
						 	   		//Energy just on layer 1, 2 or 3:
				    				//Just Enegy on layer 1:
									if (Emu1_eMF[k_eMF] != 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																								
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer1*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer1*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l1_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);
									}
									}	
									//Just Enegy on layer 2:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] != 0 && Emu3_eMF[k_eMF] == 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																	
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer2*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer2*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l2_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);								
									}
				    				}
				   			 		//Just Enegy on layer 3:
									if (Emu1_eMF[k_eMF] == 0 && Emu2_eMF[k_eMF] == 0 && Emu3_eMF[k_eMF] != 0){
									if (Ephoton1<Ephotoncut && Ephoton2<Ephotoncut){																																
										if (Eta_photon1-Eta_photon2!=0 && Phi_photon1-Phi_photon2!=0) ang_var=2*atan(Rrextendedlayer3*tan(0.141421356/2));
										if (Eta_photon1-Eta_photon2==0 || Phi_photon1-Phi_photon2==0) ang_var=2*atan(Rrextendedlayer3*tan(0.1/2));
										M_particle=sqrt(2*Ephoton1*Ephoton2*(1-cos(ang_var)));
										MASS_Particle_l3_eMF->Fill(M_particle);
										MASS_Particle_eMF->Fill(M_particle);								
									}
									}




									//AXION PARTICLE SEARCH------------------------------------
									if (ra != 0 || rbc != 0) Emm_nomuon_l12_eMF->Fill(ra,rbc);
									if (ra != 0 || rd != 0) Emm_nomuon_l13_eMF->Fill(ra,rd);
									if (rbc != 0 || rd != 0) Emm_nomuon_l23_eMF->Fill(rbc,rd);
									if (ra>0)  E_nomuon_l1_eMF->Fill(Emu1_eMF[k_eMF]);
									if (rbc>0) E_nomuon_l2_eMF->Fill(Emu2_eMF[k_eMF]);
									if (rd>0)  E_nomuon_l3_eMF->Fill(Emu3_eMF[k_eMF]);
									if (ra>0)  Emm_nomuon_l1_eMF->Fill(ra);
									if (rbc>0) Emm_nomuon_l2_eMF->Fill(rbc);
									if (rd>0)  Emm_nomuon_l3_eMF->Fill(rd);
									//---------------------------------------------------------
								}
								if (ra > 0) Ratio_E_l_eMF->Fill(ra);
								if (rbc > 0) Ratio_E_l_eMF->Fill(rbc);
								if (rd > 0) Ratio_E_l_eMF->Fill(rd);
								if (Emu1_eMF[k_eMF] > 0 && Emu2_eMF[k_eMF] > 0 && Emu3_eMF[k_eMF] > 0) E123_R3_eMF->Fill(Emu1_eMF[k_eMF],Emu2_eMF[k_eMF],Emu3_eMF[k_eMF]);
		  		 			}		
		   				}	
						if (g > 2){	
							Ejets_z1_eMF->Fill(jet_eMF[k_eMF]);
							zii_eMF++;
						}
					}	
				}

				//end muons 

				Zr_eMF[k_eMF] = zr_eMF;

				if (zr_eMF == 1) Ejets_zequal1_eMF->Fill(jet_eMF[k_eMF]);
				if ( zr_eMF > 1 ) {
					Ejets_z1_eMF->Fill(jet_eMF[k_eMF]);
					zi_eMF++;
				}
				if ( zr_eMF > 2) {
					Ejets_z2_eMF->Fill(jet_eMF[k_eMF]);
					zii_eMF++;
				}	
	
				Ejets_eMF->Fill(jet_eMF[k_eMF]);
				k_eMF++;
			}
		}//end an event
	
		for(int jpart=0;jpart<4;jpart++){
	      for(int jmodule=0;jmodule<64;jmodule++){
	        for (int jchannel=0;jchannel<48;jchannel++){
				if (ene[jpart][jmodule][jchannel] > noisecut_ene) cell_ene++;
				if (eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) cell_eMF++;
				if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
					if (ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] > 1.5 || ene[jpart][jmodule][jchannel]/eMF[jpart][jmodule][jchannel][0] < 0.5) {
						if (ene[jpart][jmodule][jchannel] > noisecut_ene) pileupcell_ene++;
						if (eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) pileupcell_eMF++;
					}
				}
			}
		  }
		}

		Dpileupcell_ene = pileupcell_ene;
		Dpileupcell_eMF = pileupcell_eMF;
		Dcell_ene = cell_ene;
		Dcell_eMF = cell_eMF;

		Cell_pileup_eMF_minus_ene->Fill(pileupcell_eMF-pileupcell_ene);
		Cell_pileup_percent_ene->Fill(100*Dpileupcell_ene/12288);
		Cell_pileup_percent_eMF->Fill(100*Dpileupcell_eMF/12288);

		if (pile == 1){
			Nmuon_COFminusOF2_pileup->Fill(muoncandidate_eMF-muoncandidate_ene);
			Njets_COFminusOF2_pileup->Fill(k_eMF-k_ene);

			if (muoncandidate_eMF>muoncandidate_ene) muoncandidate_eMF_bigger++;
			if (muoncandidate_eMF<muoncandidate_ene) muoncandidate_ene_bigger++;
			if (muoncandidate_eMF==muoncandidate_ene) muoncandidate_ene_eMF_equal++;
			if (muoncandidate_eMF != muoncandidate_ene) dif_muon++;
			if (k_eMF != k_ene) dif_jets++;	
			if (zi_ene != zi_eMF) dif_zi++;
			if (zii_ene != zii_eMF) dif_zii++;
		}

		pT_miss_eMF1 = sqrt(px_eMF*px_eMF + py_eMF*py_eMF);

		px_eMF = -1*px_eMF;
		py_eMF = -1*py_eMF;

		if (px_eMF == 0 && py_eMF > 0) phi_miss_eMF1 = PI/2;
		if (px_eMF == 0 && py_eMF < 0) phi_miss_eMF1 = 3*PI/2;
		if (px_eMF > 0 && py_eMF == 0) phi_miss_eMF1 = 0;
		if (px_eMF < 0 && py_eMF == 0) phi_miss_eMF1 = PI;
		if (px_eMF < 0 && py_eMF < 0) phi_miss_eMF1 = atan(py_eMF/px_eMF) + PI;
		if (px_eMF > 0 && py_eMF > 0) phi_miss_eMF1 = atan(py_eMF/px_eMF);
		if (px_eMF > 0 && py_eMF < 0) phi_miss_eMF1 = atan(py_eMF/px_eMF) + 2*PI;
		if (px_eMF < 0 && py_eMF > 0) phi_miss_eMF1 = atan(py_eMF/px_eMF) + PI;

		Phi_miss_eMF->Fill(phi_miss_eMF1);
		PT_miss_eMF->Fill(pT_miss_eMF1);
		eta_miss_eMF->Fill(eta_miss);
		Nmu_per_event_eMF->Fill(zmu_eMF);
		Nmu_per_event_ene_eMF->Fill(zmu_ene,zmu_eMF);
		Nmu_per_event1_eMF->Fill(zmu_eMF1);
		Nmu_per_event2_eMF->Fill(zmu_eMF2);
		Njets_per_event_eMF->Fill(k_eMF);
		Njets_per_event_z1_eMF->Fill(zi_eMF);
		Njets_per_event_zequal1_eMF->Fill(zmu1_eMF);
		Njets_per_event_z2_eMF->Fill(zii_eMF);

		//-----------------------------------------
		//Writing the events
		
		Phi_miss_comparison_ene->Fill(phi_miss_ene2,phi_miss_ene1);
		Phi_miss_comparison_eMF->Fill(phi_miss_eMF2,phi_miss_eMF1);
		PT_miss_comparison_ene->Fill(pT_miss_ene2,pT_miss_ene1);
		PT_miss_comparison_eMF->Fill(pT_miss_eMF2,pT_miss_eMF1);
		Njets_per_event_zequal1_ene_eMF->Fill(zmu1_ene,zmu1_eMF);
		Njets_per_event_ene_eMF->Fill(k_ene,k_eMF);
		Njets_per_event_z1_ene_eMF->Fill(zi_ene,zi_eMF);
		Njets_per_event_z2_ene_eMF->Fill(zii_ene,zii_eMF);

		Eta_Phi_E_ene->Reset();
		Eta_Phi_E_eMF->Reset();
		Eta_Phi_E_ene1->Reset();
		Eta_Phi_E_eMF1->Reset();
		Eta_Phi_E_ene2->Reset();
		Eta_Phi_E_eMF2->Reset();
		Eta_Phi_E_ene3->Reset();
		Eta_Phi_E_eMF3->Reset();

		if (jentry == 100000000000000000){
			cout << "Event = " << jentry << endl;
			for (p = 0; p<= k_ene; p++){
				cout << "OF2 jet nº " << p << endl;
				cout << "zr OF2 = " << zjet_ene[p] << endl;
				cout << "eta OF2 = " << njet_ene[p]*0.1-2. << "; phi OF2 = " << mjet_ene[p]*0.1 << endl;
				cout << "E OF2 = " << jet_ene[p] << endl;
				cout << "E 1st layer OF2 = " << Emu1_ene[p] << endl; 
				cout << "E 2nd layer OF2 = " << Emu2_ene[p] << endl; 
				cout << "E 3rd layer OF2 = " << Emu3_ene[p] << endl; 
				cout << "eta phi part mod channel s0 s1 s2 s3 s4 s5 s6 e_ne" << endl;
				for(int jpart=0;jpart<4;jpart++){
    			  for(int jmodule=0;jmodule<64;jmodule++){
    	   			for (int jchannel=0;jchannel<48;jchannel++){
						if (etacm_ene[p] <= ETA[jpart][jmodule][jchannel]+zjet_ene[p]*0.1 && etacm_ene[p] >= ETA[jpart][jmodule][jchannel]-zjet_ene[p]*0.1 && phicm_ene[p] <= PHI[jpart][jmodule][jchannel]+zjet_ene[p]*0.1 && phicm_ene[p] >= PHI[jpart][jmodule][jchannel]-zjet_ene[p]*0.1 && ene[jpart][jmodule][jchannel]>100) {
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) cout << "pileup channel:" << endl;
							cout << ETA[jpart][jmodule][jchannel] << " " << PHI[jpart][jmodule][jchannel] << " " << jpart << " " << jmodule << " " << jchannel << " " << sample[jpart][jmodule][jchannel][0] << " " << sample[jpart][jmodule][jchannel][1] << " " << sample[jpart][jmodule][jchannel][2] << " " << sample[jpart][jmodule][jchannel][3] << " " << sample[jpart][jmodule][jchannel][4] << " " << sample[jpart][jmodule][jchannel][5] << " " << sample[jpart][jmodule][jchannel][6] << " " << ene[jpart][jmodule][jchannel] << endl;
						}
			 		}
			 	  }
			    }
				cout << "-------------------------------------------------" << endl;
			}

			for (p = 0; p<= k_eMF; p++){
				cout << "COF jet nº " << p << endl;
				cout << "zr COF = " << zjet_eMF[p] << endl;
				cout << "eta COF = " << njet_eMF[p]*0.1-2. << "; phi COF = " << mjet_eMF[p]*0.1 << endl;
				cout << "E COF = " << jet_eMF[p] << endl;
				cout << "E 1st layer COF = " << Emu1_eMF[p] << endl; 
				cout << "E 2nd layer COF = " << Emu2_eMF[p] << endl; 
				cout << "E 3rd layer COF = " << Emu3_eMF[p] << endl;
				cout << "eta phi part mod channel s0 s1 s2 s3 s4 s5 s6 e_MF" << endl;
				for(int jpart=0;jpart<4;jpart++){
    			  for(int jmodule=0;jmodule<64;jmodule++){
    	   			for (int jchannel=0;jchannel<48;jchannel++){
						if (etacm_eMF[p] <= ETA[jpart][jmodule][jchannel]+zjet_eMF[p]*0.1 && etacm_eMF[p] >= ETA[jpart][jmodule][jchannel]-zjet_eMF[p]*0.1 && phicm_eMF[p] <= PHI[jpart][jmodule][jchannel]+zjet_eMF[p]*0.1 && phicm_eMF[p] >= PHI[jpart][jmodule][jchannel]-zjet_eMF[p]*0.1 && eMF[jpart][jmodule][jchannel][0]>100) {
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) cout << "pileup channel: " << endl;
							cout << ETA[jpart][jmodule][jchannel] << " " << PHI[jpart][jmodule][jchannel] << " " << jpart << " " << jmodule << " " << jchannel << " " << sample[jpart][jmodule][jchannel][0] << " " << sample[jpart][jmodule][jchannel][1] << " " << sample[jpart][jmodule][jchannel][2] << " " << sample[jpart][jmodule][jchannel][3] << " " << sample[jpart][jmodule][jchannel][4] << " " << sample[jpart][jmodule][jchannel][5] << " " << sample[jpart][jmodule][jchannel][6] << " " << eMF[jpart][jmodule][jchannel][0] << endl;		
						}
					 	}
				  }
				}
				cout << "-------------------------------------------------" << endl;
			}
			cout << "--------------events-data------------------------" << endl;
			cout << "phi_miss channel by channel method: OF2 = " << phi_miss_ene2 << "; COF = " << phi_miss_eMF2 << endl;
			cout << "pT_miss channel by channel method: OF2 = " << pT_miss_ene2 << "; COF = " << pT_miss_eMF2 << endl;
			cout << "phi_miss CM method: OF2 = " << phi_miss_ene1 << "; COF = " << phi_miss_eMF1 << endl;
			cout << "pT_miss CM method: OF2 = " << pT_miss_ene1 << "; COF = " << pT_miss_eMF1 << endl;
			cout << "-------------------------------------------------" << endl;
			cout << "-------------------------------------------------" << endl;
			//Reading the histogram in the rootfile (n = number of the event):
			//TFile f("Eventn.root")
			//TH1F *Eta_Phi_E = (TH1F*)f.Get("Eta_Phi_E")
		}

		k = k_ene;
		if (k_eMF>k_ene) k = k_eMF;
		for (int jn=0;jn<=k;jn++){		
			if (zjet_ene[jn] == 1 || zjet_eMF[jn] == 1){
				for(int jpart=0;jpart<4;jpart++){
	    		  for(int jmodule=0;jmodule<64;jmodule++){
	       			for (int jchannel=0;jchannel<48;jchannel++){
						if (etacm_eMF[jn] <= ETA[jpart][jmodule][jchannel]+zjet_ene[jn]*0.1 && etacm_eMF[jn] >= ETA[jpart][jmodule][jchannel]-zjet_ene[jn]*0.1 && phicm_eMF[jn] <= PHI[jpart][jmodule][jchannel]+zjet_ene[jn]*0.1 && phicm_eMF[jn] > PHI[jpart][jmodule][jchannel]-zjet_ene[jn]*0.1) {
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
								Ejets_pileup_z1_ene->Fill(jet_ene[jn]);
								Ejets_pileup_z1_eMF->Fill(jet_eMF[jn]);
							}
						}
			 		}
				  }
			   	}
			}
		}

		//Pileup events
	
		if (pile == 1){
			corresponding_jet = 0;
			k = k_ene;
			if (k_eMF>k_ene){
				k_eMF_bigger++;
				k = k_eMF;	
			}
			if (k_ene>k_eMF) k_ene_bigger++;
			if (k_ene==k_eMF) k_ene_eMF_equal++;
			jetmatch=0;
			for (int jk_eMF=0; jk_eMF<=k_eMF; jk_eMF++){
				corresponding_jet = 0;
				for (int jk_ene=0; jk_ene<=k_ene; jk_ene++){
					//Check if the COF jet match with some OF2 jet
					if (etacm_ene[jk_ene] >= etacm_eMF[jk_eMF]-0.25 && etacm_ene[jk_ene] <= etacm_eMF[jk_eMF]+0.25 && phicm_ene[jk_ene] >= phicm_eMF[jk_eMF]-0.25 && phicm_ene[jk_ene] <= phicm_eMF[jk_eMF]+0.25){
						corresponding_jet = 1; //Yes, check
						zr=Zr_eMF[jk_eMF];
						if (Zr_ene[jk_ene]>Zr_eMF[jk_eMF]) {
							zr = Zr_ene[jk_ene];
						}
						//Check if this jet has some pileup cell
						jetpile=0;
						jetmatch=0;
						for(int jpart=0;jpart<4;jpart++){
	    				  for(int jmodule=0;jmodule<64;jmodule++){
	            			for (int jchannel=0;jchannel<48;jchannel++){
								//if (jetmatch==0){
									//Condition to be pileup
									if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
										//Condition to be part of the jet
										if ((etacm_eMF[jk_eMF]+etacm_ene[jk_ene])/2 <= ETA[jpart][jmodule][jchannel]+(zr+1)*0.1 && (etacm_eMF[jk_eMF]+etacm_ene[jk_ene])/2 >= ETA[jpart][jmodule][jchannel]-(zr+1)*0.1 && (phicm_eMF[jk_eMF]+phicm_ene[jk_ene])/2 <= PHI[jpart][jmodule][jchannel]+(zr+1)*0.1 && (phicm_eMF[jk_eMF]+phicm_ene[jk_ene])/2 >= PHI[jpart][jmodule][jchannel]-(zr+1)*0.1) {
											jetpile=1;
											jetmatch=1;
											//Ratio_jet_pile_energy->Fill(jet_eMF[jk_eMF]/jet_ene[jk_ene]);
											Ejet_pileup_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
											if (jet_eMF[jk_eMF] != 0 && jet_ene[jk_ene] != 0){
												Zr_jetdif_pileup->Fill(Zr_eMF[jk_eMF]-Zr_ene[jk_ene]);

												/*
												//Check what kind of jet it is and plot the Energy Histogram
												if (Zr_ene[jk_ene] > 1 || Zr_eMF[jk_eMF] > 1){
	   			 								    Energy_jet_pile_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
													Energy_jet_pile_ene->Fill(jet_ene[jk_ene]);
													Energy_jet_pile_eMF->Fill(jet_eMF[jk_eMF]);
													//jetpile = 1;
												}
												if (Zr_ene[jk_ene] == 1 || Zr_eMF[jk_eMF] == 1){
													Energy_pile_z1_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
													Energy_pile_z1_ene->Fill(jet_ene[jk_ene]);
													Energy_pile_z1_eMF->Fill(jet_eMF[jk_eMF]);
													//z1pile = 1;
												}
												if (ismucand_ene[jk_ene] == 1 || ismucand_eMF[jk_eMF] == 1){
	    										    Energy_mucand_pile_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
	   											    Energy_mucand_pile_ene->Fill(jet_ene[jk_ene]);
											        Energy_mucand_pile_eMF->Fill(jet_eMF[jk_eMF]);
													//muonpile = 1;
												}
												*/
											}
	
											//Energy histogram of all jets, with no distinction 
											Ejets_pileup_eMF->Fill(jet_eMF[k_eMF]);
											Ejets_pileup_ene->Fill(jet_ene[k_ene]);
											Ejets_pileup_ene_eMF->Fill(jet_ene[k_ene],jet_eMF[k_eMF]);
											//Calculate the Number_of_pileup_cells in each COF and OF2 jet
											if (Zr_ene[jk_ene] == 1 && ene[jpart][jmodule][jchannel] > noisecut_ene) pile_z1_ene++; 
											if (Zr_eMF[jk_eMF] == 1 && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) pile_z1_eMF++; 
											if (Zr_ene[jk_ene] > 1 && ene[jpart][jmodule][jchannel] > noisecut_ene) jet_pile_ene++; 
											if (Zr_eMF[jk_eMF] > 1 && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) jet_pile_eMF++; 
											if (ismucand_ene[jk_ene] == 1 && ene[jpart][jmodule][jchannel] > noisecut_ene) mucand_pile_ene++; 
											if (ismucand_eMF[jk_eMF] == 1 && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) mucand_pile_eMF++;
											if (ene[jpart][jmodule][jchannel] > noisecut_ene) all_pile_ene++;
											if (eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) all_pile_eMF++;	
										}
									}
									if (jetpile==1 && corresponding_jet == 1){
										Ratio_jet_pile_energy->Fill(jet_eMF[jk_eMF]/jet_ene[jk_ene]);
										Ejet_pileup_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
									}
									if (jetpile==0 && corresponding_jet == 1){
										Ratio_jet_nopile_energy->Fill(jet_eMF[jk_eMF]/jet_ene[jk_ene]);
										Ejet_nopileup_ene_eMF->Fill(jet_ene[jk_ene],jet_eMF[jk_eMF]);
									}	
								//}
							}
	 					  } 
					  	}
					}
				}//stop trying to match a OF2 jet
				if (corresponding_jet == 1){
					Number_jet_pile_eMF->Fill(jet_pile_eMF);
					Number_jet_pile_ene->Fill(jet_pile_ene);
					Number_jet_pile_ene_eMF->Fill(jet_pile_ene,jet_pile_eMF);
					Number_pile_z1_ene->Fill(pile_z1_ene);
					Number_pile_z1_eMF->Fill(pile_z1_eMF);
					Number_pile_z1_ene_eMF->Fill(pile_z1_ene,pile_z1_eMF);
					Number_mucand_pile_ene->Fill(mucand_pile_ene);
					Number_mucand_pile_eMF->Fill(mucand_pile_eMF);
					Number_mucand_pile_ene_eMF->Fill(mucand_pile_ene,mucand_pile_eMF);
				}

				//Check which COF jet don't match with any OF2 jet and check if it has pileup cells 
				//-0.5 -> COF don't match with OF2
				//-1.0 -> OF2 don't match with COF
				if (corresponding_jet != 1) {
					//if (jet_eMF[jk_eMF] > 1000)cout << "COF reconstruct and OF2 don't = " << jentry << endl;
					//if (jetpile==1) Ratio_jet_pile_energy->Fill(-0.5);
					//if (jetpile==0) Ratio_jet_nopile_energy->Fill(-0.5);
					if (Zr_eMF[jk_eMF] == 1) zr1_dif++;
					zr_dif++;

					//Check if this jet has some pileup cell
					for(int jpart=0;jpart<4;jpart++){
		   		 	  for(int jmodule=0;jmodule<64;jmodule++){
	        			for (int jchannel=0;jchannel<48;jchannel++){
							//Condition to be pileup
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
								//Condition to be part of the jet
								if (etacm_eMF[jk_eMF] <= ETA[jpart][jmodule][jchannel]+(Zr_eMF[jk_eMF]+1)*0.1 && etacm_eMF[jk_eMF] >= ETA[jpart][jmodule][jchannel]-(Zr_eMF[jk_eMF]+1)*0.1 && phicm_eMF[jk_eMF] <= PHI[jpart][jmodule][jchannel]+(Zr_eMF[jk_eMF]+1)*0.1 && phicm_eMF[jk_eMF] >= PHI[jpart][jmodule][jchannel]-(Zr_eMF[jk_eMF]+1)*0.1) {

									if (jet_eMF[jk_eMF] != 0 && eMF[jpart][jmodule][jchannel][0]> noisecut_eMF){
										if (Zr_eMF[jk_eMF] > 1) {
											Energy_alone_jet_pile_eMF->Fill(jet_eMF[jk_eMF]);
											jet_pile_al_eMF++;
										}			
										if (Zr_eMF[jk_eMF] == 1){
											pile_z1_al_eMF++;
											Energy_alone_pile_z1_eMF->Fill(jet_eMF[jk_eMF]);
										}	
										if (ismucand_eMF[jk_eMF] == 1){
											mucand_pile_al_eMF++;
											Energy_alone_mucand_pile_eMF->Fill(jet_eMF[jk_eMF]);
										}
										Ejets_alone_pileup_eMF->Fill(jet_eMF[jk_eMF]);
									}
								}
				  			}
				 	   }
					  }
					}
					Number_alone_jet_pile_eMF->Fill(jet_pile_al_eMF);
					Number_alone_pile_z1_eMF->Fill(pile_z1_al_eMF);
					Number_alone_mucand_pile_eMF->Fill(mucand_pile_al_eMF);
				}
			}//finish COF jets count

	
			//Check which OF2 jet don't match with any COF jet and check if it has pileup cells 
			//-0.5 -> COF don't match with OF2
			//-1.0 -> OF2 don't match with COF
			for (int jk_ene=0; jk_ene<=k_ene; jk_ene++){
				corresponding_jet = 0;
				for (int jk_eMF=0; jk_eMF<=k_eMF; jk_eMF++){
					if (etacm_eMF[jk_eMF] >= etacm_ene[jk_ene]-0.25 && etacm_eMF[jk_eMF] <= etacm_ene[jk_ene]+0.25 && phicm_eMF[jk_eMF] >= phicm_ene[jk_ene]-0.25 && phicm_eMF[jk_eMF] <= phicm_ene[jk_ene]+0.25) {
						corresponding_jet = 1;
					}
				}
				//That's the condition to don't match
				if (corresponding_jet != 1) {
					//if (jet_ene[jk_ene] > 1000) cout << "OF2 reconstruct and COF don't = " << jentry << endl;
					//if (jetpile==1) Ratio_jet_pile_energy->Fill(-1.0);
					//if (jetpile==0) Ratio_jet_nopile_energy->Fill(-1.0);
					if (Zr_ene[jk_ene] == 1) zr1_dif++;
					Energy_alone_jet_pileup_ene->Fill(jet_ene[jk_ene]);
					zr_dif++;

					//Check if this jet has some pileup cell
					for(int jpart=0;jpart<4;jpart++){
			   	  	  for(int jmodule=0;jmodule<64;jmodule++){
	    	       	    for (int jchannel=0;jchannel<48;jchannel++){
							//Condition to be pileup
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
								//Condition to be part of the jet
								if (etacm_ene[jk_ene] <= ETA[jpart][jmodule][jchannel]+(Zr_ene[jk_ene]+1)*0.1 && etacm_ene[jk_ene] >= ETA[jpart][jmodule][jchannel]-(Zr_ene[jk_ene]+1)*0.1 && phicm_ene[jk_ene] <= PHI[jpart][jmodule][jchannel]+(Zr_ene[jk_ene]+1)*0.1 && phicm_ene[jk_ene] >= PHI[jpart][jmodule][jchannel]-(Zr_ene[jk_ene]+1)*0.1) {

									if (jet_ene[jk_ene] != 0 && ene[jpart][jmodule][jchannel]> noisecut_ene){
										if (Zr_eMF[jk_ene] > 1) {
											Energy_alone_jet_pile_ene->Fill(jet_ene[jk_ene]);
											jet_pile_al_ene++;
										}
										if (Zr_eMF[jk_ene] == 1){
											Energy_alone_pile_z1_ene->Fill(jet_ene[jk_ene]);
											pile_z1_al_ene++;
										}
										if (ismucand_ene[jk_ene] == 1){
											Energy_alone_mucand_pile_ene->Fill(jet_ene[jk_ene]);
											mucand_pile_al_ene++;
										}
										Ejets_alone_pileup_ene->Fill(jet_ene[jk_ene]);
									}
								}
							}
				    	}
					  }
					}

					Number_alone_jet_pile_ene->Fill(jet_pile_al_ene);
					Number_alone_pile_z1_ene->Fill(pile_z1_al_ene);
					Number_alone_mucand_pile_ene->Fill(mucand_pile_al_ene);
				}
			}

			Doubjentry = jentry;

			for(int jpart=0;jpart<4;jpart++){
    	  	  for(int jmodule=0;jmodule<64;jmodule++){
    	   	    for (int jchannel=0;jchannel<48;jchannel++){
					//Condition to be pileup
					if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
						if (ene[jpart][jmodule][jchannel] > noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) {	
						 	Energy_alone_cell_pile_ene_eMF->Fill(ene[jpart][jmodule][jchannel],eMF[jpart][jmodule][jchannel][0]);
							Eta_Phi_E_pileup_eMF->Fill(ETA[jpart][jmodule][jchannel],PHI[jpart][jmodule][jchannel],eMF[jpart][jmodule][jchannel][0]);
							Eta_Phi_E_pileup_ene->Fill(ETA[jpart][jmodule][jchannel],PHI[jpart][jmodule][jchannel],ene[jpart][jmodule][jchannel]);
						}
					}
				}
			  }
			}

			for (int n=0;n<=40;n++){
				for (int m=0;m<=64;m++){
					Energy_alone_tower_pile_ene_eMF->Fill(Eta_Phi_E_pileup_ene->GetBinContent(n+1,m+1),Eta_Phi_E_pileup_eMF->GetBinContent(n+1,m+1));
				}
			}

			//Checking the OF2 muon with COF with E_cell>0
			for (int jk_ene=0; jk_ene<=k_ene; jk_ene++){
				if (ismucand_ene[jk_ene]==1){
					Energy_alone_muon_pile_ene_eMF->Fill(jet_ene[jk_ene],Eta_Phi_E_pileup_eMF->GetBinContent(njet_ene[k_ene]+1,mjet_ene[k_ene]+1));
				}
			}

			for(int part=0; part<4;part++){
   		      for(int module=0; module<64; module++){
    	    	for(int channel=0; channel<48; channel++){
					if (tMF[part][module][channel][0] > -25 && tMF[part][module][channel][0] < 25 && pedMF[part][module][channel] < 100 && chi2MF[part][module][channel] > 0.01){
						int scinti = 0;
						if (channel == 0) scinti = 1;
						if (channel == 1) scinti = 1;
						if (channel == 12) scinti = 1;
						if (channel == 13) scinti = 1;
						if (scinti == 0) {
    	        		    if (ene[part][module][channel] > noisecut_ene && chi2[part][module][channel] < QFcut) {
								Eta_Phi_E2_ene->Fill(ETA[part][module][channel],PHI[part][module][channel],ene[part][module][channel]);	
							}
    	            		if (eMF[part][module][channel][0] > noisecut_eMF && chi2MF[part][module][channel] < QFcut) {
								Eta_Phi_E2_eMF->Fill(ETA[part][module][channel],PHI[part][module][channel],eMF[part][module][channel][0]);
							}	
						}	
					}
				}
			  }	
			}

			for(int jpart=0;jpart<4;jpart++){
	      	  for(int jmodule=0;jmodule<64;jmodule++){
	       	    for (int jchannel=0;jchannel<48;jchannel++){
					if (tMF[jpart][jmodule][jchannel][0] > -25 && tMF[jpart][jmodule][jchannel][0] < 25 && pedMF[jpart][jmodule][jchannel] < 100 && chi2MF[jpart][jmodule][jchannel] > 0.01){
						int scinti2 = 0;
						if (jchannel == 0) scinti2 = 1;
						if (jchannel == 1) scinti2 = 1;
						if (jchannel == 12) scinti2 = 1;
						if (jchannel == 13) scinti2 = 1;
						if (scinti2 == 0) {
							//Condition to be pileup
							n = 10*(2+ETA[jpart][jmodule][jchannel]);
							m = 10*(PHI[jpart][jmodule][jchannel]);
							cellpile=0;
							if (sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][4] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][2] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][3]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][1] || sample[jpart][jmodule][jchannel][2]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][5] || sample[jpart][jmodule][jchannel][4]<sample[jpart][jmodule][jchannel][6] || sample[jpart][jmodule][jchannel][1]<sample[jpart][jmodule][jchannel][0] || sample[jpart][jmodule][jchannel][5]<sample[jpart][jmodule][jchannel][6]) {
								if (ene[jpart][jmodule][jchannel]!=0){
									all_Edif0_cell_ene++;
									all_Edif0_cell_eMF++;
									if (ene[jpart][jmodule][jchannel]>=noisecut_ene && eMF[jpart][jmodule][jchannel][0]<=noisecut_ene) all_rec_cell_ene++;
									if (ene[jpart][jmodule][jchannel]<=noisecut_ene && eMF[jpart][jmodule][jchannel][0]>=noisecut_ene) all_rec_cell_eMF++;
								}
								if (ene[jpart][jmodule][jchannel]>=noisecut_ene && eMF[jpart][jmodule][jchannel][0] <= noisecut_eMF)     Ecell_pileup_OF2alone_ene_eMF->Fill(ene[jpart][jmodule][jchannel],eMF[jpart][jmodule][jchannel][0]);
								if (ene[jpart][jmodule][jchannel]<=noisecut_ene && eMF[jpart][jmodule][jchannel][0] >= noisecut_eMF)     Ecell_pileup_COFalone_eMF_ene->Fill(eMF[jpart][jmodule][jchannel][0],ene[jpart][jmodule][jchannel]);
								cellpile=1;
								//E cell with pileup
								if (ene[jpart][jmodule][jchannel] > noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF){
									Ecell_pileup_ene_eMF->Fill(ene[jpart][jmodule][jchannel],eMF[jpart][jmodule][jchannel][0]);
									Ratio_cell_pile_energy->Fill(eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]);
									Ratio_cell_all++;
									if (eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]>0 && eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]<0.5) Ratio_cell_0to05++; 
									if (eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]>0.5 && eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]<1.5) Ratio_cell_05to15++; 
									if (eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]>1.5 && eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]<3.0) Ratio_cell_15to3++; 
									if (eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]>3.0 && eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]<10.0) Ratio_cell_3to10++; 
								}

								if (Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>noisecut_ene && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)>noisecut_eMF){
									Etower_pileup_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));
									Ratio_tower_pile_energy->Fill(Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1));
									Ratio_tower_all++;
									if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>0 && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)<0.5) Ratio_tower_0to05++; 
									if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>0.5 && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)<1.5) Ratio_tower_05to15++;
									if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>1.5 && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)<3.0) Ratio_tower_15to3++; 
									if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>3.0 && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1)<10.0) Ratio_tower_3to10++; 
								}

								if (ene[jpart][jmodule][jchannel] < noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) {
									Ratio_cell_pile_energy->Fill(-0.5);
									Ratio_cell_all++;
									Ratio_cell_justCOF++;
								}
	
								if (eMF[jpart][jmodule][jchannel][0] < noisecut_eMF && ene[jpart][jmodule][jchannel] > noisecut_ene) {
									Ratio_cell_pile_energy->Fill(-1.0);
									Ratio_cell_all++;
									Ratio_cell_justOF2++;
								}

								if (Eta_Phi_E2_ene->GetBinContent(n+1,m+1) < noisecut_ene && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1) > noisecut_eMF) {
									Ratio_tower_pile_energy->Fill(-0.5);
									Ratio_tower_all++;
									Ratio_tower_justCOF++;
								}

								if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1) < noisecut_eMF && Eta_Phi_E2_ene->GetBinContent(n+1,m+1) > noisecut_ene) {
									Ratio_tower_pile_energy->Fill(-1.0);
									Ratio_tower_all++;
									Ratio_tower_justOF2++;
								}

								if (Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>0 || Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)>0){

									if (ene[jpart][jmodule][jchannel] < 0 && eMF[jpart][jmodule][jchannel][0] > 0) only_COF++;
									if (ene[jpart][jmodule][jchannel] > 0 && eMF[jpart][jmodule][jchannel][0] < 0) only_OF2++;

									if (ene[jpart][jmodule][jchannel] < noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF && chi2MF[jpart][jmodule][jchannel] < QFcut) {	

										Energy_alone_tower2_pile_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));
										Energy_alone_tower2_pileCOF_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));

										if (ene[jpart][jmodule][jchannel]<0) minus0_ene++;
										if (ene[jpart][jmodule][jchannel]>0 && ene[jpart][jmodule][jchannel]<100) I0to100_ene++;
										if (ene[jpart][jmodule][jchannel]>100 && ene[jpart][jmodule][jchannel]<200) I100to200_ene++;
										if (ene[jpart][jmodule][jchannel]>200 && ene[jpart][jmodule][jchannel]<300) I200to300_ene++;

										if (eMF[jpart][jmodule][jchannel][0]==Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)) just1cell_eMF++;
										if (eMF[jpart][jmodule][jchannel][0]!=Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)) more1cell_eMF++;
									}

									if (ene[jpart][jmodule][jchannel] > noisecut_ene && eMF[jpart][jmodule][jchannel][0] < noisecut_eMF && Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>0 && chi2[jpart][jmodule][jchannel] < QFcut) {
										n = 10*(2+ETA[jpart][jmodule][jchannel]);
										m = 10*(PHI[jpart][jmodule][jchannel]);	
										if (ene[jpart][jmodule][jchannel]==Eta_Phi_E2_ene->GetBinContent(n+1,m+1)) just1cell_ene++;
										if (ene[jpart][jmodule][jchannel]!=Eta_Phi_E2_ene->GetBinContent(n+1,m+1)) more1cell_ene++;
										if (eMF[jpart][jmodule][jchannel][0]<0) minus0_eMF++;
										if (eMF[jpart][jmodule][jchannel][0]>0 && eMF[jpart][jmodule][jchannel][0]<100) I0to100_eMF++;
										if (eMF[jpart][jmodule][jchannel][0]>100 && eMF[jpart][jmodule][jchannel][0]<200) I100to200_eMF++;
										if (eMF[jpart][jmodule][jchannel][0]>200 && eMF[jpart][jmodule][jchannel][0]<300) I200to300_eMF++;
										Energy_alone_tower2_pileOF2_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));
										Energy_alone_tower2_pile_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));
									}
								}
							}

							if (cellpile==0) {

								if (ene[jpart][jmodule][jchannel] > noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF){
									Ecell_nopileup_ene_eMF->Fill(ene[jpart][jmodule][jchannel],eMF[jpart][jmodule][jchannel][0]);
									Ratio_cell_nopile_energy->Fill(eMF[jpart][jmodule][jchannel][0]/ene[jpart][jmodule][jchannel]);
								}

								if (Eta_Phi_E2_ene->GetBinContent(n+1,m+1)>noisecut_ene && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)>noisecut_eMF){
									Etower_nopileup_ene_eMF->Fill(Eta_Phi_E2_ene->GetBinContent(n+1,m+1),Eta_Phi_E2_eMF->GetBinContent(n+1,m+1));
									Ratio_tower_nopile_energy->Fill(Eta_Phi_E2_eMF->GetBinContent(n+1,m+1)/Eta_Phi_E2_ene->GetBinContent(n+1,m+1));
								}

								if (ene[jpart][jmodule][jchannel] < noisecut_ene && eMF[jpart][jmodule][jchannel][0] > noisecut_eMF) Ratio_cell_nopile_energy->Fill(-0.5);
								if (eMF[jpart][jmodule][jchannel][0] < noisecut_eMF && ene[jpart][jmodule][jchannel] > noisecut_ene) Ratio_cell_nopile_energy->Fill(-1.0);
								if (Eta_Phi_E2_ene->GetBinContent(n+1,m+1) < noisecut_ene && Eta_Phi_E2_eMF->GetBinContent(n+1,m+1) > noisecut_eMF)  Ratio_tower_nopile_energy->Fill(-0.5);
								if (Eta_Phi_E2_eMF->GetBinContent(n+1,m+1) < noisecut_eMF && Eta_Phi_E2_ene->GetBinContent(n+1,m+1) > noisecut_ene)  Ratio_tower_nopile_energy->Fill(-1.0);
							}
						}
					}
				}
			  }
			}


			Eta_Phi_E_pileup_ene->Reset();
			Eta_Phi_E_pileup_eMF->Reset();
			Eta_Phi_E2_ene->Reset();
			Eta_Phi_E2_eMF->Reset();
		}
		//End made by Joao Pedro!
	}//end all events

	
    //Create the legend
    leg2 = new TLegend(0.7,0.4,0.8,0.65);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg2->SetBorderSize(0);
    leg2->AddEntry(Ejets_eMF,"COF","l");
    leg2->AddEntry(Ejets_ene,"OF2","l");
    leg2->SetTextSize(0.035);

    leg7 = new TLegend(0.7,0.4,0.8,0.55);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg7->SetBorderSize(0);
    leg7->AddEntry(Ejets_eMF,"COF","l");
    leg7->AddEntry(Ejets_ene,"OF2","l");
    leg7->SetTextSize(0.035);

    leg8 = new TLegend(0.3,0.6,0.45,0.8);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg8->SetBorderSize(0); 
    leg8->AddEntry(Ejets_eMF,"COF","l");
    leg8->AddEntry(Ejets_ene,"OF2","l");
    leg8->SetTextSize(0.035);

    leg9 = new TLegend(0.4,0.75,0.6,0.85);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg9->SetBorderSize(0);
    leg9->AddEntry(Ejets_eMF,"COF","l");
    leg9->AddEntry(Ejets_ene,"OF2","l");
    leg9->SetTextSize(0.035);

    leg4 = new TLegend(0.3,0.7,0.28,0.8);
    //leg = new TLegend(0.10.75,0.3,0.9);
    leg4->SetBorderSize(0);
    leg4->AddEntry(White,"COF noise cut = "+Snoisecut_eMF+"MeV","l");
    leg4->AddEntry(White,"OF2 noise cut = 300MeV","l");
    leg4->SetTextSize(0.035);

    leg5 = new TLegend(0.5,0.7,0.8,0.8);
    //leg = new TLegend(0.10.75,0.3,0.9);
    leg5->SetBorderSize(0);
    leg5->AddEntry(White,"-0.5 = COF alone","l");
    leg5->AddEntry(White,"-1.0 = OF2 alone","l");
    leg5->SetTextSize(0.035);

	

	/*
	Double_t Dzr1_dif = zr1_dif;
	Double_t Dzr_dif = zr_dif;
	Double_t Donly_COF = only_COF;
	Double_t Donly_OF2 = only_OF2;
	Double_t Djust1cell_ene = just1cell_ene;
	Double_t Dmore1cell_ene = more1cell_ene;
	Double_t Djust1cell_eMF = just1cell_eMF;
	Double_t Dmore1cell_eMF = more1cell_eMF;
	Double_t Dminus0_eMF = minus0_eMF;
	Double_t D0to100_eMF = I0to100_eMF;
	Double_t D100to200_eMF = I100to200_eMF;
	Double_t D200to300_eMF = I200to300_eMF;
	Double_t Dminus0_ene = minus0_ene;
	Double_t D0to100_ene = I0to100_ene;
	Double_t D100to200_ene = I100to200_ene;
	Double_t D200to300_ene = I200to300_ene;

	cout << "Onde COF reconstrói e OF2 não, ou vice-versa, quantos são Zr=1: " << 100*Dzr1_dif/Dzr_dif << "%" << endl;
	cout << "Where COF reconstruct but OF2 don't " << 100*Donly_COF/(Donly_COF+Donly_OF2) << "%" << endl;
	cout << "Where OF2 reconstruct but COF don't " << 100*Donly_OF2/(Donly_COF+Donly_OF2) << "%" << endl;

	cout << "OF2 alone - tower with just one cell = " << 100*Djust1cell_ene/(Djust1cell_ene+Dmore1cell_ene) << "%" << endl;
	cout << "OF2 alone - tower more than one cell = " << 100*Dmore1cell_ene/(Djust1cell_ene+Dmore1cell_ene) << "%" << endl;
	cout << "COF alone - tower with just one cell = " << 100*Djust1cell_eMF/(Djust1cell_eMF+Dmore1cell_eMF) << "%" << endl;
	cout << "COF alone - tower more than one cell = " << 100*Dmore1cell_eMF/(Djust1cell_eMF+Dmore1cell_eMF) << "%" << endl;

	cout << "When only OF2 reconstruct the tower: " << endl;
	cout <<"E_COF < 0 = " << 100*Dminus0_eMF/(Dminus0_eMF+D0to100_eMF+D100to200_eMF+D200to300_eMF) << "%" << endl;
	cout <<"0 < E_COF < 100 = " << 100*D0to100_eMF/(Dminus0_eMF+D0to100_eMF+D100to200_eMF+D200to300_eMF) << "%" << endl;
	cout <<"100 < E_COF < 200 = " << 100*D100to200_eMF/(Dminus0_eMF+D0to100_eMF+D100to200_eMF+D200to300_eMF) << "%" << endl;
	cout <<"200 < E_COF < 300 = " << 100*D200to300_eMF/(Dminus0_eMF+D0to100_eMF+D100to200_eMF+D200to300_eMF) << "%" << endl;

	cout << "When only COF reconstruct the tower: " << endl;
	cout <<"E_OF2 < 0 = " << 100*Dminus0_ene/(Dminus0_ene+D0to100_ene+D100to200_ene+D200to300_ene) << "%" << endl;
	cout <<"0 < E_OF2 < 100 = " << 100*D0to100_ene/(Dminus0_ene+D0to100_ene+D100to200_ene+D200to300_ene) << "%" << endl;
	cout <<"100 < E_OF2 < 200 = " << 100*D100to200_ene/(Dminus0_ene+D0to100_ene+D100to200_ene+D200to300_ene) << "%" << endl;
	cout <<"200 < E_OF2 < 300 = " << 100*D200to300_ene/(Dminus0_ene+D0to100_ene+D100to200_ene+D200to300_ene) << "%" << endl;

	pol1_1 = new TF1("pol1_1","pol1",0,10);
	pol1_2 = new TF1("pol1_2","pol1",0,50);
	pol1_3 = new TF1("pol1_3","pol1",0,100);

	TCanvas *pl237 = new TCanvas("pl237","pl237",900,700);
  	Energy_alone_cell_pile_ene_eMF->Draw();
  	//Energy_alone_cell_pile_ene_eMF->Fit(pol1_1);
    pl237->SaveAs("E_cell_alone_pile_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
    pl237->Close();

	TCanvas *pl231 = new TCanvas("pl231","pl231",900,700);
  	Energy_alone_tower_pile_ene_eMF->Draw();
  	//Energy_alone_tower_pile_ene_eMF->Fit(pol1_2);
    pl231->SaveAs("E_tower_alone_pile_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
    pl231->Close();

	TCanvas *pl323 = new TCanvas("pl323","pl323",900,700);
  	Energy_alone_muon_pile_ene_eMF->Draw();
  	//Energy_alone_cell_pile_ene_eMF->Fit(pol1_1);
    pl323->SaveAs("E_muon_alone_pile_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
	pl323->Close();

	TCanvas *pl2316 = new TCanvas("pl2316","pl2316",900,700);
  	Energy_alone_tower2_pile_ene_eMF->Draw();
  	//Energy_alone_tower_pile_ene_eMF->Fit(pol1_2);
    pl2316->SaveAs("E_tower_with_cell_al_pile_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
    pl2316->Close();

	TCanvas *pl2317 = new TCanvas("pl2317","pl2317",900,700);
  	Energy_alone_tower2_pileCOF_ene_eMF->Draw();
  	//Energy_alone_tower_pile_ene_eMF->Fit(pol1_2);
    pl2317->SaveAs("E_tower_with_cell_al_pileCOF_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
    pl2317->Close();

	TCanvas *pl2318 = new TCanvas("pl2318","pl2318",900,700);
  	Energy_alone_tower2_pileOF2_ene_eMF->Draw();
  	//Energy_alone_tower_pile_ene_eMF->Fit(pol1_2);
    pl2318->SaveAs("E_tower_with_cell_al_pileOF2_noisecutCOF"+Snoisecut_eMF+"_eneversuseMFQF" + SQFcut + ".pdf","pdf");
    pl2318->Close();

	*/
	/*

	cout << "" << endl;
	cout << "Jets Zr=1:" << endl;
	TCanvas *pl121 = new TCanvas("pl121","pl21",900,700);
  	Number_pile_z1_ene_eMF->Draw();
	//Number_pile_z1_ene_eMF->Fit(pol1_1);
    pl121->SaveAs("N_pile_z1_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl121->Close();

	cout << "" << endl;
	cout << "Muon candidate: " << endl;
	TCanvas *pl122 = new TCanvas("pl122","pl122",900,700);
  	Number_mucand_pile_ene_eMF->Draw();
  	//Number_mucand_pile_ene_eMF->Fit(pol1_2);
    pl122->SaveAs("N_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl122->Close();

	cout << "" << endl;
	cout << "Jets Zr>1: " << endl;
	TCanvas *pl123 = new TCanvas("pl123","pl123",900,700);
  	Number_jet_pile_ene_eMF->Draw();
  	//Number_jet_pile_ene_eMF->Fit(pol1_3);
    pl123->SaveAs("N_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl123->Close();

	TCanvas *pl124 = new TCanvas("pl124","pl24",900,700);
  	Number_pile_z1_eMF->SetTitle("Number of cells with pileup in each jet (Zr=1)(QF<" + SQFcut + ")");
  	Number_pile_z1_eMF->Draw();
  	Number_pile_z1_ene->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
	pl124->SetLogy();
    pl124->SaveAs("N_pile_z1_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl124->Close();

	TCanvas *pl725 = new TCanvas("pl725","pl725",900,700);
  	Number_mucand_pile_ene->SetTitle("Number of cells with pileup in each muon candidate (QF<" + SQFcut + ")");
  	Number_mucand_pile_ene->Draw();
  	Number_mucand_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
	pl725->SetLogy();
    pl725->SaveAs("N_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl725->Close();

	TCanvas *pl726 = new TCanvas("pl726","pl726",900,700);
  	Number_jet_pile_ene->SetTitle("Number of cells with pileup in each jet (Zr>1) (QF<" + SQFcut + ")");
  	Number_jet_pile_ene->Draw();
  	Number_jet_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
	pl726->SetLogy();
    pl726->SaveAs("N_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl726->Close();

	TCanvas *pl224 = new TCanvas("pl224","pl224",900,700);
  	Number_alone_pile_z1_eMF->Draw();
  	leg4->Draw("SAME");
	pl224->SetLogy();
    pl224->SaveAs("N_pile_z1_al_eMF_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl224->Close();

	TCanvas *pl234 = new TCanvas("pl234","pl234",900,700);
  	Number_alone_pile_z1_ene->Draw();
  	leg4->Draw("SAME");
	pl234->SetLogy();
    pl234->SaveAs("N_pile_z1_al_ene_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl234->Close();

	TCanvas *pl324 = new TCanvas("pl324","pl324",900,700);
  	Number_alone_jet_pile_eMF->Draw();
  	leg4->Draw("SAME");
	pl324->SetLogy();
    pl324->SaveAs("N_jet_pile_al_eMF_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl324->Close();

	TCanvas *pl434 = new TCanvas("pl434","pl434",900,700);
  	Number_alone_jet_pile_ene->Draw();
  	leg4->Draw("SAME");
	pl434->SetLogy();
    pl434->SaveAs("N_jet_pile_al_ene_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl434->Close();

	TCanvas *pl354 = new TCanvas("pl354","pl354",900,700);
  	Number_alone_mucand_pile_eMF->Draw();
  	leg4->Draw("SAME");
	pl354->SetLogy();
    pl354->SaveAs("N_mucand_pile_al_eMF_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl354->Close();

	TCanvas *pl464 = new TCanvas("pl464","pl464",900,700);
  	Number_alone_mucand_pile_ene->Draw();
  	leg4->Draw("SAME");
	pl464->SetLogy();
    pl464->SaveAs("N_mucand_pile_al_ene_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl464->Close();

	TCanvas *pl125 = new TCanvas("pl125","pl125",900,700);
  	Number_mucand_pile_ene->SetTitle("Number ofof cells with pileup in each muon candidate (QF<" + SQFcut + ")");
  	Number_mucand_pile_ene->Draw();
  	Number_mucand_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
	pl125->SetLogy();
    pl125->SaveAs("N_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl125->Close();

	TCanvas *pl126 = new TCanvas("pl126","pl126",900,700);
  	Number_jet_pile_ene->SetTitle("Number of cells with pileup in each jet (Zr>1) (QF<" + SQFcut + ")");
  	Number_jet_pile_ene->Draw();
  	Number_jet_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
	pl126->SetLogy();
    pl126->SaveAs("N_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl126->Close();


	pol1_1 = new TF1("pol1_1","pol1",300,10000);
	pol1_2 = new TF1("pol1_2","pol1",300,5000);
	pol1_3 = new TF1("pol1_3","pol1",300,100000);

	cout << "" << endl;
	cout << "Jets Zr=1:" << endl;
	TCanvas *pl21 = new TCanvas("pl21","pl21",900,700);
  	Energy_pile_z1_ene_eMF->Draw();
	//Energy_pile_z1_ene_eMF->Fit(pol1_1);
    pl21->SaveAs("E_pile_z1_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl21->Close();

	cout << "" << endl;
	cout << "Muon candidate: " << endl;
	TCanvas *pl22 = new TCanvas("pl22","pl22",900,700);
  	Energy_mucand_pile_ene_eMF->Draw();
  	//Energy_mucand_pile_ene_eMF->Fit(pol1_2);
    pl22->SaveAs("E_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl22->Close();

	cout << "" << endl;
	cout << "Jets Zr>1: " << endl;
	TCanvas *pl23 = new TCanvas("pl23","pl23",900,700);
  	Energy_jet_pile_ene_eMF->Draw();
  	//Energy_jet_pile_ene_eMF->Fit(pol1_3);
    pl23->SaveAs("E_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl23->Close();

	TCanvas *pl24 = new TCanvas("pl24","pl24",900,700);
  	Energy_pile_z1_ene->SetTitle("Energy of jets (Zr=1) with pile up (QF<" + SQFcut + ")");
  	Energy_pile_z1_ene->Draw();
  	Energy_pile_z1_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl24->SaveAs("E_pile_z1_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl24->Close();

	TCanvas *pl25 = new TCanvas("pl25","pl25",900,700);
  	Energy_mucand_pile_ene->SetTitle("Energy of muon candidates with pile up (QF<" + SQFcut + ")");
  	Energy_mucand_pile_ene->Draw();
  	Energy_mucand_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl25->SaveAs("E_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl25->Close();

	TCanvas *pl26 = new TCanvas("pl26","pl26",900,700);
  	Energy_jet_pile_ene->SetTitle("Energy of jets (Zr>1) with pile up (QF<" + SQFcut + ")");
  	Energy_jet_pile_ene->Draw();
  	Energy_jet_pile_eMF->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl26->SaveAs("E_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl26->Close();


	TCanvas *pl27 = new TCanvas("pl27","pl27",900,700);
  	Energy_alone_pile_z1_eMF->SetTitle("Energy of alone jets (Zr=1) with pile up (QF<" + SQFcut + ")");
  	Energy_alone_pile_z1_eMF->Draw();
  	Energy_alone_pile_z1_ene->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl27->SaveAs("E_alone_pile_z1_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl27->Close();

	TCanvas *pl28 = new TCanvas("pl28","pl28",900,700);
  	Energy_alone_mucand_pile_eMF->SetTitle("Energy of alone muon candidates with pile up (QF<" + SQFcut + ")");
  	Energy_alone_mucand_pile_eMF->Draw();
  	Energy_alone_mucand_pile_ene->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl28->SaveAs("E_alone_mucand_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl28->Close();

	TCanvas *pl29 = new TCanvas("pl29","pl29",900,700);
  	Energy_alone_jet_pile_eMF->SetTitle("Energy of alone jets (Zr>1) with pile up (QF<" + SQFcut + ")");
  	Energy_alone_jet_pile_eMF->Draw();
  	Energy_alone_jet_pile_ene->Draw("SAME");
  	leg4->Draw("SAME");
  	leg2->Draw("SAME");
    pl29->SaveAs("E_jet_pile_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
    pl29->Close();


	TCanvas *pl20 = new TCanvas("pl20","pl20",900,700);
	Cell_pileup_eMF_minus_ene->Draw();
	leg4->Draw("SAME");
	pl20->Close();	

	TCanvas *pl19 = new TCanvas("pl19","pl19",900,700);
	Cell_pileup_percent_ene->SetTitle("Percent of pileup cells per event (QF<"+SQFcut+")");
	Cell_pileup_percent_ene->Draw();
	Cell_pileup_percent_eMF->Draw("SAME");
	leg2->Draw("SAME");
	pl19->Close();

	TCanvas *pl17 = new TCanvas("pl17","pl17",900,700);
	Energy_alone_jet_pileup_ene->Draw();
	pl17->Close();

	TCanvas *plr18 = new TCanvas("plr18","plr18",900,700);
	Energy_alone_jet_pileup_eMF->Draw();
	plr18->Close();

	TCanvas *pl151 = new TCanvas("pl151","pl151",900,700);
	Ratio_jet_pile_energy->Draw();
	leg5->Draw();
	pl151->SaveAs("Ratio_jet_pile_energy_noisecutCOF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
	pl151->SetLogy();
	pl151->SaveAs("Ratio_jet_pile_energy_noisecut_LOG_COF"+Snoisecut_eMF+"QF" + SQFcut + ".pdf","pdf");
	pl151->Close();

	TCanvas *pl16 = new TCanvas("pl16","pl16",900,700);
	Zr_jetdif_pileup->Draw();
	pl16->Close();

    TCanvas *plr12 = new TCanvas("plr12","pl12",900,700);
	Nmuon_COFminusOF2_pileup->Draw();
	plr12->Close();

    TCanvas *pl13 = new TCanvas("pl13","pl13",900,700);
	Njets_COFminusOF2_pileup->Draw();
	pl13->Close();

    TCanvas *pl14 = new TCanvas("pl14","pl14",900,700);
	Njets_z1_COFminusOF2_pileup->Draw();
	pl14->Close();


	Double_t Dnopile = nopile; 
	Double_t Dpileevent = pileevent;
	Double_t Ddif_jets = dif_jets;
	Double_t Ddif_z1 = dif_z1;
	Double_t Ddif_zi = dif_zi;
	Double_t Ddif_zii = dif_zii;
	Double_t Ddif_muon = dif_muon;
	Double_t Dk_eMF_bigger = k_eMF_bigger;
	Double_t Dk_ene_bigger = k_ene_bigger;
	Double_t Dk_ene_eMF_equal = k_ene_eMF_equal;
	Double_t Dmuoncandidate_eMF_bigger = muoncandidate_eMF_bigger;
	Double_t Dmuoncandidate_ene_bigger = muoncandidate_ene_bigger;
	Double_t Dmuoncandidate_ene_eMF_equal = muoncandidate_ene_eMF_equal;
	Double_t Djet_pile_ene = jet_pile_ene;
	Double_t Dpile_z1_ene = pile_z1_ene;
	Double_t Dmucand_pile_ene = mucand_pile_ene;
	Double_t Djet_pile_eMF = jet_pile_eMF;
	Double_t Dpile_z1_eMF = pile_z1_eMF;
	Double_t Dmucand_pile_eMF = mucand_pile_eMF;
	Double_t Dall_pile_ene = all_pile_ene;
	Double_t Dall_pile_eMF = all_pile_eMF;

	cout << "Energy pile up cells" << endl;
	Double_t Doubpileevent = pileevent;
	cout << "Number of event with pile up = " << 100*Doubpileevent/Doubjentry << "%" << endl;
	cout << "Total cells = " << OF2COFall << endl;
	cout << "Total pile up cells = " << OF2COFpileup << endl;
	cout << "Events where COF has more jets than OF2 " << 100*Dk_eMF_bigger/(Dk_eMF_bigger+Dk_ene_bigger+Dk_ene_eMF_equal) << "%" << endl;
	cout << "Events where OF2 has more jets than COF " << 100*Dk_ene_bigger/(Dk_eMF_bigger+Dk_ene_bigger+Dk_ene_eMF_equal) << "%" << endl;
	cout << "Events where COF and OF2 have the same number of jets " << 100*Dk_ene_eMF_equal/(Dk_eMF_bigger+Dk_ene_bigger+Dk_ene_eMF_equal) << "%" << endl;
	cout << "Events where COF has more muon candidates than OF2 " << 100*Dmuoncandidate_eMF_bigger/(Dmuoncandidate_eMF_bigger+Dmuoncandidate_ene_bigger+Dmuoncandidate_ene_eMF_equal) << "%" << endl;
	cout << "Events where OF2 has more muon candidates than COF " << 100*Dmuoncandidate_ene_bigger/(Dmuoncandidate_eMF_bigger+Dmuoncandidate_ene_bigger+Dmuoncandidate_ene_eMF_equal) << "%" << endl;
	cout << "Events where COF and OF2 have the same number of muon candidates " << 100*Dmuoncandidate_ene_eMF_equal/(Dmuoncandidate_eMF_bigger+Dmuoncandidate_ene_bigger+Dmuoncandidate_ene_eMF_equal) << "%" << endl;
	cout << " " << endl;
	cout << "Events without pileup cells " << 100*Dnopile/(Dnopile + Dpileevent) << "%" << endl;	
	cout << "Events with pileup cells " << 100*Dpileevent/(Dnopile + Dpileevent) << "%"<< endl;
	cout << "Total = " << nopile + pileevent << endl;	 
	cout << " " << endl;
	cout << "Eventos pileup que diferem em nº de múons: " << 100*Ddif_muon/(Dnopile + Dpileevent) << "%" << endl;
	cout << "Eventos pileup que diferem em nº de jatos zr=1: " << 100*Ddif_z1/(Dnopile + Dpileevent) << "%" << endl;
	cout << "Eventos pileup que diferem em nº de jatos: " << 100*Ddif_jets/(Dnopile + Dpileevent) << "%" << endl;
	cout << "Eventos pileup que diferem em nº de jatos zr>1: " << 100*Ddif_zi/(Dnopile + Dpileevent) << "%" << endl;
	cout << "Eventos pileup que diferem em nº de jatos zr>2: " << 100*Ddif_zii/(Dnopile + Dpileevent) << "%" << endl;
	cout << " " << endl;


	cout << "noisecut COF: " << noisecut_eMF << " MeV ; noise cut OF2: " << noisecut_ene << " MeV" << endl;
	cout << "OF2 - Where the jets (Zr>1) have pile cells: " << Djet_pile_ene << " cells (all events), " << 100*Djet_pile_ene/Dall_pile_ene << "% of cells with energy > " << noisecut_ene << " in all events" << endl;
	cout << "COF - Where the jets (Zr>1) have pile cells: " << Djet_pile_eMF << " cells (all events), " << 100*Djet_pile_eMF/Dall_pile_eMF << "% of cells with energy > " << noisecut_eMF << " in all events" << endl;
	cout << "OF2 - Where the jets (Zr=1) have pile cells: " << Dpile_z1_ene << " cells (all events), " << 100*Dpile_z1_ene/Dall_pile_ene << "% of cells with energy > " << noisecut_ene << " in all events" << endl;
	cout << "COF - Where the jets (Zr=1) have pile cells: " << Dpile_z1_eMF << " cells (all events), " << 100*Dpile_z1_eMF/Dall_pile_eMF << "% of cells with energy > " << noisecut_eMF << " in all events" << endl;
	cout << "OF2 - Where the muon candidates have pile cells: " << Dmucand_pile_ene << " cells (all events), " << 100*Dmucand_pile_ene/Dall_pile_ene << "% of cells with energy > " << noisecut_ene << " in all events" << endl;
	cout << "COF - Where the muon candidates have pile cells: " << Dmucand_pile_eMF << " cells (all events), " << 100*Dmucand_pile_eMF/Dall_pile_eMF << "% of cells with energy > " << noisecut_eMF << " in all events" << endl;


	Double_t pileup_ = OF2COFpileup;
	Double_t ratioup1down2 = OF2COFup1down2;
	Double_t ratioup2down3 = OF2COFup2down3;
	Double_t ratioup3 = OF2COFup3;
	Double_t ratioup0down1 = OF2COFup0down1;
	Double_t ratioup0down05 = OF2COFup0down05;
	Double_t ratioup05down1 = OF2COFup05down1;
	Double_t ratioup_1down0 = OF2COFup_1down0;
	Double_t ratioup_2down_1 = OF2COFup_2down_1;
	Double_t ratioup_3down_2 = OF2COFup_3down_2;
	Double_t ratiodown_3 = OF2COFdown_3;
	Double_t ratioup05down15 = OF2COFup05down15;

	cout << " " << endl;
	cout << "OF2/COF > 0.5 and OF2/COF < 1.5 = " << 100*ratioup05down15/pileup_ << "%" << endl;
	cout << "OF2/COF > 1 and OF2/COF < 2 = " << 100*ratioup1down2/pileup_ << "%" << endl;
	cout << "OF2/COF > 2 and OF2/COF < 3 = " << 100*ratioup2down3/pileup_ << "%" << endl;
	cout << "OF2/COF > 3 = " << 100*ratioup3/pileup_ << "%" << endl;
	cout << "OF2/COF > 0 and OF2/COF < 1 = " << 100*ratioup0down1/pileup_ << "%" << endl;
	cout << "OF2/COF > 0 and OF2/COF < 0.5 = " << 100*ratioup0down05/pileup_ << "%" << endl;
	cout << "OF2/COF > 0.5 and OF2/COF < 1 = " << 100*ratioup05down1/pileup_ << "%" << endl;
	cout << "OF2/COF > -1 and OF2/COF < 0 = " << 100*ratioup_1down0/pileup_ << "%" << endl;
	cout << "OF2/COF > -2 and OF2/COF < -1 = " << 100*ratioup_2down_1/pileup_ << "%" << endl;
	cout << "OF2/COF > -3 and OF2/COF < -2 = " << 100*ratioup_3down_2/pileup_ << "%" << endl;
	cout << "OF2/COF < -3 = " << 100*ratiodown_3/pileup_ << "%" << endl;

	cout << "  " << endl;
	Double_t Dall_rec_cell_ene = all_rec_cell_ene;
	Double_t Dall_Edif0_cell_ene = all_Edif0_cell_ene;
	Double_t Dall_rec_cell_eMF = all_rec_cell_eMF;
	Double_t Dall_Edif0_cell_eMF = all_Edif0_cell_eMF;
	cout << "Do total de células com EOF2!=0, quantos % são E_OF2 > " << noisecut_ene << " e E_COF < " << noisecut_eMF << " = " << 100*Dall_rec_cell_ene/Dall_Edif0_cell_ene << "%" << endl;
	cout << "Do total de células com EOF2!=0, quantos % são E_COF > " << noisecut_eMF << " e E_OF2 < " << noisecut_ene << " = " << 100*Dall_rec_cell_eMF/Dall_Edif0_cell_eMF << "%" << endl;
	cout << "  " << endl;
	
	Double_t DRatio_cell_all = Ratio_cell_all;
	Double_t DRatio_cell_0to05 = Ratio_cell_0to05;
	Double_t DRatio_cell_05to15 = Ratio_cell_05to15;
	Double_t DRatio_cell_15to3 = Ratio_cell_15to3;	
	Double_t DRatio_cell_3to10 = Ratio_cell_3to10;
	Double_t DRatio_cell_justCOF = Ratio_cell_justCOF; 
	Double_t DRatio_cell_justOF2 = Ratio_cell_justOF2; 

	Double_t DRatio_tower_all = Ratio_tower_all;
	Double_t DRatio_tower_0to05 = Ratio_tower_0to05;
	Double_t DRatio_tower_05to15 = Ratio_tower_05to15;
	Double_t DRatio_tower_15to3 = Ratio_tower_15to3;	
	Double_t DRatio_tower_3to10 = Ratio_tower_3to10;
	Double_t DRatio_tower_justCOF = Ratio_tower_justCOF; 
	Double_t DRatio_tower_justOF2 = Ratio_tower_justOF2; 

	cout << "No histograma Ratio_cell:" << endl;

	cout << "ECOF/EOF2 >= 0 && ECOF/EOF2 < 0.5 = " << 100*DRatio_cell_0to05/DRatio_cell_all << "%" << endl;
	cout << "ECOF/EOF2 >= 0.5 && ECOF/EOF2 < 1.5 = " << 100*DRatio_cell_05to15/DRatio_cell_all << "%" << endl;
	cout << "ECOF/EOF2 >= 1.5 && ECOF/EOF2 < 3.0 = " << 100*DRatio_cell_15to3/DRatio_cell_all << "%" << endl;
	cout << "ECOF/EOF2 >= 3.0 && ECOF/EOF2 < 10.0 = " << 100*DRatio_cell_3to10/DRatio_cell_all << "%" << endl;
	cout << "Just COF reconstruct = " << 100*DRatio_cell_justCOF/DRatio_cell_all << "%" << endl;
	cout << "Just OF2 reconstruct = " << 100*DRatio_cell_justOF2/DRatio_cell_all << "%" << endl;

	cout << " " << endl;

	cout << "No histograma Ratio_tower:" << endl;

	cout << "ECOF/EOF2 >= 0 && ECOF/EOF2 < 0.5 = " << 100*DRatio_tower_0to05/DRatio_tower_all << "%" << endl;
	cout << "ECOF/EOF2 >= 0.5 && ECOF/EOF2 < 1.5 = " << 100*DRatio_tower_05to15/DRatio_tower_all << "%" << endl;
	cout << "ECOF/EOF2 >= 1.5 && ECOF/EOF2 < 3.0 = " << 100*DRatio_tower_15to3/DRatio_tower_all << "%" << endl;
	cout << "ECOF/EOF2 >= 3.0 && ECOF/EOF2 < 10.0 = " << 100*DRatio_tower_3to10/DRatio_tower_all << "%" << endl;
	cout << "Just COF reconstruct = " << 100*DRatio_tower_justCOF/DRatio_tower_all << "%" << endl;
	cout << "Just OF2 reconstruct = " << 100*DRatio_tower_justOF2/DRatio_tower_all << "%" << endl;


	//Pileup Analysis

    TCanvas *pile1 = new TCanvas("pile1","pile1",900,700);
    Ecell_pileup_ene_eMF->Draw(); 
    pile1->SaveAs("Ecell_pileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile1->Close();

    TCanvas *pile2 = new TCanvas("pile2","pile2",900,700);
    Ecell_nopileup_ene_eMF->Draw(); 
    pile2->SaveAs("Ecell_nopileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile2->Close();

    TCanvas *pile3 = new TCanvas("pile3","pile3",900,700);
    Etower_pileup_ene_eMF->Draw(); 
    pile3->SaveAs("Etower_pileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile3->Close();

    TCanvas *pile4 = new TCanvas("pile4","pile4",900,700);
    Etower_nopileup_ene_eMF->Draw(); 
    pile4->SaveAs("Etower_nopileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile4->Close();

    TCanvas *pile5 = new TCanvas("pile5","pile5",900,700);
    Ejet_pileup_ene_eMF->Draw(); 
    pile5->SaveAs("Ejet_pileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile5->Close();

    TCanvas *pile6 = new TCanvas("pile6","pile6",900,700);
    Ejet_nopileup_ene_eMF->Draw(); 
    pile6->SaveAs("Ejet_nopileup_ene_eMF_QF" + SQFcut + ".pdf","pdf");
	//pile6->Close();

    TCanvas *pile7 = new TCanvas("pile7","pile7",900,700);
    Ratio_cell_pile_energy->Draw();
    leg5->Draw("SAME");
    pile7->SaveAs("Ratio_cell_pile_energy_QF" + SQFcut + ".pdf","pdf");
	//pile7->Close();

    TCanvas *pile8 = new TCanvas("pile8","pile8",900,700);
    Ratio_cell_nopile_energy->Draw();
    leg5->Draw("SAME");
    pile8->SaveAs("Ratio_cell_nopile_energy_QF" + SQFcut + ".pdf","pdf");
    //pile8->Close();

    TCanvas *pile9 = new TCanvas("pile9","pile9",900,700);
    Ratio_tower_pile_energy->Draw();
    leg5->Draw("SAME");
    pile9->SaveAs("Ratio_tower_pile_energy_QF" + SQFcut + ".pdf","pdf");
	//pile9->Close();

    TCanvas *pile10 = new TCanvas("pile10","pile10",900,700);
    Ratio_tower_nopile_energy->Draw();
    leg5->Draw("SAME");
    pile10->SaveAs("Ratio_tower_nopile_energy_QF" + SQFcut + ".pdf","pdf");
	//pile10->Close();

    TCanvas *pile11 = new TCanvas("pile11","pile11",900,700);
    Ratio_jet_pile_energy->Draw();
    leg5->Draw("SAME");
    pile11->SaveAs("Ratio_cell_pile_energy_QF" + SQFcut + ".pdf","pdf");
	//pile11->Close();

    TCanvas *pile123 = new TCanvas("pile123","pile123",900,700);
    Ratio_jet_nopile_energy->Draw();
    leg5->Draw("SAME");
    pile123->SaveAs("Ratio_jet_nopile_energy_QF" + SQFcut + ".pdf","pdf");
	//pile123->Close();


    TCanvas *pile133 = new TCanvas("pile133","pile133",900,700);
    Ecell_pileup_OF2alone_ene_eMF->Draw();
    pile133->SaveAs("Ecell_pileup_OF2alone_ene_eMF_QF" + SQFcut + ".pdf","pdf");
    //pile133->Close();

    TCanvas *pile143 = new TCanvas("pile123","pile133",900,700);
    Ecell_pileup_COFalone_eMF_ene->Draw();
    pile143->SaveAs("Ecell_pileup_COFalone_eMF_ene_QF" + SQFcut + ".pdf","pdf");
    //pile143->Close();
	*/

	/*
   	int bin1_cellpile = Ratio_cell_pile_energy->FindFirstBinAbove(Ratio_cell_pile_energy->GetMaximum()/2);
   	int bin2_cellpile = Ratio_cell_pile_energy->FindLastBinAbove(Ratio_cell_pile_energy->GetMaximum()/2);
   	double fwhm_cellpile = Ratio_cell_pile_energy->GetBinCenter(bin2_cellpile) - Ratio_cell_pile_energy->GetBinCenter(bin1_cellpile);
   	int bin1_cellnopile = Ratio_cell_nopile_energy->FindFirstBinAbove(Ratio_cell_nopile_energy->GetMaximum()/2);
   	int bin2_cellnopile = Ratio_cell_nopile_energy->FindLastBinAbove(Ratio_cell_nopile_energy->GetMaximum()/2);
   	double fwhm_cellnopile = Ratio_cell_nopile_energy->GetBinCenter(bin2_cellnopile) - Ratio_cell_nopile_energy->GetBinCenter(bin1_cellnopile);
	cout << " " << endl;
	Ratio_cell_pile_energy->GetXaxis()->SetRange(15,115);
	Ratio_cell_nopile_energy->GetXaxis()->SetRange(15,115);
	cout << "Ratio_cell_pile_energy FWHM = " << fwhm_cellpile << "; GetMaximum =  " << Ratio_cell_pile_energy->GetMaximumBin() << "; MOP = " << Ratio_cell_pile_energy->GetBinCenter(Ratio_cell_pile_energy->GetMaximumBin()) << endl;
	cout << "Ratio_cell_nopile_energy FWHM = " << fwhm_cellnopile << "; GetMaximum =  " << Ratio_cell_nopile_energy->GetMaximumBin() << "; MOP = " << Ratio_cell_nopile_energy->GetBinCenter(Ratio_cell_nopile_energy->GetMaximumBin()) << endl;
	*/


	/*
    TCanvas *pl15 = new TCanvas("pl15","pl15",900,700);
    Ejets_pileup_eMF->SetTitle("Energy of jets - Pileup Events (QF<"+SQFcut+")");
    Ejets_pileup_eMF->Draw();
    Ejets_pileup_ene->Draw("SAME");
    leg2->Draw("SAME");
    leg4->Draw("SAME");
    pl15->SaveAs("Ejets_pileup_noisecutCOF"+Snoisecut_eMF+"_QF" + SQFcut + ".pdf","pdf");
    pl15->Close();

    TCanvas *pl222 = new TCanvas("pl222","pl222",900,700);
    Ejets_pileup_ene_eMF->SetTitle("Energy of jets - Pileup Events (noisecutCOF="+Snoisecut_eMF+"QF<"+SQFcut+")");
    Ejets_pileup_ene_eMF->Draw();
    pl222->SaveAs("Ejets_pileup_noisecutCOF"+Snoisecut_eMF+"_QF" + SQFcut + "2.pdf","pdf");
    pl222->Close();

    TCanvas *pl333 = new TCanvas("pl333","pl333",900,700);
    Ejets_alone_pileup_eMF->SetTitle("Energy of alone jets - Pileup Events (QF<"+SQFcut+")");
    Ejets_alone_pileup_eMF->Draw();
    Ejets_alone_pileup_ene->Draw("SAME");
    leg2->Draw("SAME");
    leg4->Draw("SAME");
    pl333->SaveAs("Ejets_alone_pileup_noisecutCOF"+Snoisecut_eMF+"_QF" + SQFcut + ".pdf","pdf");
    pl333->Close();

    TCanvas *plr16 = new TCanvas("plr16","plr16",900,700);
    Ejets_pileup_muons_eMF->SetTitle("Energy muons candidates - Pileup (QF<"+SQFcut+")");
    Ejets_pileup_muons_eMF->Draw();
    Ejets_pileup_muons_ene->Draw("SAME");
    leg2->Draw("SAME");
    plr16->SaveAs("Ejets_pileupall_muons_QF" + SQFcut + ".pdf","pdf");
    plr16->Close();

    TCanvas *plr25 = new TCanvas("plr25","plr25",900,700);
    Ejets_pileup_z1_eMF->SetTitle("Energy jets Zr = 1 - Pileup (QF<"+SQFcut+")");
    Ejets_pileup_z1_eMF->Draw();
    Ejets_pileup_z1_ene->Draw("SAME");
    leg2->Draw("SAME");
    plr25->SaveAs("Ejets_pileupall_z1_QF" + SQFcut + ".pdf","pdf");
    plr25->Close();

    TCanvas *plr17 = new TCanvas("plr17","plr17",900,700);
    Ejets_pileup_muons->SetTitle("Energy muons candidates - Pileup (QF<"+SQFcut+")");
    Ejets_pileup_muons->Draw();
    plr17->SaveAs("Ejets_pileupall2_muons_QF" + SQFcut + ".pdf","pdf");
    plr17->Close();

    TCanvas *pl18 = new TCanvas("pl18","pl18",900,700);
    Ejets_pileup_muons_enexeMF->SetTitle("Energy muons candidates - Pileup (QF<"+SQFcut+")");
    Ejets_pileup_muons_enexeMF->Draw();
    pl18->SaveAs("Ejets_pileupall3_muons_QF" + SQFcut + ".pdf","pdf");
    pl18->Close();

	*/
	/*
	//Plotting Emu_ene/Emu_eMF in a TGraph

    TCanvas *p13 = new TCanvas("p13","p13",200,10,700,500);
	Double_t bin_ene[100], bin_eMF[100], bin[100], bin_ene_error[100], bin_eMF_error[100], error[100], ratio[100];

    //p13->SetFillColor(42);
    p13->SetGrid();
    //p13->GetFrame()->SetFillColor(21);
    p13->GetFrame()->SetBorderSize(12);
	int cont;
	//cout << "error bin number: 3, " << Ejets_mu_cand_ene->GetBinError(3) << "." << endl;
	for (cont=0; cont < 101; cont++){
		bin_ene[cont] = Ejets_mu_cand_ene->GetBinContent(cont+1);
		bin_eMF[cont] = Ejets_mu_cand_eMF->GetBinContent(cont+1);
		//cout << bin_ene[cont] << endl;
		bin_ene_error[cont] = sqrt(bin_ene[cont]);
		bin_eMF_error[cont] = sqrt(bin_eMF[cont]);
		//cout << bin_ene_error[cont] << endl;
		if (bin_eMF[cont] != 0 && bin_ene[cont] != 0) ratio[cont] = bin_ene[cont]/bin_eMF[cont];
		if (bin_eMF[cont] == 0 || bin_ene[cont] == 0) ratio[cont] = 0;
		error[cont] = ratio[cont];
		if (bin_eMF[cont] != 0 && bin_ene[cont] != 0) error[cont] = ratio[cont]*sqrt(1/bin_ene[cont]+1/bin_eMF[cont]);
		//cout << error[cont] << endl;
		//cout << "error bin number: " << cont << ", " << error[cont] << endl;
		//Energy_mu_ene_eMF->Fill(count,bin_ene[cont]_eMF[cont]);	
	}

	bin[0]=0;
	bin[1]=100;
	bin[2]=200;
	bin[3]=300;
	bin[4]=400;
	bin[5]=500;
	bin[6]=600;
	bin[7]=700;
	bin[8]=800;
	bin[9]=900;
	bin[10]=1000;
	bin[11]=1100;
	bin[12]=1200;
	bin[13]=1300;
	bin[14]=1400;
	bin[15]=1500;
	bin[16]=1600;
	bin[17]=1700;
	bin[18]=1800;
	bin[19]=1900;
	bin[20]=2000;
	bin[21]=2100;
	bin[22]=2200;
	bin[23]=2300;
	bin[24]=2400;
	bin[25]=2500;
	bin[26]=2600;
	bin[27]=2700;
	bin[28]=2800;
	bin[29]=2900;
	bin[30]=3000;
	bin[31]=3100;
	bin[32]=3200;
	bin[33]=3300;
	bin[34]=3400;
	bin[35]=3500;

	Double_t bin_e[35] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	//error[10] = 0.1;
	//error[11] = 0.25;
	for (count=36; cont < 100; cont++){
		bin[cont] = cont*100;
		//bin_e[cont] = 0;
	}

	

	cont = 100;
  	//TGraph *gr = new TGraph(cont,bin,ratio);
 	TGraphErrors *gr = new TGraphErrors(cont,bin,ratio,bin_e,error);
	gr->GetXaxis()->SetTitle("Energy (MeV)");
	gr->GetYaxis()->SetTitle("OF2/COF");
  	gr->SetTitle("Energy muon candidates OF2/COF bin per bin");
  	gr->SetMarkerColor(4);
   	gr->SetMarkerStyle(21);
    //TFile *rat = new TFile("EMu_ratio_QF" + SQFcut + ".root","RECREATE");
  	//gr->Write();
	//rat->Close();
	gr->Draw("ap");
    //Energy_mu_ene_eMF->Draw(); 
	//p13->Update();
    //p13->SaveAs("E.pdf","pdf");
    //p13->Close();
	*/

    //Create the legend
    leg = new TLegend(0.7,0.65,0.8,0.8);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg->SetBorderSize(0);
    //leg->SetHeader("Legend");
    leg->AddEntry(Ejets_eMF,"COF","l");
    leg->AddEntry(Ejets_ene,"OF2","l");
    leg->SetTextSize(0.035);

	//Fit the Energy_muon, Energy_muon pileup and Energy/cm plots (Gauss, Landau and Pol9)
/*
    TCanvas *f719 = new TCanvas("f719","f719",900,700);
    Emm_all_nocut_ene->SetTitle("Energy per cm (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Emm_all_nocut_ene->Draw(); 
    Emm_all_nocut_eMF->Draw("SAME"); 
	f719->SaveAs("QF" + SQFcut + "/Emm_all_nocut_QF" + SQFcut + ".pdf","pdf");
	f719->SaveAs("QF" + SQFcut + "/Emm_all_nocut_QF" + SQFcut + ".jpeg","jpeg");
    f719->Close();

*/	//AXION PARTICLE SEARCH------------------------------------ 
	/*Ephoton<" + SEphoton + "MeV,*/

/*
    TCanvas *f741 = new TCanvas("f741","f741",900,700);
    MASS_Particle_ene->SetTitle("Invariant mass gamma+gamma  (Ephoton<" + SEphoton + "MeV, E just layer 1, 2 or 3, Layer point "+layerpoint+" Zr=1+1cell) QF<" + SQFcut + ")");
    MASS_Particle_ene->Draw(); 
    MASS_Particle_eMF->Draw("SAME"); 
    leg->Draw("SAME");
	f741->SaveAs("QF" + SQFcut + "/Axion/MASS_ParticleEphoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".pdf","pdf");
	f741->SaveAs("QF" + SQFcut + "/Axion/MASS_ParticleEphoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".jpeg","jpeg");
    //f741->Close();

    //Ephoton<" + SEphoton + "MeV, 

    TCanvas *f7411 = new TCanvas("f7411","f7411",900,700);
    MASS_Particle_l1_ene->SetTitle("Invariant mass gamma+gamma (Ephoton<" + SEphoton + "MeV, E just layer 1, Layer point "+layerpoint+" Zr=1+1cell) QF<" + SQFcut + ")");
    MASS_Particle_l1_ene->Draw(); 
    MASS_Particle_l1_eMF->Draw("SAME"); 
    leg->Draw("SAME");
	f7411->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l1Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".pdf","pdf");
	f7411->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l1Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".jpeg","jpeg");
    //f7411->Close();

    TCanvas *f7412 = new TCanvas("f7412","f7412",900,700);
    MASS_Particle_l2_ene->SetTitle("Invariant mass gamma+gamma (Ephoton<" + SEphoton + "MeV, E just layer 2, Layer point "+layerpoint+" Zr=1+1cell) QF<" + SQFcut + ")");
    MASS_Particle_l2_ene->Draw(); 
    MASS_Particle_l2_eMF->Draw("SAME"); 
    leg->Draw("SAME");
	f7412->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l2Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".pdf","pdf");
	f7412->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l2Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".jpeg","jpeg");
    //f7412->Close();

    TCanvas *f7413 = new TCanvas("f7413","f7413",900,700);
    MASS_Particle_l3_eMF->SetTitle("Invariant mass gamma+gamma (Ephoton<" + SEphoton + "MeV, E just layer 3, Layer point "+layerpoint+" Zr=1+1cell) QF<" + SQFcut + ")");
    MASS_Particle_l3_eMF->Draw(); 
    MASS_Particle_l3_ene->Draw("SAME"); 
    leg->Draw("SAME");
	f7413->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l3Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".pdf","pdf");
	f7413->SaveAs("QF" + SQFcut + "/Axion/MASS_Particle_l3Ephoton<" + SEphoton + "Layerpoint"+layerpoint+"QF" + SQFcut + ".jpeg","jpeg");
    //f7413->Close();
*/
/*
	//-------------------------------------------------------------------------------------------------------------------------------------
    TCanvas *f731 = new TCanvas("f731","f731",900,700);
    N_channel_per_jet_Energy_ene->Draw(); 
	f731->SaveAs("N_channel_per_jet_Energy_ene_QF" + SQFcut + ".pdf","pdf");
	f731->SaveAs("N_channel_per_jet_Energy_ene_QF" + SQFcut + ".jpeg","jpeg");
    f731->Close();

    TCanvas *f732 = new TCanvas("f732","f732",900,700);
    N_channel_per_jet_Energy_eMF->Draw(); 
	f732->SaveAs("N_channel_per_jet_Energy_eMF_QF" + SQFcut + ".pdf","pdf");
	f732->SaveAs("N_channel_per_jet_Energy_eMF_QF" + SQFcut + ".jpeg","jpeg");
    f732->Close();

    TCanvas *f721 = new TCanvas("f721","f721",900,700);
    N_channel_per_jet_eMF->SetTitle("Number of channels with energy per jet (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    N_channel_per_jet_eMF->Draw(); 
    N_channel_per_jet_ene->Draw("SAME"); 
	f721->SaveAs("N_channel_per_jetQF" + SQFcut + ".pdf","pdf");
	f721->SaveAs("N_channel_per_jetQF" + SQFcut + ".jpeg","jpeg");
    f721->Close();

    TCanvas *f711 = new TCanvas("f711","f711",900,700);
    Emm_nomuon_l12_ene->Draw(); 
	f711->SaveAs("Emm_nomuon_l12_ene_QF" + SQFcut + ".pdf","pdf");
	f711->SaveAs("Emm_nomuon_l12_ene_QF" + SQFcut + ".jpeg","jpeg");
    //f711->Close();

    TCanvas *f712 = new TCanvas("f712","f712",900,700);
    Emm_nomuon_l12_eMF->Draw(); 
	f712->SaveAs("Emm_nomuon_l12_eMF_QF" + SQFcut + ".pdf","pdf");
	f712->SaveAs("Emm_nomuon_l12_eMF_QF" + SQFcut + ".jpeg","jpeg");
    f712->Close();

    TCanvas *f713 = new TCanvas("f713","f713",900,700);
    Emm_nomuon_l13_ene->Draw(); 
	f713->SaveAs("Emm_nomuon_l13_ene_QF" + SQFcut + ".pdf","pdf");
	f713->SaveAs("Emm_nomuon_l13_ene_QF" + SQFcut + ".jpeg","jpeg");
    f713->Close();

    TCanvas *f714 = new TCanvas("f714","f714",900,700);
    Emm_nomuon_l13_eMF->Draw(); 
	f714->SaveAs("Emm_nomuon_l13_eMF_QF" + SQFcut + ".pdf","pdf");
	f714->SaveAs("Emm_nomuon_l13_eMF_QF" + SQFcut + ".jpeg","jpeg");
    f714->Close();

    TCanvas *f715 = new TCanvas("f715","f715",900,700);
    Emm_nomuon_l23_ene->Draw(); 
	f715->SaveAs("Emm_nomuon_l23_ene_QF" + SQFcut + ".pdf","pdf");
	f715->SaveAs("Emm_nomuon_l23_ene_QF" + SQFcut + ".jpeg","jpeg");
    f715->Close();

    TCanvas *f716 = new TCanvas("f716","f716",900,700);
    Emm_nomuon_l23_eMF->Draw(); 
	f716->SaveAs("Emm_nomuon_l23_eMF_QF" + SQFcut + ".pdf","pdf");
	f716->SaveAs("Emm_nomuon_l23_eMF_QF" + SQFcut + ".jpeg","jpeg");
    f714->Close();

    TCanvas *f71 = new TCanvas("f71","f71",900,700);
    E_nomuon_l1_eMF->SetTitle("Energy of no-muons candidates layer 1 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    E_nomuon_l1_eMF->Draw(); 
    E_nomuon_l1_ene->Draw("SAME"); 
	f71->SaveAs("E_nomuon_l1_QF" + SQFcut + ".pdf","pdf");
	f71->SaveAs("E_nomuon_l1_QF" + SQFcut + ".jpeg","jpeg");
    f71->Close();

    TCanvas *f72 = new TCanvas("f72","f72",900,700);
    E_nomuon_l2_ene->SetTitle("Energy of no-muons candidates layer 2 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    E_nomuon_l2_ene->Draw(); 
    E_nomuon_l2_eMF->Draw("SAME"); 
	f72->SaveAs("E_nomuon_l2_QF" + SQFcut + ".pdf","pdf");
	f72->SaveAs("E_nomuon_l2_QF" + SQFcut + ".jpeg","jpeg");
    f72->Close();

    TCanvas *f73 = new TCanvas("f73","f73",900,700);
    E_nomuon_l3_eMF->SetTitle("Energy of no-muons candidates layer 3 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    E_nomuon_l3_eMF->Draw(); 
    E_nomuon_l3_ene->Draw("SAME"); 
	f73->SaveAs("E_nomuon_l3_QF" + SQFcut + ".pdf","pdf");
	f73->SaveAs("E_nomuon_l3_QF" + SQFcut + ".jpeg","jpeg");
    f73->Close();

    TCanvas *f74 = new TCanvas("f74","f74",900,700);
    Emm_nomuon_l1_ene->SetTitle("Energy per cm of no-muons candidates layer 1 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Emm_nomuon_l1_ene->Draw(); 
    Emm_nomuon_l1_eMF->Draw("SAME"); 
	f74->SaveAs("Emm_nomuon_l1_QF" + SQFcut + ".pdf","pdf");
	f74->SaveAs("Emm_nomuon_l1_QF" + SQFcut + ".jpeg","jpeg");
    f74->Close();

    TCanvas *f75 = new TCanvas("f75","f75",900,700);
    Emm_nomuon_l2_ene->SetTitle("Energy per cm of no-muons candidates layer 2 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Emm_nomuon_l2_ene->Draw(); 
    Emm_nomuon_l2_eMF->Draw("SAME"); 
	f75->SaveAs("Emm_nomuon_l2_QF" + SQFcut + ".pdf","pdf");
	f75->SaveAs("Emm_nomuon_l2_QF" + SQFcut + ".jpeg","jpeg");
    f75->Close();

    TCanvas *f76 = new TCanvas("f76","f76",900,700);
    Emm_nomuon_l3_eMF->SetTitle("Energy per cm of no-muons candidates layer 3 (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Emm_nomuon_l3_eMF->Draw(); 
    Emm_nomuon_l3_ene->Draw("SAME"); 
	f76->SaveAs("Emm_nomuon_l3_QF" + SQFcut + ".pdf","pdf");
	f76->SaveAs("Emm_nomuon_l3_QF" + SQFcut + ".jpeg","jpeg");
    f76->Close();
	//---------------------------------------------------------
*/
	/*
	//MASS_PARTICLE FIT - AXION PARTICLE
	landau2_ene = new TF1("landau2_ene","landau",0,100);
	gauss3_ene = new TF1("gauss3_ene","gaus",0,800);
	landau2_eMF = new TF1("landau2_eMF","landau",0,100);
	gauss3_eMF = new TF1("gauss3_eMF","gaus",0,800);
	total2_ene = new TF1("mstotal","landau(0)+gaus(3)",0,800);
	total2_eMF = new TF1("mstotal","landau(0)+gaus(3)",0,800);
	landau2_ene->SetLineColor(1);
	gauss3_ene->SetLineColor(1);
	landau2_eMF->SetLineColor(1);
	gauss3_eMF->SetLineColor(1);


    TCanvas *fit11 = new TCanvas("fit11","fit11",900,700);
    MASS_Particle_ene->SetTitle("OF2 - Mass of particle (E just layer 1, 2 towers, theta = "+Sang_var+") QF<" + SQFcut + ")");
    MASS_Particle_ene->Draw(); 
	MASS_Particle_ene->Fit(landau2_ene,"R"); //$
	//MASS_Particle_ene->Fit(gauss3_ene,"R+"); //$
	Double_t par_ene6[14];
	//landau2_ene->GetParameters(&par_ene6[0]);
	gauss3_ene->GetParameters(&par_ene6[0]);
	total2_ene->SetParameters(par_ene6);
	//Ejets_mu_cand_ene->Fit(total_ene); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit11->SaveAs("MASS_Particle_landau-gauss-110-800_eneQF" + SQFcut + ".pdf","pdf");
	fit11->SaveAs("MASS_Particle_landau-gauss-110-800_eneQF" + SQFcut + ".jpeg","jpeg");
    fit11->Close();

    TCanvas *fit12 = new TCanvas("fit12","fit12",900,700);
    MASS_Particle_eMF->SetTitle("COF - Mass of particle (E just layer 1, 2 towers, theta = "+Sang_var+") QF<" + SQFcut + ")");
    MASS_Particle_eMF->Draw(); 
	MASS_Particle_eMF->Fit(landau2_eMF,"R"); //$
	MASS_Particle_eMF->Fit(gauss3_eMF,"R+"); //$
	Double_t par_eMF6[14];
	//landau2_eMF->GetParameters(&par_eMF6[0]);
	gauss3_eMF->GetParameters(&par_eMF6[0]);
	total2_eMF->SetParameters(par_eMF6);
	//Ejets_mu_cand_ene->Fit(total_ene); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit12->SaveAs("MASS_Particle_landau-gauss-110-800_eMFQF" + SQFcut + ".pdf","pdf");
	fit12->SaveAs("MASS_Particle_landau-gauss-110-800_eMFQF" + SQFcut + ".jpeg","jpeg");
    fit12->Close();

	//ENERGY_MUON LANDAU+2GAUSS

	landau_ene = new TF1("landau_ene","landau",0,600);
	gauss1_ene = new TF1("gauss1_ene","gaus",600,1000);
	gauss2_ene = new TF1("gauss2_ene","gaus",1000,3500);
	total_ene = new TF1("mstotal","landau(0)+gaus(3)+gaus(6)",0,3500);
	landau_eMF = new TF1("landau_eMF","landau",0,600);
	gauss1_eMF = new TF1("gauss1_eMF","gaus",600,1000);
	gauss2_eMF = new TF1("gauss2_eMF","gaus",1000,3500);
	total_eMF = new TF1("mstotal","landau(0)+gaus(3)+gaus(6)",0,3000);
	landau_ene->SetLineColor(1);
	gauss1_ene->SetLineColor(1);
	gauss2_ene->SetLineColor(1);
	total_ene->SetLineColor(1);
	landau_eMF->SetLineColor(1);
	gauss1_eMF->SetLineColor(1);
	gauss2_eMF->SetLineColor(1);
	total_eMF->SetLineColor(1);

    TCanvas *fit7 = new TCanvas("fit7","fit7",900,700);
    Ejets_mu_cand_ene->SetTitle("OF2 - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_ene->Draw(); 
	Ejets_mu_cand_ene->Fit(landau_ene,"R"); //$
	Ejets_mu_cand_ene->Fit(gauss1_ene,"R+"); //$
	Ejets_mu_cand_ene->Fit(gauss2_ene,"R+"); //$
	Double_t par_ene2[14];
	landau_ene->GetParameters(&par_ene2[0]);
	gauss1_ene->GetParameters(&par_ene2[3]);
	gauss2_ene->GetParameters(&par_ene2[6]);
	total_ene->SetParameters(par_ene2);
	//Ejets_mu_cand_ene->Fit(total_ene); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit7->SaveAs("Ejets_mucand_landau-2gauss-600-1000_eneQF" + SQFcut + ".pdf","pdf");
	fit7->SaveAs("Ejets_mucand_landau-2gauss-600-1000_eneQF" + SQFcut + ".jpeg","jpeg");
        fit7->Close();

    TCanvas *fit8 = new TCanvas("fit8","fit8",900,700);
    Ejets_mu_cand_eMF->SetTitle("COF - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_eMF->Draw(); 
	Ejets_mu_cand_eMF->Fit(landau_eMF,"R"); //$
	Ejets_mu_cand_eMF->Fit(gauss1_eMF,"R+"); //$
	Ejets_mu_cand_eMF->Fit(gauss2_eMF,"R+"); //$
	Double_t par_eMF2[14];
	landau_eMF->GetParameters(&par_eMF2[0]);
	gauss1_eMF->GetParameters(&par_eMF2[3]);
	gauss2_eMF->GetParameters(&par_eMF2[6]);
	total_eMF->SetParameters(par_eMF2);
	//Ejets_mu_cand_eMF->Fit(total_eMF); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit8->SaveAs("Ejets_mucand_landau-2gauss-600-1000_eMFQF" + SQFcut + ".pdf","pdf");
	fit8->SaveAs("Ejets_mucand_landau-2gauss-600-1000_eMFQF" + SQFcut + ".jpeg","jpeg");
    fit8->Close();

	*/
	/*
	//Energy_muon PILEUP landau + 2gauss
	landau_ene = new TF1("landau_ene","landau",0,600);
	gauss1_ene = new TF1("gauss1_ene","gaus",600,900);
	gauss2_ene = new TF1("gauss2_ene","gaus",900,3500);
	total_ene = new TF1("mstotal","landau(0)+gaus(3)+gaus(6)",0,3500);
	landau_eMF = new TF1("landau_eMF","landau",0,600);
	gauss1_eMF = new TF1("gauss1_eMF","gaus",600,900);
	gauss2_eMF = new TF1("gauss2_eMF","gaus",900,3500);
	total_eMF = new TF1("mstotal","landau(0)+gaus(3)+gaus(6)",0,3500);
	landau_ene->SetLineColor(1);
	gauss1_ene->SetLineColor(1);
	gauss2_ene->SetLineColor(1);
	total_ene->SetLineColor(1);
	landau_eMF->SetLineColor(1);
	gauss1_eMF->SetLineColor(1);
	gauss2_eMF->SetLineColor(1);
	total_eMF->SetLineColor(1);

    TCanvas *fit15 = new TCanvas("fit15","fit15",900,700);
    Ejets_pileup_muons_ene->SetTitle("OF2 - Energy of muons candidates pileup signal (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_pileup_muons_ene->Draw(); 
	Ejets_pileup_muons_ene->Fit(landau_ene,"R"); //$
	Ejets_pileup_muons_ene->Fit(gauss1_ene,"R+"); //$
	Ejets_pileup_muons_ene->Fit(gauss2_ene,"R+"); //$
	Double_t par_ene[14];
	landau_ene->GetParameters(&par_ene[0]);
	gauss1_ene->GetParameters(&par_ene[3]);
	gauss2_ene->GetParameters(&par_ene[6]);
	total_ene->SetParameters(par_ene);
	//Ejets_pileup_muons_ene->Fit(total_ene); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit15->SaveAs("Ejets_pileup_mucand_landau-2gauss-600-900_eneQF" + SQFcut + ".pdf","pdf");
	fit15->SaveAs("Ejets_pileup_mucand_landau-2gauss-600-900_eneQF" + SQFcut + ".jpeg","jpeg");
    fit15->Close();

    TCanvas *fit11 = new TCanvas("fit11","fit11",900,700);
    Ejets_pileup_muons_eMF->SetTitle("COF - Energy of muons candidates pileup signal (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_pileup_muons_eMF->Draw(); 
	Ejets_pileup_muons_eMF->Fit(landau_eMF,"R"); //$
	Ejets_pileup_muons_eMF->Fit(gauss1_eMF,"R+"); //$
	Ejets_pileup_muons_eMF->Fit(gauss2_eMF,"R+"); //$
	Double_t par_eMF[14];
	landau_eMF->GetParameters(&par_eMF[0]);
	gauss1_eMF->GetParameters(&par_eMF[3]);
	gauss2_eMF->GetParameters(&par_eMF[6]);
	total_eMF->SetParameters(par_eMF);
	//Ejets_pileup_muons_eMF->Fit(total_eMF); 
	//Comment the last line to run landau+2gauss separated
	//Comment $'s lines to run landau+2gauss together
	fit11->SaveAs("Ejets_pileup_mucand_landau-2gauss-600-900_eMFQF" + SQFcut + ".pdf","pdf");
	fit11->SaveAs("Ejets_pileup_mucand_landau-2gauss-600-900_eMFQF" + SQFcut + ".jpeg","jpeg");
    fit11->Close();


	//Energy_muon PILEUP landau + pol9
	landau2_ene = new TF1("landau2_ene","landau",0,600);
	pol9_ene = new TF1("pol9_ene","pol9",600,3500);
	total2_ene = new TF1("mstotal","landau(0)+pol9(3)",0,3500);
	landau2_eMF = new TF1("g1_eMF","landau",0,600);
	pol9_eMF = new TF1("pol9_eMF","pol9",600,3500);
	total2_eMF = new TF1("mstotal","landau(0)+pol9(3)",0,3500);
	landau2_ene->SetLineColor(1);
    pol9_ene->SetLineColor(1);
	landau2_eMF->SetLineColor(1);
    pol9_eMF->SetLineColor(1);
	total2_ene->SetLineColor(1);
    total2_eMF->SetLineColor(1);

    TCanvas *fit12 = new TCanvas("fit12","fit12",900,700);
    Ejets_pileup_muons_ene->SetTitle("OF2 - Energy of muons candidates pileup signal (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_pileup_muons_ene->Draw(); 
	Ejets_pileup_muons_ene->Fit(landau2_ene,"R"); 
	Ejets_pileup_muons_ene->Fit(pol9_ene,"R+"); 
	fit12->SaveAs("Ejets_pileup_mucand_landau-pol-600_eneQF" + SQFcut + ".pdf","pdf");
	fit12->SaveAs("Ejets_pileup_mucand_landau-pol-600_eneQF" + SQFcut + ".jpeg","jpeg");
    //fit12->Close();

    TCanvas *fit13 = new TCanvas("fit13","fit13",900,700);
    Ejets_pileup_muons_eMF->SetTitle("COF - Energy of muons candidates pileup signal (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_pileup_muons_eMF->Draw(); 
	Ejets_pileup_muons_eMF->Fit(landau2_eMF,"R"); 
	Ejets_pileup_muons_eMF->Fit(pol9_eMF,"R+"); 
	fit13->SaveAs("Ejets_pileup_mucand_landau-pol-600_eMFQF" + SQFcut + ".pdf","pdf");
	fit13->SaveAs("Ejets_pileup_mucand_landau-pol-600_eMFQF" + SQFcut + ".jpeg","jpeg");
    fit13->Close();
	*/


	//ENERGY_MUON LANDAU+POL9
	landau2_ene = new TF1("landaue_ene","landau",0,600);
	pol9_ene = new TF1("pol9_ene","pol9",600,3500);
	total2_ene = new TF1("mstotal_ene","landau(0)+pol9(3)",0,3500);
	landau2_eMF = new TF1("landaue_eMF","landau",0,600);
	pol9_eMF = new TF1("pol9_eMF","pol9",600,3500);
	total2_eMF = new TF1("mstotal_eMF","landau(0)+pol9(3)",0,3500);
	landau2_ene->SetLineColor(1);
    pol9_ene->SetLineColor(1);
	landau2_eMF->SetLineColor(1);
    pol9_eMF->SetLineColor(1);
	total2_ene->SetLineColor(1);
    total2_eMF->SetLineColor(1);
/*
    TCanvas *fit9 = new TCanvas("fit9","fit9",900,700);
    Ejets_mu_cand_ene->SetTitle("OF2 - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_ene->Draw(); 
	Ejets_mu_cand_ene->Fit(landau2_ene,"R"); 
	Ejets_mu_cand_ene->Fit(pol9_ene,"R+"); 
    Ejets_mu_cand_ene->Draw(); 
	fit9->SaveAs("Slide/Ejets_mucand_landau-pol-600_eneQF" + SQFcut + ".pdf","pdf");
	fit9->SaveAs("Slide/Ejets_mucand_landau-pol-600_eneQF" + SQFcut + ".jpeg","jpeg");
    //fit9->Close();


    TCanvas *fit10 = new TCanvas("fit10","fit10",900,700);
    Ejets_mu_cand_eMF->SetTitle("COF - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_eMF->Draw(); 
	Ejets_mu_cand_eMF->Fit(landau2_eMF,"R"); 
	Ejets_mu_cand_eMF->Fit(pol9_eMF,"R+"); 
    Ejets_mu_cand_eMF->Draw(); 
	fit10->SaveAs("Slide/Ejets_mucand_landau-pol-600_eMFQF" + SQFcut + ".pdf","pdf");
	fit10->SaveAs("Slide/Ejets_mucand_landau-pol-600_eMFQF" + SQFcut + ".jpeg","jpeg");
    //fit10->Close();
*/

	/*
	//Ratio E/mm Fit landau
	landau3_ene = new TF1("landau3_ene","landau",4,100);
	landau3_eMF = new TF1("landau3_eMF","landau",4,100);
	landau3_ene->SetLineColor(1);
	landau3_eMF->SetLineColor(1);

	TCanvas *fit16 = new TCanvas("fit16","fit16",900,700);
	Ratio_E_l_ene->SetTitle("OF2 - Energy/length in each layer (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
	Ratio_E_l_ene->Draw(); 
	Ratio_E_l_ene->Fit(landau3_ene);
	fit16->SaveAs("Ratio_E_l-landaufit-eneQF" + SQFcut + ".pdf","pdf");
	fit16->SaveAs("Ratio_E_l-landaufit-eneQF" + SQFcut + ".jpeg","jpeg");
    fit16->Close();

	TCanvas *fit17 = new TCanvas("fit17","fit17",900,700);
	Ratio_E_l_eMF->SetTitle("COF - Energy/length in each layer (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
	Ratio_E_l_eMF->Draw(); 
	Ratio_E_l_eMF->Fit(landau3_eMF);
	fit17->SaveAs("Ratio_E_l-landaufit-eMFQF" + SQFcut + ".pdf","pdf");
	fit17->SaveAs("Ratio_E_l-landaufit-eMFQF" + SQFcut + ".jpeg","jpeg");
    fit17->Close();
	//--------------------




    TCanvas *fit6 = new TCanvas("fit6","fit6",900,700);
    Ejets_mu_cand_eMF->SetTitle("COF - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_eMF->Draw(); 
	Ejets_mu_cand_eMF->Fit(g1_eMF);
	Ejets_mu_cand_eMF->Fit(g2_eMF);
	//Ejets_mu_cand_eMF->Fit(g3_eMF);
	g1_eMF->GetParameters(&par[0]);
	g2_eMF->GetParameters(&par[3]);
	//g3_eMF->GetParameters(&par[6]);
	total_eMF->SetParameters(par);
	Ejets_mu_cand_eMF->Fit(total_eMF);


    fit6->SaveAs("QF" + SQFcut + "/Ejets_mu_cand_eMFQF" + SQFcut + ".pdf","pdf");
    TFile *fit1 = new TFile("QF" + SQFcut + "/Ejets_mu_cand_eMFQF" + SQFcut + ".root","RECREATE");
    Ejets_mu_cand_eMF->Write();
    fit6->Write();
    //fit6->Close();
    //fit1->Close();



    TCanvas *fit7 = new TCanvas("fit7","fit7",900,700);
    Ejets_mu_cand_ene->SetTitle("OF2 - Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_ene->Draw(); 
	Ejets_mu_cand_ene->Fit(g1_ene,"R");
	Ejets_mu_cand_ene->Fit(g2_ene,"R+");
	//Ejets_mu_cand_ene->Fit(g3_ene);
	g1_ene->GetParameters(&par[0]);
	g2_ene->GetParameters(&par[3]);
	//g3_ene->GetParameters(&par[6]);
	total_ene->SetParameters(par);
	//Ejets_mu_cand_ene->Fit(total_ene,"R");

    fit7->SaveAs("QF" + SQFcut + "/Ejets_mu_cand_eneQF" + SQFcut + ".pdf","pdf");
    TFile *fit2 = new TFile("QF" + SQFcut + "/Ejets_mu_cand_eneQF" + SQFcut + ".root","RECREATE");
    Ejets_mu_cand_ene->Write();
    fit7->Write();
    //fit7->Close();
    //fit2->Close();


    TCanvas *p12 = new TCanvas("p12","p12",900,700);
    EneXQF_30->Draw(); 
    //p1->SetLogy()
    p12->SaveAs("TProfiles/EneXQF_30.pdf","pdf");
    p12->Close();

    TCanvas *p11 = new TCanvas("p11","p11",900,700);
    EMFXQF_30->Draw(); 
    //p5->SetLogy();
    p11->SaveAs("TProfiles/EMFXQF_30.pdf","pdf");
    p11->Close();

    TCanvas *p9 = new TCanvas("p9","p9",900,700);
    EneXQF_50->Draw(); 
    //p1->SetLogy();
    p9->SaveAs("TProfiles/EneXQF_50.pdf","pdf");
    p9->Close();

    TCanvas *p10 = new TCanvas("p10","p10",900,700);
    EMFXQF_50->Draw(); 
    //p5->SetLogy();
    p10->SaveAs("TProfiles/EMFXQF_50.pdf","pdf");
    p10->Close();

    TCanvas *p1 = new TCanvas("p1","p1",900,700);
    EneXQF_100->Draw(); 
    //p1->SetLogy();
    p1->SaveAs("TProfiles/EneXQF_100.pdf","pdf");
    p1->Close();

    TCanvas *p5 = new TCanvas("p5","p5",900,700);
    EMFXQF_100->Draw(); 
    //p5->SetLogy();
    p5->SaveAs("TProfiles/EMFXQF_100.pdf","pdf");
    p5->Close();

    TCanvas *b1 = new TCanvas("b1","b1",900,700);
    EnexQF->Draw(); 
    //e1->SetLogy();
    b1->SaveAs("TProfiles/EnexQF.pdf","pdf");
    b1->Close();

    TCanvas *b2 = new TCanvas("b2","b2",900,700);
    EMFxQF->Draw(); 
    //e1->SetLogy();
    b2->SaveAs("TProfiles/EMFxQF.pdf","pdf");
    b2->Close();

    TCanvas *b3 = new TCanvas("b3","b3",900,700);
    EnexQF2->Draw(); 
    //e1->SetLogy();
    b3->SaveAs("TProfiles/EnexQF2.pdf","pdf");
    b3->Close();

    TCanvas *b4 = new TCanvas("b4","b4",900,700);
    EMFxQF2->Draw(); 
    //e1->SetLogy();
    b4->SaveAs("TProfiles/EMFxQF2.pdf","pdf");
    b4->Close();

    TCanvas *b5 = new TCanvas("b5","b5",900,700);
    EnexQF3->Draw(); 
    //e1->SetLogy();
    b5->SaveAs("TProfiles/EnexQF3.pdf","pdf");
    b5->Close();

    TCanvas *b6 = new TCanvas("b6","b6",900,700);
    EMFxQF3->Draw(); 
    //e1->SetLogy();
    b6->SaveAs("TProfiles/EMFxQF3.pdf","pdf");
    b6->Close();
	*/

	/*
    TCanvas *t1 = new TCanvas("t1","t1",900,700);
    Emu_layer1_ene->SetTitle("Energy of muons in the first layer QF<"+SQFcut);
    Emu_layer1_ene->Draw();
    Emu_layer1_eMF->Draw("SAME");
    leg->Draw("SAME");
    t1->SaveAs("Emu_1layerQF" + SQFcut + ".pdf","pdf");
    t1->Close();

    TCanvas *t2 = new TCanvas("t2","t2",900,700);
    Emu_layer2_eMF->SetTitle("Energy of muons in the second layer QF<"+SQFcut);
    Emu_layer2_eMF->Draw();
    Emu_layer2_ene->Draw("SAME");
    leg->Draw("SAME");
    t2->SaveAs("Emu_2layerQF" + SQFcut + ".pdf","pdf");
    t2->Close();

    TCanvas *t3 = new TCanvas("t3","t3",900,700);
    Emu_layer3_ene->SetTitle("Energy of muons in the third layer QF<"+SQFcut);
    Emu_layer3_ene->Draw();
    Emu_layer3_eMF->Draw("SAME");
    leg->Draw("SAME");
    t3->SaveAs("Emu_3layerQF" + SQFcut + ".pdf","pdf");
    t3->Close();
	*/

	/*
    TCanvas *lll4 = new TCanvas("lll4","lll4",900,700);
    lll4->SetLogy();
    Elayersall_ene->SetTitle("Jets with energy in all layers (Zr=1, QF<"+SQFcut+")");
    Elayersall_ene->Draw();
    Elayersall_eMF->Draw("SAME");
    leg->Draw("SAME");
    lll4->SaveAs("Elayersall_QF" + SQFcut + ".pdf","pdf");
    lll4->Close();

    TCanvas *l1 = new TCanvas("l1","l1",900,700);
    l1->SetLogx();
    Elayers12_ene->Draw("LEGO");
    l1->SaveAs("Elayers12ene_QF" + SQFcut + ".pdf","pdf");
    l1->Close();

    TCanvas *ll1 = new TCanvas("ll1","ll1",900,700);
    ll1->SetLogx();
    Elayers12_eMF->Draw("LEGO");
    ll1->SaveAs("Elayers12eMF_QF" + SQFcut + ".pdf","pdf");
    ll1->Close();

    TCanvas *l2 = new TCanvas("l2","l2",900,700);
    l2->SetLogx();
    Elayers13_ene->Draw("LEGO");
    l2->SaveAs("Elayers13ene_QF" + SQFcut + ".pdf","pdf");
    l2->Close();

    TCanvas *ll2 = new TCanvas("ll2","ll2",900,700);
    ll2->SetLogx();
    Elayers13_eMF->Draw("LEGO");
    ll2->SaveAs("Elayers13eMF_QF" + SQFcut + ".pdf","pdf");
    ll2->Close();

    TCanvas *l3 = new TCanvas("l3","l3",900,700);
    l3->SetLogx();
    Elayers23_ene->Draw("LEGO");
    l3->SaveAs("Elayers23ene_QF" + SQFcut + ".pdf","pdf");
    l3->Close();

    TCanvas *ll3 = new TCanvas("ll3","ll3",900,700);
    ll3->SetLogx();
    Elayers23_eMF->Draw("LEGO");
    ll3->SaveAs("Elayers23eMF_QF" + SQFcut + ".pdf","pdf");
    ll3->Close();

    TCanvas *l4 = new TCanvas("l4","l4",900,700);
    l4->SetLogx();
    Elayers1_ene->Draw("LEGO");
    l4->SaveAs("Elayers1ene_QF" + SQFcut + ".pdf","pdf");
    l4->Close();

    TCanvas *ll4 = new TCanvas("ll4","ll4",900,700);
    ll4->SetLogx();
    Elayers1_eMF->Draw("LEGO");
    ll4->SaveAs("Elayers1eMF_QF" + SQFcut + ".pdf","pdf");
    ll4->Close();

    TCanvas *l5 = new TCanvas("l5","l5",900,700);
    l5->SetLogx();
    Elayers2_ene->Draw("LEGO");
    l5->SaveAs("Elayers2ene_QF" + SQFcut + ".pdf","pdf");
    l5->Close();

    TCanvas *ll5 = new TCanvas("ll5","ll5",900,700);
    ll5->SetLogx();
    Elayers2_eMF->Draw("LEGO");
    ll5->SaveAs("Elayers2eMF_QF" + SQFcut + ".pdf","pdf");
    ll5->Close();

    TCanvas *l6 = new TCanvas("l6","l6",900,700);
    l6->SetLogx();
    Elayers3_ene->Draw("LEGO");
    l6->SaveAs("Elayers3ene_QF" + SQFcut + ".pdf","pdf");
    l6->Close();

    TCanvas *ll6 = new TCanvas("ll6","ll6",900,700);
    ll6->SetLogx();
    Elayers3_eMF->Draw("LEGO");
    ll6->SaveAs("Elayers3eMF_QF" + SQFcut + ".pdf","pdf");
    ll6->Close();
	*/
	/*
    TCanvas *h1 = new TCanvas("h1","h1",900,700);
    //h1->SetLogy();
    Ejets_all3layers_eMF->SetTitle("Jets all layers (Zr=1, QF<"+SQFcut+")");
    Ejets_all3layers_eMF->Draw();
    Ejets_all3layers_ene->Draw("SAME");
    leg->Draw("SAME");
    h1->SaveAs("E_single_jet_all_layers(2)QF" + SQFcut + ".pdf","pdf");
    h1->Close();

    TCanvas *h3 = new TCanvas("h3","h3",900,700);
    //h3->SetLogy();
    Ejets_layers_eMF->SetTitle("Single cell jet all TileCal layers (Zr=1, QF<"+SQFcut+")");
    Ejets_layers_eMF->Draw();
    Ejets_layers_ene->Draw("SAME");
    leg->Draw("SAME");
    h3->SaveAs("E_single_jet_all_layersQF" + SQFcut + ".pdf","pdf");
    h3->Close();
	
    TCanvas *h2 = new TCanvas("h2","h2",900,700);
    //h2->SetLogy();
    Ejets_1layers_eMF->SetTitle("Jets with energy only in the first layer (Zr=1 or Zr=1+1cell, QF<"+SQFcut+")");
    Ejets_1layers_eMF->Draw();
    Ejets_1layers_ene->Draw("SAME");
    leg->Draw("SAME");
    h2->SaveAs("E_single_jet_just_layer1QF" + SQFcut + "NEW.pdf","pdf");
    h2->Close();

    TCanvas *h4 = new TCanvas("h4","h4",900,700);	
    //h4->SetLogy();
    Ejets_2layers_ene->SetTitle("Jets with energy only in the second layer (Zr=1 or Zr=1+1cell, QF<"+SQFcut+")");
    Ejets_2layers_ene->Draw();
    Ejets_2layers_eMF->Draw("SAME");
    leg->Draw("SAME");
    h4->SaveAs("E_single_jet_just_layer2QF" + SQFcut + ".pdf","pdf");
    h4->Close();

    TCanvas *h5 = new TCanvas("h4","h4",900,700);	
    //h5->SetLogy();
    Ejets_3layers_eMF->SetTitle("Jets with energy only in the third layer (Zr=1 or Zr=1+1cell, QF<"+SQFcut+")");
    Ejets_3layers_eMF->Draw();
    Ejets_3layers_ene->Draw("SAME");
    leg->Draw("SAME");
    h5->SaveAs("E_single_jet_just_layer3QF" + SQFcut + ".pdf","pdf");
    h5->Close();

    TCanvas *h21 = new TCanvas("h21","h21",900,700);
    //h21->SetLogy();
    Ejets_12layers_eMF->SetTitle("Jets with energy only in the first and second layer (Zr=1 or Zr=1+1cell, QF<"+SQFcut+")");
    Ejets_12layers_eMF->Draw();
    Ejets_12layers_ene->Draw("SAME");
    leg->Draw("SAME");
    h21->SaveAs("E_single_jet_just_layer12QF" + SQFcut + "NEW.pdf","pdf");
    h21->Close();

    TCanvas *h42 = new TCanvas("h42","h42",900,700);	
    //h42->SetLogy();
    Ejets_23layers_ene->SetTitle("Jets with energy only in the second and third layer (Zr=1 or Zr=1+1cell, QF<"+SQFcut+")");
    Ejets_23layers_ene->Draw();
    Ejets_23layers_eMF->Draw("SAME");
    leg->Draw("SAME");
    h42->SaveAs("E_single_jet_just_layer23QF" + SQFcut + ".pdf","pdf");
    h42->Close();


	
    TCanvas *e1 = new TCanvas("e1","e1",900,700);
    PT_miss_eMF2->SetTitle("pT miss (QF<" + SQFcut + ")");
    PT_miss_eMF2->Draw(); 
    PT_miss_ene2->Draw("SAME");
    leg->Draw("SAME");
    //e1->SetLogy();
    e1->SaveAs("pT_miss_ene_eMF2QF" + SQFcut + ".pdf","pdf");
    TFile *eb1 = new TFile("pT_miss2_QF" + SQFcut + ".root","RECREATE");
    PT_miss_eMF2->Write();
    PT_miss_ene2->Write();
    eb1->Write();
    eb1->Close();
    e1->Close();
	*/
    //Create the legend
    leg2 = new TLegend(0.7,0.4,0.8,0.65);
    //leg = new TLegend(0.1,0.75,0.3,0.9);
    leg2->SetBorderSize(0);
    //leg->SetHeader("Legend");
    leg2->AddEntry(Ejets_eMF,"COF","l");
    leg2->AddEntry(Ejets_ene,"OF2","l");
    leg2->SetTextSize(0.035);
	/*
    TCanvas *e5 = new TCanvas("e5","e5",900,700);
    Phi_miss_eMF2->SetTitle("phi miss (QF<" + SQFcut + ")");
    Phi_miss_eMF2->Draw(); 
    Phi_miss_ene2->Draw("SAME");
    leg7->Draw("SAME");
    //e5->SetLogy();
    e5->SaveAs("phi_miss_ene_eMF2QF" + SQFcut + ".pdf","pdf");
    e5->SaveAs("phi_miss_ene_eMF2QF" + SQFcut + ".jpeg","jpeg");
    TFile *eb5 = new TFile("phi_miss2_QF" + SQFcut + ".root","RECREATE");
    Phi_miss_eMF2->Write();
    Phi_miss_ene2->Write();
    eb5->Write();
    eb5->Close();
    e5->Close();
	*/

	/*
    TCanvas *g6 = new TCanvas("g6","g6",900,700);
    //g3->SetLogy();
    E123_R1_ene->Draw("LEGO"); 
    g6->SaveAs("E123_R1_eneQF" + SQFcut + ".pdf","pdf");
    TFile *gb6 = new TFile("E123_R1_eneQF" + SQFcut + ".root","RECREATE");
    E123_R1_ene->Write();
    gb6->Write();
    gb6->Close();
    g6->Close();

    TCanvas *g7 = new TCanvas("g7","g7",900,700);
    //g3->SetLogy();
    E123_R1_eMF->Draw("LEGO"); 
    g7->SaveAs("E123_R1_eMFQF" + SQFcut + ".pdf","pdf");
    TFile *gb7 = new TFile("E123_R1_eMFQF" + SQFcut + ".root","RECREATE");
    E123_R1_eMF->Write();
    gb7->Write();
    gb7->Close();
    g7->Close();

    TCanvas *g8 = new TCanvas("g8","g8",900,700);
    //g3->SetLogy();
    E123_R2_ene->Draw("LEGO"); 
    g8->SaveAs("E123_R2_eneQF" + SQFcut + ".pdf","pdf");
    TFile *gb8 = new TFile("E123_R2_eneQF" + SQFcut + ".root","RECREATE");
    E123_R2_ene->Write();
    gb8->Write();
    gb8->Close();
    g8->Close();

    TCanvas *g9 = new TCanvas("g9","g9",900,700);
    //g3->SetLogy();
    E123_R2_eMF->Draw("LEGO"); 
    g9->SaveAs("E123_R2_eMFQF" + SQFcut + ".pdf","pdf");
    TFile *gb9 = new TFile("E123_R2_eMFQF" + SQFcut + ".root","RECREATE");
    E123_R2_eMF->Write();
    gb9->Write();
    gb9->Close();
    g9->Close();

    TCanvas *g10 = new TCanvas("g10","g10",900,700);
    //g3->SetLogy();
    E123_R3_ene->Draw("LEGO"); 
    g10->SaveAs("E123_R3_eneQF" + SQFcut + ".pdf","pdf");
    TFile *gb10 = new TFile("E123_R3_eneQF" + SQFcut + ".root","RECREATE");
    E123_R3_ene->Write();
    gb10->Write();
    gb10->Close();
    g10->Close();

    TCanvas *g11 = new TCanvas("g11","g11",900,700);
    //g3->SetLogy();
    E123_R3_eMF->Draw("LEGO"); 
    g11->SaveAs("E123_R3_eMFQF" + SQFcut + ".pdf","pdf");
    TFile *gb11 = new TFile("E123_R3_eMFQF" + SQFcut + ".root","RECREATE");
    E123_R3_eMF->Write();
    gb11->Write();
    gb11->Close();
    g11->Close();

	*/
	//TFile *pTmissCOF = new TFile("pTmissCOF.root","NEW");
	//PT_miss_eMF->Write("pT_miss_eMF"+SQFcut,TObject::kWriteDelete);
	//PT_miss_eMF->Write();
	//pTmissCOF->Close();

	//TFile *pTmissOF2 = new TFile("pTmissOF2.root","NEW");
	//PT_miss_ene->Write("pT_miss_ene"+SQFcut,TObject::kWriteDelete);
	//PT_miss_ene->Write();
	//pTmissOF2->Close();

	//TFile *phimissCOF = new TFile("phimissCOF.root","NEW");
	//Phi_miss_eMF->Write("phi_miss_eMF"+SQFcut,TObject::kWriteDelete);
	//phi_miss_eMF->Write();
	//phimissCOF->Close();

	//TFile *phimissOF2 = new TFile("phimissOF2.root","NEW");
	//Phi_miss_ene->Write("phi_miss_ene"+SQFcut,TObject::kWriteDelete);
	//phi_miss_ene->Write();
	//phimissOF2->Close();

/*

	//Everything here is alright
	
    //------------------------------------------JETS----------------------------------------

	//Energy of jets
	//All jets
    TCanvas *c1 = new TCanvas("c1","c1",900,700);
    Ejets_eMF->SetTitle("Energy of the jets of all events (QF<" + SQFcut + ")");
    Ejets_eMF->Draw(); 
    Ejets_ene->Draw("SAME");
    leg->Draw("SAME");
    c1->SetLogy();
    if (QFcut == 50) c1->SaveAs("QF" + SQFcut + "/Jets/Ejets_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c1->SaveAs("QF" + SQFcut + "/Jets/Ejets_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fb1 = new TFile("QF" + SQFcut + "/Jets/Ejets_QF" + SQFcut + ".root","RECREATE");
    Ejets_eMF->Write();
    Ejets_ene->Write();
    Ejets_eMF->SetTitle("Jet Energy - All events");
    fb1->Write();
    fb1->Close();
    c1->Close();

	//Zr>1
    TCanvas *c2 = new TCanvas("c2","c2",900,700);
    Ejets_z1_eMF->SetTitle("Energy of the jets (Zr>1, QF< " + SQFcut + ")");
    Ejets_z1_eMF->Draw(); 
    Ejets_z1_ene->Draw("SAME");
    leg->Draw("SAME");
    c2->SetLogy();
    if (QFcut == 50) c2->SaveAs("QF" + SQFcut + "/Jets/Ejets_z1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c2->SaveAs("QF" + SQFcut + "/Jets/Ejets_z1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fb2 = new TFile("QF" + SQFcut + "/Jets/Ejets_z1_QF" + SQFcut + ".root","RECREATE");
    Ejets_z1_eMF->Write();
    Ejets_z1_ene->Write();
    fb2->Write();
    fb2->Close();
    c2->Close();

	//Zr>2
    TCanvas *c3 = new TCanvas("c3","c3",900,700);
    Ejets_z2_eMF->SetTitle("Energy of the jets (Zr>2, QF<" + SQFcut + ")");
    Ejets_z2_eMF->Draw(); 
    Ejets_z2_ene->Draw("SAME");
    leg->Draw("SAME");
    c3->SetLogy();
    if (QFcut == 50) c3->SaveAs("QF" + SQFcut + "/Jets/Ejets_z2_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c3->SaveAs("QF" + SQFcut + "/Jets/Ejets_z2_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fb3 = new TFile("QF" + SQFcut + "/Jets/Ejets_z2_QF" + SQFcut + ".root","RECREATE");
    Ejets_z2_eMF->Write();
    Ejets_z2_ene->Write();
    fb3->Write();
    fb3->Close();
    c3->Close();

	//Energy of Zr=1 jets
    TCanvas *c29 = new TCanvas("c29","c29",900,700);
    Ejets_zequal1_eMF->SetTitle("Energy of the jets (Zr=1, QF< " + SQFcut + ")");
    Ejets_zequal1_eMF->Draw(); 
    Ejets_zequal1_ene->Draw("SAME");
    leg->Draw("SAME");
    c29->SetLogy();
    if (QFcut == 50) c29->SaveAs("QF" + SQFcut + "/Ejets_zequal1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c29->SaveAs("QF" + SQFcut + "/Ejets_zequal1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fb29 = new TFile("QF" + SQFcut + "/Ejets_zequal1_QF" + SQFcut + ".root","RECREATE");
    Ejets_zequal1_eMF->Write();
    Ejets_zequal1_ene->Write();
    fb29->Write();
    fb29->Close();
    //c29->Close();


	//Number of jets per event
	//All jets
    TCanvas *c8 = new TCanvas("c8","c8",900,700);
    Njets_per_event_ene->SetTitle("Number of jets per event (QF<" + SQFcut + ")");
    Njets_per_event_ene->Draw(); 
    Njets_per_event_eMF->Draw("SAME");
    leg->Draw("SAME");
    if (QFcut == 50) c8->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c8->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    c8->Close();
    TCanvas *c15 = new TCanvas("c15","c15",900,700);
    Njets_per_event_ene_eMF->Draw();
    if (QFcut == 50) c15->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_ene_versus_eMF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c15->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_ene_versus_eMF" + SQFcut + ".jpeg","jpeg");
    c15->Close();
    TFile *fb8 = new TFile("QF" + SQFcut + "/Jets/Njets_per_event_QF" + SQFcut + ".root","RECREATE");
    Njets_per_event_ene_eMF->Write();
    Njets_per_event_eMF->Write();
    Njets_per_event_ene->Write();
    fb8->Write();
    fb8->Close();

	//Zr>1
    TCanvas *c9 = new TCanvas("c9","c9",900,700);
    Njets_per_event_z1_ene->SetTitle("Number of jets per event (Zr>1, QF<" + SQFcut + ")");
    Njets_per_event_z1_ene->Draw(); 
    Njets_per_event_z1_eMF->Draw("SAME");
    leg->Draw("SAME");
    if (QFcut == 50) c9->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c9->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    c9->Close();
    TCanvas *c16 = new TCanvas("c16","c16",900,700);
    Njets_per_event_z1_ene_eMF->Draw();
    if (QFcut == 50) c16->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z1_ene_versus_eMF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c16->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z1_ene_versus_eMF" + SQFcut + ".jpeg","jpeg");
    c16->Close();
    TFile *fb9 = new TFile("QF" + SQFcut + "/Jets/Njets_per_event_z1_QF" + SQFcut + ".root","RECREATE");
    Njets_per_event_z1_eMF->Write();
    Njets_per_event_z1_ene_eMF->Write();
    Njets_per_event_z1_ene->Write();
    fb9->Write();
    fb9->Close();

	//Zr>2
    TCanvas *c10 = new TCanvas("c10","c10",900,700);
    Njets_per_event_z2_ene->SetTitle("Number of jets per event (Zr>2, QF<" + SQFcut + ")");
    Njets_per_event_z2_ene->Draw(); 
    Njets_per_event_z2_eMF->Draw("SAME");
    leg->Draw("SAME");
    if (QFcut == 50) c10->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z2_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c10->SaveAs("QF" + SQFcut + "Jets//Njets_per_event_z2_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    c10->Close();
    TCanvas *c17 = new TCanvas("c17","c17",900,700);
    Njets_per_event_z2_ene_eMF->Draw();
    if (QFcut == 50) c17->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z2_ene_versus_eMF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c17->SaveAs("QF" + SQFcut + "/Jets/Njets_per_event_z2_ene_versus_eMF" + SQFcut + ".jpeg","jpeg");
    c17->Close();
    TFile *fb10 = new TFile("QF" + SQFcut + "Jets//Njets_per_event_z2_QF" + SQFcut + ".root","RECREATE");
    Njets_per_event_z2_eMF->Write();
    Njets_per_event_z2_ene_eMF->Write();
    Njets_per_event_z2_ene->Write();
    fb10->Write();
    fb10->Close();	

	//Number of Zr=1 jets
    TCanvas *c99 = new TCanvas("c99","c99",900,700);
    Njets_per_event_zequal1_ene->SetTitle("Number of jets per event (Zr=1, QF<" + SQFcut + ")");
    Njets_per_event_zequal1_ene->Draw(); 
    Njets_per_event_zequal1_eMF->Draw("SAME");
    leg->Draw("SAME");
    if (QFcut == 50) c99->SaveAs("QF" + SQFcut + "/Njets_per_event_zequal1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c99->SaveAs("QF" + SQFcut + "/Njets_per_event_zequal1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    //c99->Close();
    TCanvas *c169 = new TCanvas("c169","c169",900,700);
    Njets_per_event_zequal1_ene_eMF->Draw();
    if (QFcut == 50) c169->SaveAs("QF" + SQFcut + "/Njets_per_event_zequal1_ene_versus_eMF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c169->SaveAs("QF" + SQFcut + "/Njets_per_event_zequal1_ene_versus_eMF" + SQFcut + ".jpeg","jpeg");
    //c169->Close();
    TFile *fb99 = new TFile("QF" + SQFcut + "/Njets_per_event_zequal1_QF" + SQFcut + ".root","RECREATE");
    Njets_per_event_zequal1_eMF->Write();
    Njets_per_event_zequal1_ene_eMF->Write();
    Njets_per_event_zequal1_ene->Write();
    fb99->Write();
    fb99->Close();

   	//DeltaR distribuition
    TCanvas *z2 = new TCanvas("z2","z2",900,700);
    z2->SetLogy();
    DeltaR_dist_eMF->SetTitle("Delta R distribution QF<"+SQFcut);
    DeltaR_dist_eMF->Draw();
    DeltaR_dist_ene->Draw("SAME");
    leg->Draw("SAME");
    if (QFcut == 50) z2->SaveAs("QF" + SQFcut + "/DeltaR_distQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) z2->SaveAs("QF" + SQFcut + "/DeltaR_distQF" + SQFcut + ".jpeg","jpeg");
    // z2->Close();

	//Zr Distribuition
    TCanvas *z1 = new TCanvas("z1","z1",900,700);
    z1->SetLogy();
    Zr_dist_eMF->SetTitle("Zr distribuition QF<"+SQFcut);
    Zr_dist_eMF->Draw();
    Zr_dist_ene->Draw("SAME");
    leg->Draw("SAME");
    z1->SaveAs("Zr_distQF" + SQFcut + ".pdf","pdf");
    //z1->Close();

	//Ratio of Energy/mm
    TCanvas *g3 = new TCanvas("g3","g3",900,700);
    Ratio_E_l_eMF->SetTitle("Energy/length of each jet in each layer (Zr = 1 or 1+1cell, QF<" + SQFcut + ")");
    Ratio_E_l_eMF->Draw(); 
    Ratio_E_l_ene->Draw("SAME");
    leg->Draw("SAME");
    //g3->SetLogy();
    if (QFcut == 50) g3->SaveAs("QF" + SQFcut + "/Ratio_E_l_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) g3->SaveAs("QF" + SQFcut + "/Ratio_E_l_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *gb3 = new TFile("QF" + SQFcut + "/Ratio_E_l_QF" + SQFcut + ".root","RECREATE");
    Ratio_E_l_ene->Write();
    Ratio_E_l_eMF->Write();
    gb3->Write();
    gb3->Close();
    // g3->Close();



	//-----------------------------------------E miss--------------------------------------------------------------

	//pTmiss events pile_ev
    TCanvas *e46 = new TCanvas("e46","e46",900,700);
    PT_miss_pileev_eMF2->SetTitle("pT miss");
    PT_miss_pileev_eMF2->Draw(); 
    //PT_miss_ene2->SetLineColor(2);
    PT_miss_pileev_ene2->Draw("SAME");
    leg->Draw("SAME");
    //e46->SetLogy();
    e46->SaveAs("QF" + SQFcut + "/E_miss/pT_miss_pilev_ene_eMFQF" + SQFcut + ".pdf","pdf");
    TFile *eb46 = new TFile("QF" + SQFcut + "/E_miss/pT_miss_pilev_QF" + SQFcut + ".root","RECREATE");
    PT_miss_pileev_eMF2->Write();
    PT_miss_pileev_ene2->Write();
    PT_miss_pileev_eMF2->SetTitle("pT miss");
    e46->SaveAs("QF" + SQFcut + "/E_miss/pT_miss_pilev_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    eb46->Write();
    eb46->Close();
    //e46->Close();
*/

	//pTmiss 
    TCanvas *z109 = new TCanvas("z109","z109",900,700);
    PT_miss_ene_eMF->SetTitle("pT miss (QF_OF2<" + SQFcut_OF2 + ")");
    PT_miss_ene_eMF->Draw();
    //z109->SaveAs("QF" + SQFcut + "/E_miss/pT_miss_ene_eMF_2D" + SQFcut + ".pdf","pdf");
    //z109->SaveAs("QF" + SQFcut + "/E_miss/pT_miss_ene_eMF_2D" + SQFcut + ".jpeg","jpeg");
    z109->SaveAs("QFanalysis/pT_miss_ene_eMF_2D" + SQFcut_OF2 + ".pdf","pdf");
    z109->SaveAs("QFanalysis/pT_miss_ene_eMF_2D" + SQFcut_OF2 + ".jpeg","jpeg");
    //z109->Close();

    TCanvas *e4 = new TCanvas("e4","e4",900,700);
    PT_miss_eMF->SetTitle("pT miss (QF_OF2<" + SQFcut_OF2 + ")");
    PT_miss_eMF->Draw(); 
    PT_miss_ene->Draw("SAME");
    leg->Draw("SAME");
    //e4->SetLogy();
    //e4->SaveAs("QF" + SQFcut + "/pT_miss_ene_eMFQF" + SQFcut + ".pdf","pdf");
    e4->SaveAs("QFanalysis/pT_miss_ene_eMFQF" + SQFcut_OF2 + ".pdf","pdf");
    e4->SaveAs("QFanalysis/pT_miss_ene_eMFQF" + SQFcut_OF2 + ".jpeg","jpeg");
    TFile *eb4 = new TFile("QFanalysis/pT_miss_QF" + SQFcut_OF2 + ".root","RECREATE");
    PT_miss_eMF->Write();
    PT_miss_ene->Write();
    //e4->SaveAs("QF" + SQFcut + "/E_miss/pT_miss_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    e4->SaveAs("QFanalysis/pT_miss_ene_eMFQF" + SQFcut_OF2 + ".jpeg","jpeg");
    eb4->Write();
    eb4->Close();
    //e4->Close();
/*
	//phi_miss
    TCanvas *z10 = new TCanvas("z10","z10",900,700);
    Phi_miss_ene_eMF->Draw();
    z10->SaveAs("QF" + SQFcut + "/E_miss/phi_miss_ene_eMF_2D" + SQFcut + ".pdf","pdf");
    z10->SaveAs("QF" + SQFcut + "/E_miss/phi_miss_ene_eMF_2D" + SQFcut + ".jpeg","jpeg");
    //z10->Close();

    TCanvas *e6 = new TCanvas("e6","e6",900,700);
    Phi_miss_ene->SetTitle("phi miss (QF<" + SQFcut + ")");
    Phi_miss_ene->Draw(); 
    Phi_miss_eMF->Draw("SAME");
    leg7->Draw("SAME");
    //e6->SetLogy();
    e6->SaveAs("QF" + SQFcut + "/E_miss/phi_miss_ene_eMFQF" + SQFcut + ".pdf","pdf");
    e6->SaveAs("QF" + SQFcut + "/E_miss/phi_miss_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *eb6 = new TFile("QF" + SQFcut + "/E_miss/phi_miss_QF" + SQFcut + ".root","RECREATE");
    Phi_miss_eMF->Write();
    Phi_miss_ene->Write();
    eb6->Write();
    eb6->Close();
    //e6->Close();

	//phi_cm
    TCanvas *e3 = new TCanvas("e3","e3",900,700);
    phi_cm_ene->SetTitle("PHI_CM per jet (QF<" + SQFcut + ")");
    phi_cm_ene->Draw(); 
    phi_cm_eMF->Draw("SAME");
    leg9->Draw("SAME");
    //e3->SetLogy();
    e3->SaveAs("QF" + SQFcut + "/E_miss/phi_cm_ene_eMFQF" + SQFcut + ".pdf","pdf");
    e3->SaveAs("QF" + SQFcut + "/E_miss/phi_cm_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *eb3 = new TFile("QF" + SQFcut + "/E_miss/phi_cm_QF" + SQFcut + ".root","RECREATE");
    phi_cm_eMF->Write();
    phi_cm_ene->Write();
    eb3->Write();
    eb3->Close();
    //e3->Close();

	//eta_cm
    TCanvas *e2 = new TCanvas("e2","e2",900,700);
    eta_cm_eMF->SetTitle("ETA_CM per jet (QF<" + SQFcut + ")");
    eta_cm_eMF->Draw(); 
    eta_cm_ene->Draw("SAME");
    leg8->Draw("SAME");
    //e2->SetLogy();
    e2->SaveAs("QF" + SQFcut + "/E_miss/eta_cm_ene_eMFQF" + SQFcut + ".pdf","pdf");
    e2->SaveAs("QF" + SQFcut + "/E_miss/eta_cm_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *eb2 = new TFile("QF" + SQFcut + "/E_miss/eta_cm_QF" + SQFcut + ".root","RECREATE");
    eta_cm_eMF->Write();
    eta_cm_ene->Write();
    eb2->Write();
    eb2->Close();
    //e2->Close();


	//--------------------------------------Energy of muons candidates---------------------------------------------

	//Zr = 1 or 1 + 1 cell
    TCanvas *j6 = new TCanvas("j6","j6",900,700);
    Ejets_mu_cand_eMF->SetTitle("Energy of muons candidates (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand_eMF->Draw(); 
    Ejets_mu_cand_ene->Draw("SAME");
    leg->Draw("SAME");
    j6->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j6->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    j6->SetLogy();
    j6->SaveAs("QF" + SQFcut + "/Muons/QF" + SQFcut + "/Ejets_mu_cand_log_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j6->SaveAs("QF" + SQFcut + "/Muons/QF" + SQFcut + "/Ejets_mu_cand_log_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fj6 = new TFile("QF" + SQFcut + "/Muons/Ejets_mu_cand_QF" + SQFcut + ".root","RECREATE");
    Ejets_mu_cand_eMF->Write();
    Ejets_mu_cand_ene->Write();
    fj6->Write();
    fj6->Close();
    // j6->Close();

	//Zr = 1
    TCanvas *j2 = new TCanvas("j2","j2",900,700);
    Ejets_mu_cand1_eMF->SetTitle("Energy of muons candidates (Zr = 1, QF<" + SQFcut + ")");
    Ejets_mu_cand1_eMF->Draw(); 
    Ejets_mu_cand1_ene->Draw("SAME");
    leg->Draw("SAME");
    j2->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j2->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    //j2->SetLogy();
    j2->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand1_log_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j2->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand1_log_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fj2 = new TFile("QF" + SQFcut + "/Muons/Ejets_mu_cand1_QF" + SQFcut + ".root","RECREATE");
    Ejets_mu_cand1_eMF->Write();
    Ejets_mu_cand1_ene->Write();
    fj2->Write();
    fj2->Close();
    //j2->Close();

	//Zr = 1 + 1 cell
    TCanvas *j4 = new TCanvas("j4","j4",900,700);
    Ejets_mu_cand2_ene->SetTitle("Energy of muons candidates (Zr = 1 + 1 cell, QF<" + SQFcut + ")");
    Ejets_mu_cand2_ene->Draw(); 
    Ejets_mu_cand2_eMF->Draw("SAME");
    leg->Draw("SAME");
    gPad->Update();
    j4->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand2_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j4->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand2_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    j4->SetLogy();
    gPad->Update();
    j4->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand2_log_ene_eMFQF" + SQFcut + ".pdf","pdf");
    j4->SaveAs("QF" + SQFcut + "/Muons/Ejets_mu_cand2_log_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fj4 = new TFile("QF" + SQFcut + "/MuonsEjets_mu_cand2_QF" + SQFcut + ".root","RECREATE");
    Ejets_mu_cand2_eMF->Write();
    Ejets_mu_cand2_ene->Write();
    fj4->Write();
    fj4->Close();
    // j4->Close();


	//Number of muons per event
	//Zr = 1 or 1 + 1 cell
    TCanvas *j7 = new TCanvas("j7","j7",900,700);
    Nmu_per_event_eMF->SetTitle("Number of muons candidates per event (Zr = 1 or 1 + 1 cell, QF<" + SQFcut + ")");
    Nmu_per_event_eMF->Draw(); 
    Nmu_per_event_ene->Draw("SAME");
    leg->Draw("SAME");
    //j2->SetLogy();
    if (QFcut == 50) j7->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) j7->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TCanvas *c159 = new TCanvas("c159","c159",900,700);
    Nmu_per_event_ene_eMF->Draw();
    if (QFcut == 50) c159->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_per_event_ene_versus_eMF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) c159->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_per_event_ene_versus_eMF" + SQFcut + ".jpeg","jpeg");
    //c159->Close();
    TFile *fj7 = new TFile("QF" + SQFcut + "/Muons/Njets_mu_cand_QF" + SQFcut + ".root","RECREATE");
    Nmu_per_event_eMF->Write();
    Nmu_per_event_ene->Write();
    fj7->Write();
    fj7->Close();
    //j7->Close();

	//Zr = 1
    TCanvas *j3 = new TCanvas("j3","j3",900,700);
    Nmu_per_event1_eMF->SetTitle("Number of muons candidates per event (Zr = 1, QF<" + SQFcut + ")");
    Nmu_per_event1_eMF->Draw(); 
    Nmu_per_event1_ene->Draw("SAME");
    leg->Draw("SAME");
    //j2->SetLogy();
    if (QFcut == 50) j3->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand1_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) j3->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand1_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fj3 = new TFile("QF" + SQFcut + "/Muons/Njets_mu_cand1_QF" + SQFcut + ".root","RECREATE");
    Nmu_per_event1_eMF->Write();
    Nmu_per_event1_ene->Write();
    fj3->Write();
    fj3->Close();
    //j3->Close();

	//Zr = 1 + 1 cell
    TCanvas *j5 = new TCanvas("j5","j5",900,700);
    Nmu_per_event2_ene->SetTitle("Number of muons candidates per event (Zr = 1 + 1 cell, QF<" + SQFcut + ")");
    Nmu_per_event2_ene->Draw(); 
    Nmu_per_event2_eMF->Draw("SAME");
    leg->Draw("SAME");
    j5->SetLogy();
    if (QFcut == 50) j5->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand2_ene_eMFQF" + SQFcut + ".pdf","pdf");
    if (QFcut == 50) j5->SaveAs("QF" + SQFcut + "/Muons/Njets_mu_cand2_ene_eMFQF" + SQFcut + ".jpeg","jpeg");
    TFile *fj5 = new TFile("QF" + SQFcut + "/Muons/Njets_mu_cand2_QF" + SQFcut + ".root","RECREATE");
    Nmu_per_event2_eMF->Write();
    Nmu_per_event2_ene->Write();
    fj5->Write();
    fj5->Close();
    //j5->Close();
*/
    //---------------------------------------------QF_analysis------------------------------------------------------
    TCanvas *QF1 = new TCanvas("QF1","QF1",900,700);
    QFOF2xQFCOF_1->Draw(); 
    QF1->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_1_OF2.pdf","pdf");
    QF1->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_1_OF2.jpeg","jpeg");
    //QF1->Close();

    TCanvas *QF2 = new TCanvas("QF2","QF2",900,700);
    QFOF2xQFCOF_2->Draw(); 
    QF2->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_2_OF2.pdf","pdf");
    QF2->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_2_OF2.jpeg","jpeg");
    //QF2->Close();

    TCanvas *QF3 = new TCanvas("QF3","QF3",900,700);
    QFOF2xQFCOF_3->Draw(); 
    QF3->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_3_OF2.pdf","pdf");
    QF3->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_3_OF2.jpeg","jpeg");
    //QF3->Close();

    TCanvas *QF4 = new TCanvas("QF4","QF4",900,700);
    QFOF2xQFCOF_4->Draw(); 
    QF4->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_4_OF2.pdf","pdf");
    QF4->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_4_OF2.jpeg","jpeg");
    //QF4->Close();

    TCanvas *QF5 = new TCanvas("QF5","QF5",900,700);
    QFOF2xQFCOF_5->Draw(); 
    QF5->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_5_OF2.pdf","pdf");
    QF5->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_5_OF2.jpeg","jpeg");
    //QF5->Close();

    TCanvas *QF6 = new TCanvas("QF6","QF6",900,700);
    QFOF2xQFCOF_6->Draw(); 
    QF6->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_6_OF2.pdf","pdf");
    QF6->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_6_OF2.jpeg","jpeg");
    //QF6->Close();

    TCanvas *QF7 = new TCanvas("QF7","QF7",900,700);
    QFOF2xQFCOF_7->Draw(); 
    QF7->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_7_OF2.pdf","pdf");
    QF7->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_7_OF2.jpeg","jpeg");
    //QF7->Close(); 

    TCanvas *QF8 = new TCanvas("QF8","QF8",900,700);
    QFOF2xQFCOF_8->Draw(); 
    QF8->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_8_OF2.pdf","pdf");
    QF8->SaveAs("QFanalysis/chi2_COF_VERSUS_chi2_8_OF2.jpeg","jpeg");
    //QF8->Close();    

	TCanvas *QF9 = new TCanvas("QF9","QF9",900,700);
	E_OF2_chi2->Draw();
	E_COF_chi2->Draw("SAME");
	leg->Draw("SAME");
	QF9->SaveAs("QFanalysis/Energy_channel_superimp_QF" + SQFcut_OF2 + ".pdf","pdf");
	QF9->SaveAs("QFanalysis/Energy_channel_superimp_QF" + SQFcut_OF2 + ".jpeg","jpeg");
    TFile *QF91 = new TFile("QFanalysis/Energy_channel_QF" + SQFcut_OF2 + ".root","RECREATE");
    E_OF2_chi2->Write();
    E_COF_chi2->Write();
    QF91->Write();
    QF91->Close();
    QF9->Close();

	TCanvas *QF10 = new TCanvas("QF10","QF10",900,700);
    E_OF2_versus_E_COF_chi2->Draw();
	QF10->SaveAs("QFanalysis/Energy_channel_QF" + SQFcut_OF2 + ".pdf","pdf");
	QF10->SaveAs("QFanalysis/Energy_channel_QF" + SQFcut_OF2 + ".jpeg","jpeg");
	QF10->Close();


}
