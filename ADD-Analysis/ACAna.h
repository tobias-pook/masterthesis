#ifndef ACAna_h
#define ACAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TEventList.h>
#include <TRandom3.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include "HistClass.h"

//~ #include "Combinations.h"

#include "TreeContent.h"
//#include <cppunit/extensions/HelperMacros.h>
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
//#include "LumiReweightingStandAlone.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "RunLumiRanges.h"
using std::cout;
using std::endl;
using namespace std;
//using namespace reweight;
//using namespace edm;

class ACAna : public TreeContent {
  
public :
	TString trigname;
    ACAna(TTree *tree=0);
    void Loop(TString fout, bool debug, TString *type);

    bool DEBUG;

    void init(int stages,TString *type);
    void write(TString fout, int nstages , TString *type);

    bool find_duplicate(int run, int evt, double x1, double x2);
    void BasicDump(int i);
	void TriggerDump(TString sel);
	void MuonDump(bool full=1);
	void PFJetDump();
	void TruthJetDump();
	void TruthDump();
	void VertexDump();
	void METDump();
	void SCDump();
	void EleDump(bool full=1);
	void PFEleDump(bool full=1);

    void MuonCuts();
    void Muon2Cuts(std::vector<int> goodMuons);
    void ElectronCutter();
    Double_t DeltaPhi(double a, double b);
    Double_t mT(double et1, double phi1, double et2, double phi2);

    Double_t AlphaT(const std::vector<TLorentzVector> & objects);
    Double_t MinDEt(const std::vector<TLorentzVector> & objects, 
          std::vector<UInt_t> * lista = NULL, 
          std::vector<UInt_t> * listb = NULL);
    Double_t SumET(const std::vector<TLorentzVector> & objects);
    Double_t MT(const std::vector<TLorentzVector> & objects);
    Double_t MT(Double_t*Et, Double_t* px, Double_t* py);
    Double_t Minv(Double_t*E, Double_t* px, Double_t* py,Double_t* pz);
    Double_t Minv(Double_t E, Double_t E1, Double_t px,Double_t px1, Double_t py,Double_t py1,  Double_t pz,Double_t pz1);
    Double_t Minv(TLorentzVector *particle[],Int_t Number);
    Double_t MT(TLorentzVector *particle[],Int_t Number);
    
    bool filterNoise();
    void readCut(Double_t* y_data, int ndata_max);
    
    void FillHist(Int_t hcut, TString *type);

    void FillHist2Mu(Int_t hcut,std::vector<int> goodMuons, TString *type, bool save_bool=0);
    void FillHistTruth(Int_t hcut,std::vector<int> goodMuons, double mZ);

    void FillHistNm1(Int_t hcut, TLorentzVector *imuon, Int_t jmuon, TString *type);

    
    bool MatchMuons(Int_t index, TLorentzVector *kmuon );

    void doCounter(bool cutbool2[], int hcut2);
    void doNm1Cuts(bool cutbool2[],double weight, TString *type);

    
    bool opositeSign(Int_t imuon,Int_t jmuon);
    typedef std::pair< pair<int,int> , pair<double,double> > Key;
    typedef std::set<Key> KeySet;
    typedef KeySet::const_iterator KeyIter;
    TRandom3 *rmd;
    JetCorrectionUncertainty *m_jecUnc;
    
  //~ edm::Lumi3DReWeighting LumiWeights_;
  //~ edm::Lumi3DReWeighting LumiWeights_up; 
  //~ edm::Lumi3DReWeighting LumiWeights_down;
	lumi::RunLumiRanges *runcfg; 
	
	edm::LumiReWeighting LumiWeights2;
	TString triggerRun;
   bool isTrigger; 
    Double_t muonmass;
    Double_t weight;
    Double_t vtxweight;
    Double_t pileupweightup;
    Double_t pileupweightdown;
    Double_t mZ;
    TFile* inputFile_;
    TH2* lut_;
    Int_t goodvtx;
    Int_t globalMuN;
    Int_t btagedJets;
    Int_t NJetspt50;
    std::vector<Int_t>  indexBtagedJets;
    Int_t btagedJets2;
    TString pickedEvents[300];
    Int_t pickedEventsCounter;
    Double_t deltaRmuon[100];
    bool cutmet[100];
    bool cutPt[100];
    //bool cutEta[100];
    bool cutCosmicAngle[100];

	bool cut_mu_ptrans[100];
	bool cut_mu_eta[100];
	bool cut_mu_id1[100];
	bool cut_muo_trkiso_zu_pt[100];
	bool cut_muo_TrackerLayersMeasTk[100];
	bool cut_muo_TrackerLayersMeasTkTight[100];
	bool cut_muo_StationsMatched[100];
	bool cut_muo_ValidMuonHitsCm[100];
	bool cut_muo_ValidPixelHitsCm[100];
	bool cut_trigger[100];
	bool cut_muo_d0[100];
	bool cut_muo_dzTk[100];
    
    bool cut_eta_2mu;
    bool cut_trigger_2mu;
    bool cut_cosmic_2mu;
    bool cut_charge_2mu;
    
    bool filter_passed[6];
    bool filter_tracking;
    bool filter_rest;
    
    Double_t cut[13];

	bool isTriggerDimuon;
	bool isTriggerSingle;
	bool isInAcceptance;
    bool isFullSelected;
 
    KeySet _keys;
    Int_t Zmuo[4];
 
    Int_t jetN;
    Int_t jetSoftN;
    Double_t sumJetPt;
    Double_t pxJet;
    Double_t pyJet;

    //pair <TString, pair<Double_t, Double_t> > ra[200];
    //std::vector<pair <TLorentzVector*, Int_t> > muon;
    //std::vector<pair <TLorentzVector*, Int_t> > muonAfterCut;
    //std::vector<pair <TLorentzVector*, Int_t> > muonForAnalyse;

    //Double_t met[3];
    TVector3* dEle[100];
    TLorentzVector* muon[100];
    TLorentzVector* muon_smeared[100];
    TLorentzVector* muon_scaled[100];
    TLorentzVector* muonAfterCut[100];
    TLorentzVector* muonAfterCutNodR[100];
    TLorentzVector *hardJet[100];
    std::vector< pair<int,int> > dRJet;
    //TLorentzVector* elec[100];

	bool filled[30];
    
    TRandom3 generator;
    
    std::vector<HistClass> Allhists;
    std::vector<HistClass> Allhists2D;


    std::vector<TH1F*> h1_nvtx;
    std::vector<TH1F*> h1_nvtx_weighted;    
    
    std::vector<TH1F*> h1_mu_pt;
    std::vector<TH1F*> h1_mu_pt_weighted;
    std::vector<TH1F*> h1_mu_phi;
    std::vector<TH1F*> h1_mu_phi_weighted;
    std::vector<TH1F*> h1_mu_eta;
    std::vector<TH1F*> h1_mu_eta_weighted;
    std::vector<TH1F*> h1_mu_n;
    
    std::vector<TH1F*> h1_2mu_n;
    std::vector<TH1F*> h1_2mu_n_weighted;
    std::vector<TH1F*> h1_2mu_mass;
    std::vector<TH1F*> h1_2mu_mass_weighted;
    std::vector<TH1F*> h1_2mu_mass_BB;
    std::vector<TH1F*> h1_2mu_mass_BB_weighted;
    std::vector<TH1F*> h1_2mu_mass_EB;
    std::vector<TH1F*> h1_2mu_mass_EB_weighted;
    std::vector<TH1F*> h1_2mu_mass_EE;
    std::vector<TH1F*> h1_2mu_mass_EE_weighted;
    std::vector<TH1F*> h1_2mu_mass_smeared;
    std::vector<TH1F*> h1_2mu_mass_smeared_weighted;
    std::vector<TH1F*> h1_2mu_mass_scaled;
    std::vector<TH1F*> h1_2mu_mass_scaled_weighted;
    std::vector<TH1F*> h1_2mu_eta;
    std::vector<TH1F*> h1_2mu_eta_weighted;
    std::vector<TH1F*> h1_2mu_phi;
    std::vector<TH1F*> h1_2mu_phi_weighted;  
    std::vector<TH1F*> h1_2mu_barrelendcap;
    std::vector<TH1F*> h1_2mu_barrelendcap_weighted;  
    std::vector<TH1F*> h1_2mu_pt;
    std::vector<TH1F*> h1_2mu_pt_weighted;  
    std::vector<TH1F*> h1_mZ; 
    
    std::vector<TH1F*> h1_2mu_muon_eta_truth;
    std::vector<TH1F*> h1_2mu_mass_truth;
    std::vector<TH1F*> h1_2mu_mass_truth_weighted;
};
#endif // #ifdef ACAna_cxx
