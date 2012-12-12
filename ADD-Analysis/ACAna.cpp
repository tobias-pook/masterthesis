 #define ACAna_cxx
#include "ACAna.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TSystem.h"


#include "SUSYDump.h"

using std::cout;
using std::endl;
using namespace std;

ACAna::ACAna(TTree *tree) : TreeContent(tree) {

  if (tree ==  0) {
    cerr << "ACAna: ERROR - Tree empty !!!" << endl;
    exit(0);
  }

}

// main event loop
void ACAna::Loop(TString fout, bool debug, TString *type) {
    
  DEBUG = debug;

  ProcInfo_t info;

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  
  cout << "ACAna: Running over " << nentries << " events" << endl;

  if (DEBUG) cout << "ACAna: DEBUG modus " << endl;

  if (type->Contains("data")) 
    cout << "ACAna: You are running with flag 'data'" << endl;
  else if (type->Contains("mc") || type->Contains("dataDriven"))
    cout << "ACAna: You are running with flag 'mc' "<<*type << endl;
  else 
    cout << "ACAna: No flag specified !" << endl;

  TString trigger;  

  //define Number of cuts for histograms
  Int_t ncut = 26;
  //To get JetCorrectionUncertainty:
  
  init(ncut,type);



  

        //

  //

        
  Long64_t nbytes = 0, nb = 0;
  //Int_t SCMatchCounter=0;


  
 
  
  // main event loop
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) cout<<"all wrong"<<endl;//break;
    nb = fChain->GetEntry(jentry);nbytes += nb;
      
    if (DEBUG && jentry == 100000) break;
    
    if (fmod((double)jentry,(double)10000)==0) 
      cout << "Event " << jentry << endl;

    if (fmod((double)jentry,(double)50000)==0){
      gSystem->GetProcInfo(&info);
      cout << " -> Memory in MB : " << info.fMemResident/1000. << " (resident) " 
	   << info.fMemVirtual/1000. << " (virtual) " << endl;
    }
    
    
     if(type->Contains("data") && type->Contains("golden")){
     //check if we shall analyze this event
         lumi::ID run = global_run;
         lumi::ID LS  = lumi_section;

             if( ! (*runcfg).check( run, LS ) ) {
                //cout<<"skipping event"<<endl;
                continue;
             }
        }
          
   //~ cout << "N-filters " << eventfilter_n << endl;
   filter_rest=1;
   filter_tracking=1;
   //~ if (type->Contains("data")){
	   //~ for (int o=0; o<eventfilter_n; o++){
			//cout << o << " " << unpack((int*)eventfilter_names[o]) << eventfilter_results[o] << endl;
			//~ filter_passed[o]=eventfilter_results[o];
			//~ if (o<(eventfilter_n-1) && filter_passed[o]==0) filter_rest=0; 
			//~ if (o==(eventfilter_n-1) && filter_passed[o]==0) filter_tracking=0; 
			//~ 
	   //~ }
	//~ }
	if (type->Contains("data")){
		for(Int_t o=0; o<eventfilter_n;++o){

			if(!eventfilter_results[o] && o!=5 ){

				filter_rest=0;
			}
		}
	}
    //if(!isMETFILTER && type.Contains("filter") ){
    
    
   for (Int_t jj=0; jj<ncut; jj++) {
	 filled[jj]=false;
	}
    weight=1;
    pileupweightup=1.;
    pileupweightdown=1.;
    
	if(type->Contains("pileup")){
        vtxweight= LumiWeights2.weight((float)pu_TrueNrInter);
        weight=weight*vtxweight;
        
    }
	double mZmax=-1.;
	double mZ=-1.;
	for (int tri=0; tri<truth_n; tri++){
		if (abs(truth_pdgid[tri])==22 || abs(truth_pdgid[tri])==23) {
			mZ=truth_m[tri]; 
			if (mZ > mZmax) mZmax=mZ;
		}	
	}
	mZ=mZmax;	
	
	h1_mZ[0]->Fill(mZ);
	if (type->Contains("DY")){
					if (type->Contains("DY20") && !(type->Contains("DY200")) && !(type->Contains("DY2000"))){
						if (mZ > 120.) continue;
					}
					if (type->Contains("DY120")){
						if (mZ < 120. || mZ > 200) continue;
					}		
					if (type->Contains("DY200") && !(type->Contains("DY2000"))){
						if (mZ < 200. || mZ > 500) continue;
					}
					if (type->Contains("DY500")){
						if (mZ < 500. || mZ > 800) continue;
					}												
					if (type->Contains("DY800")){
						if (mZ < 800. || mZ > 1000) continue;
					}					
					if (type->Contains("DY1000")){
						if (mZ < 1000. || mZ > 1500) continue;
					}
					if (type->Contains("DY1500")){
						if (mZ < 1500. || mZ > 2000) continue;
					}					
					if (type->Contains("DY2000")){
						if (mZ < 2000.) continue;
					}
					h1_mZ[1]->Fill(mZ);											
	}	 



	if (type->Contains("pythia")){
					if (type->Contains("pythia20") && !(type->Contains("pythia200"))){
						if (mZ > 120.) continue;
					}
					if (type->Contains("pythia120")){
						if (mZ < 120. || mZ > 200) continue;
					}		
					if (type->Contains("pythia200")){
						if (mZ < 200. || mZ > 500) continue;
					}
					if (type->Contains("pythia500")){
						if (mZ < 500. || mZ > 800) continue;
					}												
					if (type->Contains("pythia800")){
						if (mZ < 800. || mZ > 1300) continue;
					}					
					if (type->Contains("pythia1300")){
						if (mZ < 1300. || mZ > 1600) continue;
					}	
					if (type->Contains("pythia1600")){
						if (mZ < 1600.) continue;
					}
					h1_mZ[1]->Fill(mZ);											
	}	






	
	if (type->Contains("Tau")){
					if (type->Contains("Tau20") && !(type->Contains("Tau200"))){
						if (mZ > 100.) continue;
					}
					if (type->Contains("Tau100")){
						if (mZ < 100. || mZ > 200) continue;
					}		
					if (type->Contains("Tau200")){
						if (mZ < 200. || mZ > 400) continue;
					}
					if (type->Contains("Tau400")){
						if (mZ < 400. || mZ > 800) continue;
					}												
					if (type->Contains("Tau800")){
						if (mZ < 800.) continue;
					}					

					h1_mZ[1]->Fill(mZ);											
	}	 

	
    goodvtx=0;
    
	for (int bb=0; bb<vtx_n; bb++){
		  if(vtx_fake[bb]==0){
			  if(vtx_ndof[bb]>4.){
				  if(fabs(vtx_z[bb])<24.){
					goodvtx++;	
						   }
				  }
		   }

	}   
    
    FillHist(21, type);
    
    triggerRun="HLT_Mu40";
    isTrigger=false;
    isTriggerDimuon=false;
    isTriggerSingle=false;
	isInAcceptance=false;
	isFullSelected=false;
	for(Int_t i=0; i<trig_n;++i){
        trigger = unpack((int*)trig_name[i]);
        //~ cout << trigger << endl;
        //~ if(trigger.Contains("HLT_Mu17")) cout << trigger << " " << trig_L1prescale[i] << " " << trig_HLTprescale[i]<< endl;
		//~ if(trigger.Contains("HLT_Mu17_Mu8")) {
		if(trigger.Contains("HLT_Mu17_TkMu8")) {
            isTriggerDimuon=true;
			//~ cout << "found dimoun" << endl;
		}    
		if(trigger.Contains(triggerRun.Data())) {
			isTriggerSingle=true;
            isTrigger=true;

		}    
    }
    //~ isTrigger=true;
	for(int i=0; i<100;++i){
		muon[i]= new TLorentzVector();
		muon_smeared[i]= new TLorentzVector();
		muon_scaled[i]= new TLorentzVector();


		cutCosmicAngle[i] = true;

		cut_mu_ptrans[i] = false;
		cut_mu_eta[i] = false;
		cut_mu_id1[i] = false;
		cut_muo_trkiso_zu_pt[i] = false;
		cut_muo_TrackerLayersMeasTk[i] = false;
		cut_muo_TrackerLayersMeasTkTight[i] = false;
		cut_muo_StationsMatched[i] = false;
		cut_muo_ValidMuonHitsCm[i] = false;
		cut_muo_ValidPixelHitsCm[i] = false;
		cut_trigger[i] = false;
		cut_muo_d0[i] = false;
		cut_muo_dzTk[i] = false;
	}
	

	cut_eta_2mu=false;
    cut_trigger_2mu=false;
    cut_cosmic_2mu=false;
    cut_charge_2mu=false;
	MuonCuts();

    if (type->Contains("trigger")) {
		isTrigger=true;		
		for(int i=0; i<100;++i){
			cut_trigger[i]=true;
		}
	}

    

    
    std::vector<int> goodMuons;   
    std::vector<int> goodMuons_pre;   
    int nmuo=0;
    int hcut=0;
    for (int jmuon=0; jmuon<muo_n; jmuon++){
		if (cut_mu_id1[jmuon]) 	{
					muon[jmuon]->SetPtEtaPhiM(muo_Cocktail_pt[jmuon],muo_Cocktail_eta[jmuon],muo_Cocktail_phi[jmuon],0.105658);
					//~ double newx=muon[jmuon]->Px()*generator.Gaus(1.,0.1);
					//~ double newy=muon[jmuon]->Py()*generator.Gaus(1.,0.1);
					//~ muon_smeared[jmuon]->SetPxPyPzE(newx,newy,muon[jmuon].Pz(),sqrt(newx*newx+newy*newy+muon[jmuon].Pz()*muon[jmuon].Pz()+muonmass*muonmass);
					muon_smeared[jmuon]->SetPtEtaPhiM(muo_Cocktail_pt[jmuon]*generator.Gaus(1.,0.03),muo_Cocktail_eta[jmuon],muo_Cocktail_phi[jmuon],0.105658);
					double scalefac=0.05*0.001*muon[jmuon]->Pt()+1.;
					muon_scaled[jmuon]->SetPtEtaPhiM(muo_Cocktail_pt[jmuon]*scalefac,muo_Cocktail_eta[jmuon],muo_Cocktail_phi[jmuon],0.105658);
					//~ cout << generator.Gaus(1.,0.1)<<endl;	
					
			}		
		else {
					muon[jmuon]->SetPtEtaPhiM(muo_pt[jmuon],muo_eta[jmuon],muo_phi[jmuon],0.105658);
					//~ muon_smeared[jmuon]->SetPtEtaPhiM(muon[jmuon]->Pt()*generator.Gaus(1.,0.1),muon[jmuon]->Eta(),muon[jmuon]->Phi(),0.105658);
					muon_smeared[jmuon]->SetPtEtaPhiM(muon[jmuon]->Pt()*generator.Gaus(1.,0.03),muon[jmuon]->Eta(),muon[jmuon]->Phi(),0.105658);
					double scalefac=0.05*0.001*muon[jmuon]->Pt()+1.;
					muon_scaled[jmuon]->SetPtEtaPhiM(muon[jmuon]->Pt()*scalefac,muon[jmuon]->Eta(),muon[jmuon]->Phi(),0.105658);
			}		
		hcut=0;	

		FillHistNm1(hcut,muon[jmuon],jmuon,type);
		hcut++;
		
		if (isTrigger==true){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}	
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}	
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}		
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}
		hcut++;		
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}		
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon] &&
			cut_muo_ValidMuonHitsCm[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}	
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon] &&
			cut_muo_ValidMuonHitsCm[jmuon] &&
			cut_muo_ValidPixelHitsCm[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}	
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon] &&
			cut_muo_ValidMuonHitsCm[jmuon] &&
			cut_muo_ValidPixelHitsCm[jmuon] &&
			cut_muo_d0[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}	
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon] &&
			cut_muo_ValidMuonHitsCm[jmuon] &&
			cut_muo_ValidPixelHitsCm[jmuon] &&
			cut_muo_d0[jmuon] &&
			cut_muo_dzTk[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
		}			
		hcut++;
		if (isTrigger==true &&
			cut_mu_ptrans[jmuon] &&
			cut_mu_eta[jmuon] &&
			cut_mu_id1[jmuon] &&
			cut_muo_TrackerLayersMeasTk[jmuon] &&
			cut_muo_StationsMatched[jmuon] &&
			cut_muo_ValidMuonHitsCm[jmuon] &&
			cut_muo_ValidPixelHitsCm[jmuon] &&
			cut_muo_d0[jmuon] &&
			cut_muo_dzTk[jmuon]&&
			cut_muo_trkiso_zu_pt[jmuon]){
			FillHistNm1(hcut,muon[jmuon],jmuon,type);
			goodMuons_pre.push_back(jmuon);
			nmuo++;
		}			
		
	}
	hcut=12;
	
	FillHist2Mu(hcut,goodMuons,type);
	//~ if (muo_n >2 && goodMuons_pre.size()==2) goodMuons_pre.push_back(3);
	if (goodMuons_pre.size()>2){
		//~ cout << "More than 2 muons found" <<endl;
		double mass_max=0.;
		int mu1_max=-1;
		int mu2_max=-1;
		
		for (int mu1=0; mu1<(goodMuons_pre.size()-1); mu1++){
			for (int mu2=mu1+1; mu2<goodMuons_pre.size(); mu2++){

				TLorentzVector Z=(*(muon[mu1]))+(*(muon[mu2]));
				double mass_pair=Z.M();
				//~ cout << mu1 << " " << mu2 << " " << mass_pair <<endl;
				if (mass_pair>mass_max){
						mass_max=mass_pair;
						mu1_max=mu1;
						mu2_max=mu2;						
				}
			}	
		}
		goodMuons.push_back(mu1_max);
		goodMuons.push_back(mu2_max);
	}
	else goodMuons=goodMuons_pre;


	hcut++;  //13
	if (goodMuons.size()==2){
		Muon2Cuts(goodMuons);
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //14
	if (goodMuons.size()==2 &&
		cut_eta_2mu){
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //15
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu){
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //16
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu){
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //17
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu &&
		cut_cosmic_2mu){
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //18
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu &&
		cut_cosmic_2mu &&
		goodvtx > 0.){
		//~ isFullSelected=true;	
		FillHist2Mu(hcut,goodMuons,type);
	}      
	hcut++;  //19
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu &&
		cut_cosmic_2mu &&
		goodvtx > 0. &&
		filter_rest==1){
		//~ isFullSelected=true;	
		FillHist2Mu(hcut,goodMuons,type);
	}
	hcut++;  //20
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu &&
		cut_cosmic_2mu &&
		goodvtx > 0. &&
		filter_rest==1 &&
		filter_tracking==1){
		isFullSelected=true;	
		FillHist2Mu(hcut,goodMuons,type, 1);
	}	
	
	if (goodMuons.size()==2 &&
		cut_eta_2mu &&
		cut_trigger_2mu &&
		cut_charge_2mu &&
		cut_cosmic_2mu &&
		goodvtx > 0. &&
		filter_rest==1 &&
		filter_tracking==1){
		if (cut_muo_TrackerLayersMeasTkTight[goodMuons[0]] &&
		cut_muo_TrackerLayersMeasTkTight[goodMuons[1]]){
				FillHist2Mu(25,goodMuons,type);
		}		
	}	
if (type->Contains("trigger")){
		hcut++; //21
		if (goodMuons.size()==2 &&
			cut_eta_2mu &&
			cut_trigger_2mu &&
			cut_charge_2mu &&
			cut_cosmic_2mu &&
			goodvtx > 0. &&
			isTriggerSingle==true){
			FillHist2Mu(hcut,goodMuons,type);
		}     
		hcut++; //22
		if (goodMuons.size()==2 &&
			cut_eta_2mu &&
			cut_trigger_2mu &&
			cut_charge_2mu &&
			cut_cosmic_2mu &&
			goodvtx > 0. &&
			isTriggerDimuon==true){
			FillHist2Mu(hcut,goodMuons,type);
		}     

}
	// CLEANUP
    for(int i=0; i<100;++i){
		delete muon[i];
		delete muon_smeared[i];
		delete muon_scaled[i];
	}
	// TRUTH INFORMATION 
	double massmax=-1;
	int muon1=-1;
	int muon2=-1;
	vector<int> goodTruthMuons;
	vector<int> AnyTruthMuons;
	for (int tri=0; tri<truthl_n; tri++){
			for (int tri2=tri+1; tri2<truthl_n; tri2++){
				
				
				if (fabs(truthl_pdgid[tri])==13 && fabs(truthl_pdgid[tri2])==13){
					TLorentzVector mm1;
					TLorentzVector mm2;
					AnyTruthMuons.push_back(tri);
					AnyTruthMuons.push_back(tri2);
					FillHistTruth(7,AnyTruthMuons,mZ);
					mm1.SetPtEtaPhiM(truthl_pt[tri],truthl_eta[tri],truthl_phi[tri],muonmass);
					mm2.SetPtEtaPhiM(truthl_pt[tri2],truthl_eta[tri2],truthl_phi[tri2],muonmass);
					TLorentzVector Z=(mm1)+(mm2);
					if (Z.M()>massmax) {
							massmax=Z.M();
							muon1=tri;
							muon2=tri2;
					}
				}
			}
	}
			goodTruthMuons.push_back(muon1);
			goodTruthMuons.push_back(muon2);	
	
				
		if (massmax>0.){

				FillHistTruth(0,goodTruthMuons,mZ);
				if (fabs(truthl_eta[muon1])<2.4 && fabs(truthl_eta[muon2])<2.4){
						if(fabs(truthl_eta[muon1])<2.1 || fabs(truthl_eta[muon2])<2.1){
							FillHistTruth(1,goodTruthMuons,mZ);
							//~ if (truthl_pt[muon1]>45. && truthl_pt[muon2]>45.&& muo_pt[0] > 45. && muo_pt[1]>45.){
							if (truthl_pt[muon1]>45. && truthl_pt[muon2]>45.){
								FillHistTruth(2,goodTruthMuons,mZ);
								if (isTrigger) FillHistTruth(3,goodTruthMuons,mZ);
								if (isFullSelected) FillHistTruth(4,goodTruthMuons,mZ);
							}
						}
				}
				
		}
		//~ if (massmax>0.){
				//~ goodTruthMuons.push_back(muon1);
				//~ goodTruthMuons.push_back(muon2);
				//~ FillHistTruth(0,goodTruthMuons,massmax);
				//~ if (fabs(truthl_eta[muon1])<2.4 && fabs(truthl_eta[muon2])<2.4){
						//~ if(fabs(truthl_eta[muon1])<2.1 || fabs(truthl_eta[muon2])<2.1){
							//~ FillHistTruth(1,goodTruthMuons,massmax);
							//~ if (truthl_pt[muon1]>45. && truthl_pt[muon2]>45.){
								//~ FillHistTruth(2,goodTruthMuons,massmax);
								//~ if (isTrigger) FillHistTruth(3,goodTruthMuons,massmax);
							//~ }
						//~ }
				//~ }
				//~ if (isFullSelected) FillHistTruth(4,goodTruthMuons,massmax);
		//~ }				
		
	 
	  
    }  //entry loop 
      write(fout, ncut, type);
}

// duplicates finder
bool ACAna::find_duplicate(int run, int evt, double x1, double x2) {

  pair<int,int> temp1(run,evt);
  pair<double,double> temp2(x1,x2);
  Key key (temp1, temp2);
  
  KeyIter pos = _keys.find (key);
  
  if (pos == _keys.end()) {
    _keys.insert (key);
    return false;
  }
  else {
    if (DEBUG) cout << "ACAna: duplicate run " << run << " , evt " << evt << " , scale " << x1 << " , vtx_z " << x2 << endl;
    return true;
  }
    
}

// initialize histograms
void ACAna::init(int nstages, TString *type) {
    //Int_t NumberOfHist=122;
  
	muonmass=0.105658;
	 //~ if(type->Contains("Fall")){
		 //~ LumiWeights_=edm::Lumi3DReWeighting("histProbFunctionmidFall.root","histAllData.root","pileup","pileup");
			//~ if (type->Contains("sys")){
				//~ LumiWeights_up=edm::Lumi3DReWeighting("histProbFunctionmidFall.root","histAllData.root","pileup","pileup");
				//~ LumiWeights_down=edm::Lumi3DReWeighting("histProbFunctionmidFall.root","histAllData.root","pileup","pileup");
				//~ LumiWeights_up.weight3D_init(1.08);
				//~ LumiWeights_down.weight3D_init(0.92);
			//~ }
	 //~ } else {
		// LumiWeights_=edm::Lumi3DReWeighting("MC_2012.root","histAllData_Summer12_60vtx.root","pileup","pileup");
		 //~ LumiWeights_=edm::Lumi3DReWeighting("MC_2012.root","histAllData_Summer12_60vtx.root","pileup","pileup");
		// LumiWeights_=edm::Lumi3DReWeighting("MC_2012.root","data_true.root","pileup","pileup");
		// LumiWeights_=edm::Lumi3DReWeighting("PileUpSummer12.root","histAllData.root","pileup","pileup");
		  //~ if (type->Contains("sys")){
			//~ LumiWeights_up=edm::Lumi3DReWeighting("PileUpSummer12.root","histAllData.root","pileup","pileup");
			//~ LumiWeights_down=edm::Lumi3DReWeighting("PileUpSummer12.root","histAllData.root","pileup","pileup");
			//~ LumiWeights_up.weight3D_init(1.08);
			//~ LumiWeights_down.weight3D_init(0.92);
			//~ }
	 //~ }
  
  //~ LumiWeights_.weight3D_init(1.00);
	if(type->Contains("data")){
		if(type->Contains("golden")) runcfg=new lumi::RunLumiRanges("190456-204567_combined.txt");
		if(type->Contains("recover")) runcfg=new lumi::RunLumiRanges("Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON_MuonPhys.txt");
    } 
    
	LumiWeights2 = edm::LumiReWeighting("MC_2012.root","MyDataPileupHistogram.root","pileup","pileup");
	
    

    HistClass hist; 
    
    hist.init("h1_%i_mu_pt","p_{T}",8000, 0., 8000.            ,&h1_mu_pt );     Allhists.push_back(hist);
    hist.init("h1_%i_mu_pt_weighted","p_{T}",8000, 0., 8000.            ,&h1_mu_pt_weighted );     Allhists.push_back(hist);
    hist.init("h1_%i_mu_eta","#eta",800, -4., 4.            ,&h1_mu_eta );     Allhists.push_back(hist);
    hist.init("h1_%i_mu_eta_weighted","#eta",800, -4., 4.            ,&h1_mu_eta_weighted );     Allhists.push_back(hist);
    hist.init("h1_%i_mu_phi","#phi",640, -3.2, 3.2            ,&h1_mu_phi );     Allhists.push_back(hist);   
    hist.init("h1_%i_mu_phi_weighted","#phi",640, -3.2, 3.2            ,&h1_mu_phi_weighted );     Allhists.push_back(hist);   
    hist.init("h1_%i_nvtx","p_{T}",60, 0., 60.            ,&h1_nvtx);     Allhists.push_back(hist);
    hist.init("h1_%i_nvtx_weighted","p_{T}",60, 0., 60.            ,&h1_nvtx_weighted );     Allhists.push_back(hist);    
    
    hist.init("h1_%i_mu_n","N_{#mu}",30, 0., 30.            ,&h1_mu_n );     Allhists.push_back(hist);
    hist.init("h1_%i_2mu_n","N_{#mu}",30, 0., 30.            ,&h1_2mu_n );     Allhists.push_back(hist);
    hist.init("h1_%i_2mu_n_weighted","N_{#mu}  #",30, 0., 30.            ,&h1_2mu_n_weighted );     Allhists.push_back(hist);
    hist.init("h1_%i_2mu_mass","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_weighted );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_BB","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_BB );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_BB_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_BB_weighted );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_EB","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_EB );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_EB_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_EB_weighted );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_EE","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_EE );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_EE_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_EE_weighted );     Allhists.push_back(hist);    
    
    hist.init("h1_%i_2mu_mass_smeared","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_smeared );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_smeared_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_smeared_weighted );     Allhists.push_back(hist);
    hist.init("h1_%i_2mu_mass_scaled","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_scaled );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_scaled_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_scaled_weighted );     Allhists.push_back(hist);    
    
    
    hist.init("h1_%i_2mu_phi","#phi",640, -3.2, 3.2            ,&h1_2mu_phi );     Allhists.push_back(hist);     
    hist.init("h1_%i_2mu_phi_weighted","#phi",640, -3.2, 3.2            ,&h1_2mu_phi_weighted );     Allhists.push_back(hist);     
    hist.init("h1_%i_2mu_barrelendcap","#au",21, -0.5, 20.5            ,&h1_2mu_barrelendcap );     Allhists.push_back(hist);     
    hist.init("h1_%i_2mu_barrelendcap_weighted","#au",21, -0.5, 20.5            ,&h1_2mu_barrelendcap_weighted );     Allhists.push_back(hist);     
    
    hist.init("h1_%i_2mu_eta","#eta",800, -4., 4.            ,&h1_2mu_eta );     Allhists.push_back(hist);      
    hist.init("h1_%i_2mu_eta_weighted","#eta",800, -4., 4.            ,&h1_2mu_eta_weighted );     Allhists.push_back(hist);      
    hist.init("h1_%i_2mu_pt","p_{T} [GeV]",8000, 0., 8000.            ,&h1_2mu_pt );     Allhists.push_back(hist);
    hist.init("h1_%i_2mu_pt_weighted","p_{T} [GeV]",8000, 0., 8000.            ,&h1_2mu_pt_weighted );     Allhists.push_back(hist);
    
    hist.init("h1_%i_mZ","m [GeV]",8000, 0., 8000.            ,&h1_mZ );     Allhists.push_back(hist);

    hist.init("h1_%i_2mu_muon_eta_truth","#eta",800, -4., 4            ,&h1_2mu_muon_eta_truth );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_truth","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_truth );     Allhists.push_back(hist);    
    hist.init("h1_%i_2mu_mass_truth_weighted","m_{2#mu}",8000, 0., 8000.            ,&h1_2mu_mass_truth_weighted );     Allhists.push_back(hist);    
    Char_t h_name[200];
    Char_t h_title[200];
    TH1F* h_temp;
  strcpy(h_title,"\0");
  for (Int_t jj=0; jj<nstages; jj++) {
	  
	  
	   for(UInt_t ii=0; ii<Allhists.size(); ++ii){
    
		   strcpy(h_name,"\0");
			
			sprintf(h_name,Allhists[ii].HistName, jj);
			h_temp = new TH1F(h_name, h_title, Allhists[ii].nBins ,Allhists[ii].xMin , Allhists[ii].xMax );
			h_temp->Sumw2();
			h_temp->SetXTitle(Allhists[ii].XTitle);
			Allhists[ii].hist->push_back(h_temp);
    }
	
 }	  
	
  
}

// write histograms to file
void ACAna::write(TString fout, int nstages, TString *type) {
  
  TFile *f = new TFile(fout,"RECREATE");
  f->SetCompressionLevel(9);
  f->cd();
	for (Int_t jj=0; jj<nstages; jj++) {
        for(UInt_t ii=0; ii<Allhists.size(); ++ii){
            //if((*(Allhists[ii].hist))[jj]->Integral()>0.001)
                (*(Allhists[ii].hist))[jj]->Write();
        }
        //~ for(UInt_t ii=0; ii<Allhists2D.size(); ++ii){
            //~ (*(Allhists2D[ii].hist2))[jj]->Write();
        //~ }
   }
  f->Close();     
    

}

//FillHist(hcut, kmuon, muonAfterCutIndex[imuon], cutbool, muonAfterCut, qualityN, met );

void ACAna::FillHist(Int_t hcut, TString *type){
	if (filled[hcut]==0) {
			h1_nvtx[hcut]->Fill(goodvtx);
			h1_nvtx_weighted[hcut]->Fill(goodvtx,weight);
   
   
		filled[hcut]=1;
	}
}


void ACAna::FillHistNm1(Int_t hcut, TLorentzVector *imuon, Int_t jmuon, TString *type){
	  FillHist(hcut,type);
      h1_mu_pt_weighted[hcut]->Fill(imuon->Pt(),weight);
      h1_mu_pt[hcut]->Fill(imuon->Pt());
      
      h1_mu_eta_weighted[hcut]->Fill(imuon->Eta(),weight);
      h1_mu_eta[hcut]->Fill(imuon->Eta());     
      h1_mu_phi_weighted[hcut]->Fill(imuon->Phi(),weight);
      h1_mu_phi[hcut]->Fill(imuon->Phi());      
}


void ACAna::FillHist2Mu(Int_t hcut,std::vector<int> goodMuons, TString *type, bool save_bool){
   
   for (unsigned int kmuon=0; kmuon<goodMuons.size(); kmuon++){
	  FillHistNm1(hcut, muon[goodMuons[kmuon]],goodMuons[kmuon],type);
   }
   if (goodMuons.size()>=2){
	   
		int barrel=0;
		int endcap=0;
		if (fabs(muon[goodMuons[0]]->Eta())<1.2) barrel++;
		else endcap++;
		if (fabs(muon[goodMuons[1]]->Eta())<1.2) barrel++;
		else endcap++;		
 	   
		TLorentzVector* muon0=muon[goodMuons[0]];
	    TLorentzVector* muon1=muon[goodMuons[1]];
	    TLorentzVector Z=(*muon0)+(*muon1);
	    
	    TLorentzVector Z_smeared=*(muon_smeared[goodMuons[0]])+*(muon_smeared[goodMuons[1]]);
	    TLorentzVector Z_scaled=*(muon_scaled[goodMuons[0]])+*(muon_scaled[goodMuons[1]]);
	    
	    h1_2mu_n_weighted[hcut]->Fill(goodMuons.size(),weight);
	    h1_2mu_n[hcut]->Fill(goodMuons.size());
	   
	    h1_2mu_mass[hcut]->Fill(Z.M());
	    h1_2mu_mass_weighted[hcut]->Fill(Z.M(),weight);
	    
	    if(barrel==2){
			h1_2mu_mass_BB[hcut]->Fill(Z.M());
			h1_2mu_mass_BB_weighted[hcut]->Fill(Z.M(),weight);
		} else if(barrel==1 && endcap==1){
			h1_2mu_mass_EB[hcut]->Fill(Z.M());
			h1_2mu_mass_EB_weighted[hcut]->Fill(Z.M(),weight);
		} else if(endcap==2){
			h1_2mu_mass_EE[hcut]->Fill(Z.M());
			h1_2mu_mass_EE_weighted[hcut]->Fill(Z.M(),weight);
		}
	    h1_2mu_mass_smeared[hcut]->Fill(Z_smeared.M());
	    h1_2mu_mass_smeared_weighted[hcut]->Fill(Z_smeared.M(),weight);
	    h1_2mu_mass_scaled[hcut]->Fill(Z_scaled.M());
	    h1_2mu_mass_scaled_weighted[hcut]->Fill(Z_scaled.M(),weight);
	   
	    h1_2mu_eta[hcut]->Fill(Z.Eta());
	    h1_2mu_eta_weighted[hcut]->Fill(Z.Eta(),weight);   

	    h1_2mu_phi[hcut]->Fill(Z.Phi());
	    h1_2mu_phi_weighted[hcut]->Fill(Z.Phi(),weight);     

	    h1_2mu_pt[hcut]->Fill(Z.Pt());
 	    h1_2mu_pt_weighted[hcut]->Fill(Z.Pt(),weight); 

	    h1_2mu_barrelendcap[hcut]->Fill(barrel+5*endcap);
 	    h1_2mu_barrelendcap_weighted[hcut]->Fill(barrel+5*endcap,weight); 


 	    
 	    if (Z.M() > 700. && save_bool==true && type->Contains("save")){
			     FILE *pFile;
			     pFile = fopen("save.txt","append");
			     fprintf(pFile,"Run: %i  Lumisection: %i  Event: %i \n",global_run,lumi_section,global_event);
		//~ 
				 fprintf(pFile,"Z: M: %g Eta: %g Phi: %g PT: %g \n",Z.M(), Z.Eta(), Z.Phi(), Z.Pt());
				 fprintf(pFile,"Chosen Muons: %i %i\n",goodMuons[0], goodMuons[1]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexMass[goodMuons[0]][goodMuons[1]]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][5]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][4]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][3]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][2]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][1]);
				 //~ fprintf(pFile,"Inv Mass Vertex Constraint: %g\n", muo_DiMuonVertexValues[goodMuons[0]][goodMuons[1]][0]);
				 fprintf(pFile,"Vertices: %i\n",goodvtx);
				 
				 
				for (Int_t i=0; i<muo_n; i++)
					{
				
						fprintf (pFile, "   Muon %i: pt:%g charge:%g pt(ohneCocktail):%g eta:%g phi:%g TrkIso:%g global:%i tracker:%i Chi2:%g hitsTrk:%i d0Cm:%g  d0Tk:%g d0Origin:%g\n",i,muo_Cocktail_pt[i],muo_charge[i],muo_pt[i],muo_eta[i],muo_phi[i],muo_TrkIso[i],muo_ID[i][1],muo_ID[i][3],muo_TrkChiNormCm[i],muo_hitsTk[i],muo_d0Cm[i],muo_d0Tk[i],muo_d0OriginCm[i]);
						fprintf (pFile, "   GlobalTrack: %g: Cocktail: %g Inner:%g TPFMS:%g Picky:%g DYT:%g Default: %g \n",muo_TevReco_pt[i][0],muo_TevReco_pt[i][1],muo_TevReco_pt[i][2],muo_TevReco_pt[i][3],muo_TevReco_pt[i][4],muo_TevReco_pt[i][5],muo_TevReco_pt[i][6]);
						fprintf (pFile, " Chi2/ndof  GlobalTrack: %g: Cocktail: %g Inner:%g TPFMS:%g Picky:%g DYT:%g Default: %g \n",muo_TevReco_chi2[i][0]/muo_TevReco_ndof[i][0],muo_TevReco_chi2[i][1]/muo_TevReco_ndof[i][1],muo_TevReco_chi2[i][2]/muo_TevReco_ndof[i][2],muo_TevReco_chi2[i][3]/muo_TevReco_ndof[i][3],muo_TevReco_chi2[i][4]/muo_TevReco_ndof[i][4],muo_TevReco_chi2[i][5]/muo_TevReco_ndof[i][5],muo_TevReco_chi2[i][6]/muo_TevReco_ndof[i][6]);
						fprintf (pFile, "TrackIso:%g, EcalIso:%g, HCalIso:%g, RelIso:%g\n",muo_TrkIso[i],muo_ECalIso[i],muo_HCalIso[i],(muo_TrkIso[i]+muo_ECalIso[i]+muo_HCalIso[i])/muo_pt[i]);
						fprintf (pFile, "LayersMeas: %i\n",muo_TrackerLayersMeasCm[i]);
						fprintf (pFile, "PFCandidate pt: %g, eta: %g, phi: %g, DeltaR: %g pdgid: %i \n \n",sqrt(pow(muo_PFCand_px[i],2)+pow(muo_PFCand_py[i],2)),muo_PFCand_eta[i],muo_PFCand_phi[i],muo_PFCand_DeltaR[i],muo_PFCand_pfid[i]);
						fprintf (pFile, "Standard: pt: %g ", muo_pt[i]);
						fprintf (pFile, "Global: pt: %g +/- %g \n ",muo_TevReco_pt[i][0],muo_TevReco_ptError[i][0]);
						fprintf (pFile, "Cocktail: pt: %g +/- %g \n ",muo_TevReco_pt[i][1],muo_TevReco_ptError[i][1]);
						fprintf (pFile, "Inner: pt: %g +/- %g \n ",muo_TevReco_pt[i][2],muo_TevReco_ptError[i][2]);     			
						fprintf (pFile, "TPFMS: pt: %g +/- %g \n ",muo_TevReco_pt[i][3],muo_TevReco_ptError[i][3]);     			
						fprintf (pFile, "Picky: pt: %g +/- %g \n ",muo_TevReco_pt[i][4],muo_TevReco_ptError[i][4]);     			
						fprintf (pFile, "DYT: pt: %g +/- %g \n ",muo_TevReco_pt[i][5],muo_TevReco_ptError[i][5]);     		
			
					}
				
				for (Int_t i=0; i<pfjet_n;i++)
				{
					if(pfjet_pt[i]>30.) {
						fprintf (pFile,"  Jet %i: pt:%f eta:%f phi:%f\n",i,pfjet_pt[i],pfjet_eta[i],pfjet_phi[i]); 
					}	
				}
				
				
					
				fprintf(pFile,"\n   MET  Calo:%g PFMet:%g TCMet:%g \n",met_et[0],met_et[3],met_et[4]);
				fclose (pFile);	
		}
 	    
 	    
 	}       
}

void ACAna::FillHistTruth(Int_t hcut,std::vector<int> goodMuons,double mZ = -1.){
   
   h1_2mu_muon_eta_truth[hcut]->Fill(truthl_eta[goodMuons[0]]);
   h1_2mu_muon_eta_truth[hcut]->Fill(truthl_eta[goodMuons[1]]);
   
   if (goodMuons.size()>=2){
		TLorentzVector mm1;
		TLorentzVector mm2;
		mm1.SetPtEtaPhiM(truthl_pt[goodMuons[0]],truthl_eta[goodMuons[0]],truthl_phi[goodMuons[0]],muonmass);
		mm2.SetPtEtaPhiM(truthl_pt[goodMuons[1]],truthl_eta[goodMuons[1]],truthl_phi[goodMuons[1]],muonmass);
		TLorentzVector Z=(mm1)+(mm2);
	   
	   if (mZ<1.){
	    h1_2mu_mass_truth[hcut]->Fill(Z.M());
	    h1_2mu_mass_truth_weighted[hcut]->Fill(Z.M(),weight);
	}else {
	    h1_2mu_mass_truth[hcut]->Fill(mZ);
	    h1_2mu_mass_truth_weighted[hcut]->Fill(mZ,weight);
	 }  
	    
 	    
 	}       
}


// helper
Double_t ACAna::DeltaPhi(double a, double b) {

  double temp = fabs(a-b);
  if (temp <= TMath::Pi())
    return temp;
  else
    return  2.*TMath::Pi() - temp;
}


Double_t ACAna::mT(double et1, double phi1, double et2, double phi2) {

  double mm = 2 * et1 * et2 * ( 1. - cos(phi1 - phi2) );
  return sqrt(mm);

}

void ACAna::readCut(Double_t* y_data, int ndata_max){
    string x_string[ndata_max];
    string x_file;		// x value found on file
    Double_t y_file;		// y value found on file
    int ndata = 0;
    // define input stream "in", open input file
    ifstream in("cuts.txt");
    if (!in)
    {
        cout << "No File 'cuts.txt' in this folder!"<<endl;
    }
    // read all data and store numbers in arrays
    while(1)
    {
        in >> x_file >> y_file;
        if (in.eof()) break;
        if(ndata+1 <= ndata_max)
        {
            x_string[ndata]=x_file;
            y_data[ndata] = y_file;
            cout<<"eingelesen:  "<<x_file<<y_file<<endl;
        }
        ndata++;
    }
    // close input stream
    in.close();
}

Double_t ACAna::MinDEt(const std::vector<TLorentzVector> & objects, 
			 std::vector<UInt_t> * lista, 
			 std::vector<UInt_t> * listb) {
  
   //~ //Find the combination with the lowest DEt
  //~ UInt_t n = objects.size();
  //~ if (n==0) return 0.;
  //~ if (n==1) return objects[0].Et();
  //~ if (n>10)  {
    //~ cout << "MinDEt: n too big : " << n << endl;
    //~ return -1;
  //~ }
  //~ if (lista!=0 && listb!=0) { lista->clear(); listb->clear(); }
  //~ 
  //~ double mindiff = 1000000000., diff = 0.;
  //~ 
  //~ // Determine the combination that minimises the difference
  //~ std::vector< std::vector<UInt_t> > combinationset1;
  //~ std::vector< std::vector<UInt_t> > combinationset2;
  //~ Combinations::mycombinations(n, combinationset1, combinationset2);
  //~ 
  //~ if (combinationset1.size() != combinationset2.size() ) {
    //~ cout << "MinDEt: Combination set sizes to not match - something has gone wrong..." << endl;
  //~ }
  //~ 
  //~ for (UInt_t set = 0; set<combinationset1.size(); set++) {
    //~ 
    //~ std::vector<UInt_t> la = combinationset1[set]; //!< Temporary list a for calculating best combo
    //~ std::vector<UInt_t> lb = combinationset2[set]; //!< Temporary list b for calculating best combo
    //~ 
    //~ Double_t aEt = 0., bEt = 0.;
    //~ for (std::vector<UInt_t>::iterator ia=la.begin();ia!=la.end();++ia) {
            //~ cout << (*ia) << " ";
      //~ aEt += objects[ (*ia) ].Et();
    //~ }
        //~ cout << ", ";
    //~ for (std::vector<UInt_t>::iterator ib=lb.begin();ib!=lb.end();++ib) {
      //~ bEt += objects[ (*ib) ].Et();
            //~ cout << (*ib) << " ";
    //~ }
        //~ cout << endl;
    //~ diff = fabs(aEt - bEt);
        //~ cout << "Difference in Et is " << diff << endl;
    //~ if (diff < mindiff) {
      //~ mindiff = diff;
      //~ if (lista!=0 && listb!=0) { *lista = la; *listb = lb; }
    //~ }
    //~ la.clear(); lb.clear();
    //~ 
  //~ } // end of loop over combination sets
  //~ //
    //~ cout << "Minimum difference is " << mindiff << endl << endl << endl;
    //~ cout << "===========================================" << endl;
  //~ //
  //~ return mindiff;
  return 42.;
}
Double_t ACAna::AlphaT(const std::vector<TLorentzVector> & objects) {

  return 0.5*((SumET(objects) - MinDEt(objects))/(MT(objects)));

}

Double_t ACAna::SumET(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0;    
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et(); 
  }
  return sEt;

}
Double_t ACAna::MT(Double_t*Et, Double_t* px, Double_t* py) {
    Double_t sEt = 0, sPx = 0, sPy = 0;
    for (UInt_t i=0;i<sizeof(Et);i++) { 
        sEt+=Et[i];sPx+=px[i];sPy+=py[i]; 
    }
    Double_t MTsq = sEt*sEt - sPx*sPx - sPy*sPy;
    return MTsq >= 0. ? sqrt(MTsq) : -sqrt(-MTsq);

}

Double_t ACAna::MT(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0, sPx = 0, sPy = 0;
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et();sPx+=o->Px();sPy+=o->Py(); 
  }
  Double_t MTsq = sEt*sEt - sPx*sPx - sPy*sPy;
  return MTsq >= 0. ? sqrt(MTsq) : -sqrt(-MTsq);

}

Double_t ACAna::Minv(TLorentzVector *particle[],Int_t Number) {
  TLorentzVector *sparticle = new TLorentzVector(0,0,0,0);
  for (Int_t i=0;i<Number;++i) { 
      //cout<<"Energie: "<< particle[i]->E()<<" i: "<<i<<endl;
    *sparticle += *(particle[i]); 
  }
  Double_t Minv =sparticle->M(); 
  delete sparticle;
  return Minv;
  //Double_t Minvsq = sE*sE - sPx*sPx - sPy*sPy- sPz*sPz;
  //return Minvsq >= 0. ? sqrt(Minvsq) : -sqrt(-Minvsq);

}

Double_t ACAna::Minv(Double_t E, Double_t E1, Double_t px,Double_t px1, Double_t py,Double_t py1,  Double_t pz, Double_t pz1) {
    Double_t sE = E+E1, sPx = px+px1, sPy = py+py1, sPz=pz+pz1;
//     cout<<sizeof(E)<<endl;
    Double_t Minvsq = sE*sE - sPx*sPx - sPy*sPy - sPz*sPz;
    return Minvsq >= 0. ? sqrt(Minvsq) : -sqrt(-Minvsq);

}

Double_t ACAna::Minv(Double_t *E, Double_t *px,Double_t *py, Double_t *pz) {
//     cout<<sizeof(E)<<endl;
    Double_t sE = 0, sPx = 0, sPy = 0, sPz = 0;
    for (UInt_t i=0;i<sizeof(E);i++) { 
        sE+=E[i];sPx+=px[i];sPy+=py[i];sPz+=pz[i];
    }
    Double_t Minvsq = sE*sE - sPx*sPx - sPy*sPy - sPz*sPz;
    return Minvsq >= 0. ? sqrt(Minvsq) : -sqrt(-Minvsq);

}

Double_t ACAna::MT(TLorentzVector *particle[], Int_t Number) {
  TLorentzVector *sparticle = new TLorentzVector(0,0,0,0);
  for (Int_t i=0;i<Number;++i) { 
    *sparticle += *(particle[i]); 
  }
  Double_t Mt=sparticle->Mt();
  delete sparticle;
  return Mt;
  //Double_t Minvsq = sE*sE - sPx*sPx - sPy*sPy- sPz*sPz;
  //return Minvsq >= 0. ? sqrt(Minvsq) : -sqrt(-Minvsq);

}


bool ACAna::opositeSign(Int_t imuon,Int_t jmuon){
    if(muo_charge[imuon]*muo_charge[jmuon]<0)
        return true;
    return false;
}


bool ACAna::MatchMuons(Int_t index, TLorentzVector *kmuon ){
    TString triggerTemp;
    bool matchedMuons=false;
    //cout<<"haha "<<muo_trign[index]<<endl;
    
    //muo_trig[index][muo_trign[index]]
    for(int m=0; m<muo_trign[index];m++){
        
        Int_t pos=muo_trig[index][m];
        
        Double_t deltaR = sqrt(pow(muo_eta[index]-trig_eta[pos],2)+pow(DeltaPhi(muo_phi[index],trig_phi[pos]),2));
        //cout<<"haha "<<deltaR<<endl;

        triggerTemp=unpack(trig_name[pos]);
        
        if((triggerRun.Contains(triggerTemp.Data())) && (deltaR<0.02)){
            matchedMuons=true;
            }
    }
    
    return matchedMuons;
}


void ACAna::MuonCuts(){
    cut[0]=45.; //PT Threshold
    cut[1]=2.4; // ETA Threshold
    cut[2]=1.; 	// muo_ID[0][1]==1 GLOBAL
    cut[3]=0.1; // muo_TrkIso[i]/muon_pt[i]<0.1
    
    cut[4]=5; 	// muo_TrackerLayersMeasCm[i]=>cutlayer
    
    //~ cut[4]=9; 	// muo_TrackerLayersMeasCm[i]=>cutlayer
    
    cut[5]=2; // muo_ChambersMatched[i]>=2
    cut[6]=1; // muo_ValidMuonHitsCm[i]>=1
    cut[7]=1; // muo_ValidPixelHitsCm[i]>=1
    cut[8]=0.2;// fabs(muo_d0Tk[i])<0.2
	    
    
    for(Int_t jmuon=0;jmuon<muo_n;jmuon++){
       
			if(muo_pt[jmuon]>45.){ //cut2
				  cut_mu_ptrans[jmuon]=true;
				} // end cut2

				if((Double_t)fabs(muo_eta[jmuon])<cut[1]){ //cut3
				  cut_mu_eta[jmuon]=true;
				} // end cut3

				if(muo_ID[jmuon][1]==1){ //cut4  //CHANGED
				  cut_mu_id1[jmuon]=true;
				} // end cut4

				if(muo_TrkIso[jmuon]/muo_pt[jmuon]<0.1){
				  cut_muo_trkiso_zu_pt[jmuon]=true;
				} // end cut5

				if(muo_TrackerLayersMeasTk[jmuon]>5){ //cut6 //CHANGED //CHCHANGED //CUTLAYER
				  cut_muo_TrackerLayersMeasTk[jmuon]=true;
				} // end cut6
				
				if(muo_TrackerLayersMeasTk[jmuon]>8){ //cut6 //CHANGED //CHCHANGED //CUTLAYER
				  cut_muo_TrackerLayersMeasTkTight[jmuon]=true;
				}				

			 if (muo_StationsMatched[jmuon]>=2){ //cut7 //CHANGED
				  cut_muo_StationsMatched[jmuon]=true;
				}//end cut 7

				if (muo_ValidMuonHitsCm[jmuon]>=1){//cut 8 //CHANGED
				  cut_muo_ValidMuonHitsCm[jmuon]=true;
				}//end cut 8
				
			if (muo_ValidPixelHitsCm[jmuon]>=1){ //cut 9 //CHANGED
				  cut_muo_ValidPixelHitsCm[jmuon]=true;
				}//end cut 9
			bool matched=false;	
			for(int m=0; m<muo_trign[jmuon];m++){
				  Int_t pos=muo_trig[jmuon][m];
				  Double_t deltaR = sqrt(pow(muo_eta[jmuon]-trig_eta[pos],2)+pow(DeltaPhi(muo_phi[jmuon],trig_phi[pos]),2));
				  trigname=unpack(trig_name[pos]);
					if(( trigname.Contains("HLT_Mu24") || trigname.Contains("HLT_Mu30") || trigname.Contains("HLT_Mu40")) && (deltaR<0.2)){
					  matched=1;
					}
				}
				if (matched==1){//begin cut 10 //CHANGED
					cut_trigger[jmuon]=true;
				}//end cut 10

			 if (fabs(muo_d0Tk[jmuon])<0.02){//begin cut 11
					 cut_muo_d0[jmuon]=true;
			 }//cut 11
			//     //cout<<muo_d0Tk[jmuon]<<endl;

			if (fabs(muo_dzTk[jmuon])<0.5){//begin cut 11 //eigentlich muo_dzTk ; falsche bez im vorl skimmer
					cut_muo_dzTk[jmuon]=true;
			}//cut 11 
        for(Int_t iimuon=0;iimuon<muo_n;iimuon++){
            if(iimuon!=jmuon)
                if(acos(-(muo_px[iimuon]*muo_px[jmuon]+muo_py[iimuon]*muo_py[jmuon]+muo_pz[iimuon]*muo_pz[jmuon])/(muo_p[iimuon]*muo_p[jmuon])) < 0.02)
                cutCosmicAngle[jmuon]=false;
        }
        
    }
    
}

void ACAna::Muon2Cuts(std::vector<int> goodMuons){
	int muon0=goodMuons[0];
	int muon1=goodMuons[1];
	if((fabs(muo_eta[muon0]) < 2.1 ) || (fabs(muo_eta[muon1]) < 2.1 )){
			cut_eta_2mu=true;
	}
	if((cut_trigger[muon0]) || (cut_trigger[muon1])){
			cut_trigger_2mu=true;
	}
	if((cutCosmicAngle[muon0]) && (cutCosmicAngle[muon1])){
			cut_cosmic_2mu=true;
	}	
	if(muo_charge[muon0]*muo_charge[muon1]<0.){
			cut_charge_2mu=true;
	}		
	
}	
void ACAna::doCounter(bool cutbool[], int hcut2){
    //~ if(cutbool[hcut2]){
        //~ h1_counter[4] ->Fill(hcut2);
        //~ cutbool[hcut2]=false;
    //~ }
}
