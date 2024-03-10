//external includes (root, fastjet, pythia)
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TInterpreter.h"

//standard cpp libraries
#include <vector>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

//namespaces
using namespace Pythia8;
using namespace fastjet;


int main(int argc, char *argv[]) {

//need to generate these dictionaries to output a vector of vectors...	
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

//get input job number (used for labeling outputs in batch submission if needed)
  if(argc != 2)
  {
    std::cout << "Usage: ./llHigherEnergies <job number>" << std::endl;
    return 1;
  }  
  int jobNumber = std::atoi(argv[1]);

  //cut on the minimum number of charged constituents in jets before the jet is output (0 is no cut)
  const float NchCut = 0;

//*********************************************** PYTHIA SETTINGS **************************************************
  Pythia pythia;
  Event& event = pythia.event;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  //pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("WeakDoubleBoson:ffbar2gmZgmZ = on");
  pythia.readString("WeakDoubleBoson:ffbar2WW = on");
  pythia.readString("HiggsSM:ffbar2HZ = on");

  //beam intialization (11 for electrons, 13 for muons)
  //pythia.readString("Beams:idA =  11");
  //pythia.readString("Beams:idB = -11");
  //pythia.readString("Beams:idA =  13");
  //pythia.readString("Beams:idB = -13");
  
  //center of mass energy (should be >Z mass for this setup, use a different setup for Z pole events to force hadronic decays)
  //pythia.settings.parm("Beams:eCM", 163.0);//WW
  //pythia.settings.parm("Beams:eCM", 240.0);//ZHiggs
  //pythia.settings.parm("Beams:eCM", 365.0);//ttbar
  //pythia.settings.parm("Beams:eCM", 3000.0);//muon scenario 1
  //pythia.settings.parm("Beams:eCM", 10000.0);//muon scenario 2

  //Set the RNG seed to be based on the system clock (seed = 0 means use clock)
  //Should be set to another value if you want reproducable MC events
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
 
  //initialize
  pythia.init();

  
  //*********************************************** ROOT and Tree SETUP  **************************************************
  clock_t now = clock();

  //basic stuff
  float NchHadrons;  
  float weight;
  int processCode;
  int nMPI;
  std::vector< float > px;
  std::vector< float > py;
  std::vector< float > pz;
  std::vector< float > m;
  std::vector< int > pid;
  std::vector< int > chg;
  //DIS stuff
  float Q2,W2,x,y;  
  //jet stuff
  std::vector< float > genJetPt;
  std::vector< float > genJetEta;
  std::vector< float > genJetPhi;
  std::vector< int >   genJetChargedMultiplicity;
  std::vector< std::vector<int> > gendau_chg;
  std::vector< std::vector<int> > gendau_pid;
  std::vector< std::vector<float> > gendau_pt;
  std::vector< std::vector<float> > gendau_eta;
  std::vector< std::vector<float> > gendau_phi;

  TFile * f = TFile::Open(Form("llHigherEnergies_%d.root",jobNumber),"recreate");
  TTree * trackTree = new TTree("trackTree","v1");
 
  trackTree->Branch("NchHadrons",&NchHadrons);
  trackTree->Branch("weight",&weight);
  trackTree->Branch("processCode",&processCode);
  trackTree->Branch("nMPI",&nMPI);
  trackTree->Branch("px",&px);
  trackTree->Branch("py",&py);
  trackTree->Branch("pz",&pz);
  trackTree->Branch("m",&m);
  trackTree->Branch("pid",&pid);
  trackTree->Branch("chg",&chg);
  
  trackTree->Branch("Q2",&Q2);
  trackTree->Branch("W2",&W2);
  trackTree->Branch("x",&x);
  trackTree->Branch("y",&y);
  
  trackTree->Branch("genJetEta",&genJetEta);
  trackTree->Branch("genJetPt",&genJetPt);
  trackTree->Branch("genJetPhi",&genJetPhi);
  trackTree->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);
  trackTree->Branch("genDau_chg",		&gendau_chg); 
  trackTree->Branch("genDau_pid",		&gendau_pid);	 
  trackTree->Branch("genDau_pt",		&gendau_pt);
  trackTree->Branch("genDau_eta",		&gendau_eta);	 
  trackTree->Branch("genDau_phi",		&gendau_phi );

  //*************************************** FASTJET SETUP ********************************************
  //jet clustering
  // choose a jet definition
  double R = 0.8;
  JetDefinition jet_def(antikt_algorithm, R);

  //******************************************** ANALYZER ********************************************************
  // Begin event loop.
  int nEvent = 10000;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if( iEvent%1000 == 0 ) std::cout << iEvent << std::endl;
    if (!pythia.next()) continue;
    
    //dumping basic particle info
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal()){
	weight = pythia.info.weight();
	processCode = pythia.info.code();
	nMPI = pythia.info.nMPI();
        if(pythia.event[i].isCharged() && pythia.event[i].isHadron()) NchHadrons++;
        px.push_back(pythia.event[i].px());
        py.push_back(pythia.event[i].py());
        pz.push_back(pythia.event[i].pz());
        m.push_back(pythia.event[i].m());
        pid.push_back(pythia.event[i].id());
        chg.push_back(pythia.event[i].charge());
      }

    //getting pythia particle info and setting up pseudojets for clustering
    vector<PseudoJet> particles;
    int multiplicity = 0;
    for (int i = 1; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal()){
       
	//basic kinematic cuts
        //if( TMath::Abs( pythia.event[i].eta() ) > 2.4 || pythia.event[i].pT() < 0.3) continue;
        PseudoJet pj = PseudoJet(   pythia.event[i].px(),  pythia.event[i].py(),  pythia.event[i].pz(), pythia.event[i].e() );
        pj.set_user_index(i);
        particles.push_back( pj );
        if( pythia.event[i].isCharged() ) multiplicity++;
      } 
    }
    //used for filtering on multiplicity if desired 
    if(multiplicity<NchCut) continue;

    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    //loop over jets
    for (unsigned i = 0; i < jets.size(); i++) {
      //jet cut	    
      //if(jets[i].pt() < 500) continue;

      //loop over constituents
      std::vector< float > tmp_pt;
      std::vector< float > tmp_eta;
      std::vector< float > tmp_phi;
      std::vector< int > tmp_chg;
      std::vector< int > tmp_pid;
      vector<PseudoJet> constituents = jets[i].constituents();
      int chMult = 0;
      for (unsigned j = 0; j < constituents.size(); j++) {
        if( ( pythia.event[ constituents[j].user_index() ] ).isCharged()){
          chMult++;
        }
        tmp_pt.push_back(constituents[j].pt());
        tmp_eta.push_back(constituents[j].eta());
        tmp_phi.push_back(constituents[j].phi());
        tmp_chg.push_back( (pythia.event[constituents[j].user_index()]).charge() );
        tmp_pid.push_back( (pythia.event[constituents[j].user_index()]).id() );
      }
 
      //used for filtering on multiplicity if desired
      if(chMult<NchCut) continue;

      //filling output branches
      genJetPt.push_back(jets[i].pt());
      genJetEta.push_back(jets[i].eta());
      genJetPhi.push_back(jets[i].phi());
      genJetChargedMultiplicity.push_back(chMult);
      gendau_pt.push_back(tmp_pt);
      gendau_eta.push_back(tmp_eta);
      gendau_phi.push_back(tmp_phi);
      gendau_chg.push_back(tmp_chg);
      gendau_pid.push_back(tmp_pid);
    }

    //fill and clear branches
    trackTree->Fill();
    NchHadrons = 0;
    weight = 0;
    processCode = 0;
    nMPI = 0;
    px.clear();
    py.clear();
    pz.clear();
    m.clear();
    pid.clear();
    chg.clear();
    x=0;
    y=0;
    Q2=0;
    W2=0;
    genJetPt.clear();
    genJetEta.clear();
    genJetPhi.clear();
    genJetChargedMultiplicity.clear();
    gendau_chg.clear();
    gendau_pid.clear();
    gendau_pt.clear();
    gendau_eta.clear();
    gendau_phi.clear();
  }

  //end of job stuff
  pythia.stat();
  std::cout << ((float)(clock() - now))/CLOCKS_PER_SEC << " seconds" << std::endl;

  trackTree->Write();
  f->Close();
  return 0;
}
