// =====================================================================
//  Pi0EtaEfficiency.cc
//  sPHENIX Run-24 p+p  pi0/eta -> gamma gamma  reconstruction EFFICIENCY
// ---------------------------------------------------------------------
//  WHAT THIS MODULE DOES (read this first):
//
//  The efficiency is  eff(pT) = N_reco(pT) / N_truth(pT)  for a meson that
//  decays to two photons. We get it in TWO passes:
//
//   PASS 1 -- Fun4All SubsysReco pass (runs inside the G4/reco chain):
//     * Init()           : open the output ROOT file, book trees + histograms.
//     * process_event()  : per event, read the vertex, (for Au+Au) centrality,
//                          EMCal CLUSTERS (reco photons) and the G4 TRUTH
//                          particles. It fills:
//                            - _eventTree     : one entry/event of reco clusters
//                            - _truthEventTree: the truth pi0/eta whose BOTH
//                                               decay photons land in acceptance
//                                               (this is the efficiency DENOMINATOR)
//     * End()            : write the file.
//
//   PASS 2 -- standalone matching pass (called from macro/run_efficiency.C):
//     * recoTruthMatch() : reads the two trees back, pairs reco clusters into
//                          diphotons, position-matches them to the truth
//                          photons, applies the analysis cuts, and fills the
//                          NUMERATOR (h_recoPtVsCentrality) and DENOMINATOR
//                          (h_truthPtVsCentrality). Their ratio = efficiency.
//
//  The Loop()/Loop_background_event_mixing() methods are the alternative
//  invariant-mass (foreground + mixed-event background) path; not used by the
//  default efficiency driver but kept for mass-spectrum studies.
//
//  COLLISION SYSTEM (important for this space): this is PYTHIA p+p, which has
//  NO centrality. The m_isPP flag (default true) makes the centrality logic
//  degenerate to a single inclusive bin and guards the centrality node. See
//  setIsPP() and getCentralityBin().
// =====================================================================

//truth-reco matching
#include "eventTree.h"
#include "truthEventTree.h"

#include "Pi0EtaEfficiency.h"

// Utility
#include <vector>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF2.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TLorentzVector.h>
#include <algorithm>  // for max, max_element
#include <cmath>      // for abs
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <TRandom3.h>
#include <random>
#include <TStopwatch.h>
#include <deque>

// Fun4All
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// Event
#include <Event/Event.h>
#include <Event/packet.h>
#include <ffaobjects/EventHeaderv1.h>

//Trigger
#include <calotrigger/TriggerRunInfov1.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <calotrigger/MinimumBiasInfo.h>

// Global vertex
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoSimv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfov3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfov4.h>
#include <calobase/TowerInfoDefs.h>

// MBD
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

// Centrality MB
#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov2.h>

//GL1
#include <ffarawobjects/Gl1Packet.h>

// G4
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>


//____________________________________________________________________________..
Pi0EtaEfficiency::Pi0EtaEfficiency(const std::string &name, const char *outfile)
  : SubsysReco(name)
  , _caloname("CEMC")

{
  outFileName = Form("%s",outfile);  
}

//____________________________________________________________________________..
Pi0EtaEfficiency::~Pi0EtaEfficiency()
{
  delete _eventTree;
}


//____________________________________________________________________________..
int Pi0EtaEfficiency::Init(PHCompositeNode* topNode)
{
  if (topNode == 0)
    {
      std::cout << "printing topNode" << topNode << std::endl;	
    }

  if(isDATA)
    {
      InitOutputFile();
      InitTree();
      InitHistos();
    }
  
  if(isHIJING)
    {
      InitOutputFile();
      InitTree();
      InitHistos();
      InitTruthTree();
      InitTruthHistos();
    }
    
  std::cout << "Done initializing outputfile and ttree..." << std::endl;
  
  return 0;
  
}

void Pi0EtaEfficiency::InitOutputFile()
{
  std::cout << "Writing to output file : " << outFileName.c_str() << std::endl;
  anaOutFile = new TFile(outFileName.c_str(),"RECREATE");
}

void Pi0EtaEfficiency::InitTree()
{
  std::cout << "Initializing TTree..." << std::endl;
  
  _eventTree = new TTree("_eventTree", "An event level info Tree");
  
  _eventTree->Branch("_runnumber",&_runnumber);
  _eventTree->Branch("_nClusters", &_nClusters, "_nClusters/I");
 // _eventTree->Branch("_Centrality", &_Centrality, "_Centrality/F");
  _eventTree->Branch("_vertex", _vertex, "_vertex[3]/F");
  _eventTree->Branch("_towerEnergy", &_towerEnergy, "_towerEnergy[_nClusters]/F");
  _eventTree->Branch("_clusterEnergies", _clusterEnergies, "_clusterEnergies[_nClusters]/F");
  _eventTree->Branch("_clusterPts", _clusterPts, "_clusterPts[_nClusters]/F");
  _eventTree->Branch("_clusterEtas", _clusterEtas, "_clusterEtas[_nClusters]/F");
  _eventTree->Branch("_clusterPhis", _clusterPhis, "_clusterPhis[_nClusters]/F");
  _eventTree->Branch("_maxTowerEtas", _maxTowerEtas, "_maxTowerEtas[_nClusters]/I");
  _eventTree->Branch("_maxTowerPhis", _maxTowerPhis, "_maxTowerPhis[_nClusters]/I");
  _eventTree->Branch("_clusterChi2", _clusterChi2, "_clusterChi2[_nClusters]/F");
  _eventTree->Branch("_photonProb",&_photonProb,"_photonProb[_nClusters]/F");
  _eventTree->Branch("ScaledTriggerBit",ScaledTriggerBit,"ScaledTriggerBit[64]/O");
  _eventTree->Branch("LiveTriggerBit",LiveTriggerBit,"LiveTriggerBit[64]/O");
  _eventTree->Branch("Scaledowns",m_scaledowns,"Scaledowns[64]/I");
  _eventTree->Branch("count_raw",count_raw,"count_raw[64]/L");
  _eventTree->Branch("count_live",count_live,"count_live[64]/L");
  _eventTree->Branch("count_scaled",count_scaled,"count_scaled[64]/L");
      
  
}//end InitTree

void Pi0EtaEfficiency::InitTruthTree()
{
  _truthEventTree = new TTree("_truthEventTree","");
  _truthEventTree->Branch("_truthZvtx", _truthZvtx, "_truthZvtx[3]/F");
  //_truthEventTree->Branch("_npart", &_npart, "_npart/I");
  //_truthEventTree->Branch("_ncoll", &_ncoll, "_ncoll/I");
 // _truthEventTree->Branch("_bimp", &_bimp, "_bimp/F");
  _truthEventTree->Branch("_nEtas",&_nEtas,"_nEtas/I");
 // _truthEventTree->Branch("_truthCentrality", &_truthCentrality, "_truthCentrality/I");
  //_truthEventTree->Branch("_centImpactParam", &_centImpactParam, "_centImpactParam/F");
  _truthEventTree->Branch("_etaEmbedStatus",_etaEmbedStatus,"_etaEmbedStatus[_nEtas]/O");
  _truthEventTree->Branch("_NrejectedEtas",&_NrejectedEtas,"_NrejectedEtas/I");
  _truthEventTree->Branch("_allEtas",&_allEtas,"_allEtas/I");

  acceptedEtaVec = new TClonesArray("TLorentzVector");
  acceptedPhoton1Vec = new TClonesArray("TLorentzVector");
  acceptedPhoton2Vec = new TClonesArray("TLorentzVector");
  
  _truthEventTree->Branch("acceptedEtaVec", "TClonesArray", &acceptedEtaVec, 32000, 0);
  _truthEventTree->Branch("acceptedPhoton1Vec", "TClonesArray", &acceptedPhoton1Vec, 32000, 0);
  _truthEventTree->Branch("acceptedPhoton2Vec", "TClonesArray", &acceptedPhoton2Vec, 32000, 0);
  
  rejectedEtaVec = new TClonesArray("TLorentzVector");
  rejectedPhoton1Vec = new TClonesArray("TLorentzVector");
  rejectedPhoton2Vec = new TClonesArray("TLorentzVector");
  
  _truthEventTree->Branch("rejectedEtaVec", "TClonesArray", &rejectedEtaVec, 32000, 0);
  _truthEventTree->Branch("rejectedPhoton1Vec", "TClonesArray", &rejectedPhoton1Vec, 32000, 0);
  _truthEventTree->Branch("rejectedPhoton2Vec", "TClonesArray", &rejectedPhoton2Vec, 32000, 0);
  
  //_truthEventTree->Branch("_photonAcceptance", _photonAcceptance,"_photonAcceptance[_nEtas]/O");
  //_truthEventTree->Branch("_etaIsGood",_etaIsGood, "_etaIsGood[_nEtas]/O");
  //_truthEventTree->Branch("_truthEtaPseudo",_truthEtaPseudo,"_truthEtaPseudo[_nEtas]/F");
  //_truthEventTree->Branch("_truthEtaPt",_truthEtaPt,"_truthEtaPt[_nEtas]/F");
  //_truthEventTree->Branch("_truthEtaE",_truthEtaE,"_truthEtaE[_nEtas]/F");
  //_truthEventTree->Branch("_truthEtaPhi",_truthEtaPhi,"_truthEtaPhi[_nEtas]/F"); 
  //_truthEventTree->Branch("_nGoodPhotons",&_nGoodPhotons,"_nGoodPhotons/I");
  //_truthEventTree->Branch("_nBadPhotons",&_nBadPhotons,"_nBadPhotons/I");
  //_truthEventTree->Branch("_truthBadPhotonPt",_truthBadPhotonPt,"_truthBadPhotonPt[_nBadPhotons]/F");
  //_truthEventTree->Branch("_truthBadPhotonEta",_truthBadPhotonEta,"_truthBadPhotonEta[_nBadPhotons]/F");
  //_truthEventTree->Branch("_truthBadPhotonPhi",_truthBadPhotonPhi,"_truthBadPhotonPhi[_nBadPhotons]/F");
  //_truthEventTree->Branch("_truthGoodPhotonPt",_truthGoodPhotonPt,"_truthGoodPhotonPt[_nGoodPhotons]/F");
  //_truthEventTree->Branch("_truthGoodPhotonEta",_truthGoodPhotonEta,"_truthGoodPhotonEta[_nGoodPhotons]/F");
  //_truthEventTree->Branch("_truthGoodPhotonPhi",_truthGoodPhotonPhi,"_truthGoodPhotonPhi[_nGoodPhotons]/F");

  
}


void Pi0EtaEfficiency::InitHistos()
{
  
  h_diphotonFgMass = new TH1F("h_diphotonFgMass", "Foreground diphoton invariant mass;m_{#gamma#gamma} [GeV/c^{2}];Counts", 100, 0.0, 1.0);

  h_diphotonBgMass = new TH1F("h_diphotonBgMass", "Mixed-event background diphoton mass;m_{#gamma#gamma} [GeV/c^{2}];Counts", 100, 0.0, 1.0);

  h_diphotonFgMassVsPt = new TH2F("h_diphotonFgMassVsPt", "Foreground diphoton mass vs pair p_{T};m_{#gamma#gamma} [GeV/c^{2}];pair p_{T} [GeV]", 100, 0.0, 1.0, 300, 0.0, 30.0);

  h_diphotonBgMassVsPt = new TH2F("h_diphotonBgMassVsPt", "Background diphoton mass vs pair p_{T};m_{#gamma#gamma} [GeV/c^{2}];pair p_{T} [GeV]", 100, 0.0, 1.0, 300, 0.0, 30.0);

  //clusterMap = new TH2F("clusterMap","Cluster Position", 96, -0.5, 95.5, 256, -0.5, 255.5);
  //clusterMap->SetXTitle("Cluster #eta");
  //clusterMap->SetYTitle("Cluster #phi");

  h_clusterEnergy = new TH1F("h_clusterEnergy","Cluster energy;Cluster energy [GeV];Counts / (0.250 GeV)",160,0,40);

  //diPhotonMap = new TH2F("diPhotonMap","diPhoton Map",  96, -1.2, 1.2, 256, -3.2, 3.2);
 // diPhotonMap->SetXTitle("Reco #eta");
 // diPhotonMap->SetYTitle("Reco #phi");

  h_zVertex = new TH1F("h_zVertex","Event z-vertex;z-vertex [cm];Counts",600,-150,150);

  //h_eventCentrality = new TH1F("h_eventCentrality","Event centrality;Centrality class;Counts",100,0,1000);

  h_diphotonFgDeltaR = new TH1F("h_diphotonFgDeltaR","Foreground photon-pair #DeltaR;#DeltaR;Counts",150.0, 0, 1.5);

  h_diphotonBgDeltaR = new TH1F("h_diphotonBgDeltaR","Background photon-pair #DeltaR;#DeltaR;Counts",150.0, 0, 1.5);

  h_diphotonFgAsymmetry = new TH1F("h_diphotonFgAsymmetry","Foreground energy asymmetry;#alpha = |E_{1}-E_{2}|/(E_{1}+E_{2});Counts", 110.0, 0.0, 1.1);

  h_diphotonBgAsymmetry = new TH1F("h_diphotonBgAsymmetry","Background energy asymmetry;#alpha = |E_{1}-E_{2}|/(E_{1}+E_{2});Counts", 110.0, 0.0, 1.1);

  h_towerGoodMap = new TH2F("h_towerGoodMap","Good-tower occupancy map;#eta bin;#phi bin;N_{isGood}",96.0, -0.5, 95.5, 256.0, -0.5, 255.5);

  //ieta_iphi_energy = new TH3F("ieta_iphi_energy","",96, -0.5, 95.5, 256, -0.5, 255.5, 20.0, 0.0, 10.0);
  //ieta_iphi_energy->SetXTitle("#eta bin");
  //ieta_iphi_energy->SetYTitle("#phi bin");
  //ieta_iphi_energy->SetZTitle("Energy [GeV]");

  //bkg_eff = new TH1F("bkg_eff","",1,0,1);
  //bkg_eff->SetXTitle("");
  //bkg_eff->SetYTitle("Efficiency");

  h_usedEventCount = new TH1F("h_usedEventCount","",1,0,1);

  /*
  h_scaledTriggerBit = new TH1I("h_scaledTriggerBit","",63,0,63);
  h_scaledTriggerBit->GetXaxis()->SetNdivisions(64);
  h_scaledTriggerBit->GetXaxis()->SetLabelSize(0.024);

  //for different run # range that have different trigger names for certain indices
  h2_scaledTriggerBit = new TH1I("h2_scaledTriggerBit","",63,0,63);
  h2_scaledTriggerBit->GetXaxis()->SetNdivisions(64);
  h2_scaledTriggerBit->GetXaxis()->SetLabelSize(0.024);
  */

 // h_zVertexVsCentrality = new TH2F("h_zVertexVsCentrality","z-vertex vs centrality;z-vertex [cm];Centrality [%]",60,-30,30,100,0,1000);

  h_towerHotMap = new TH2F("h_towerHotMap","Hot-tower occupancy map;#eta bin;#phi bin",96, -0.5, 95.5, 256, -0.5, 255.5);


  h_diphotonFgDeltaRVsPt = new TH2F("h_diphotonFgDeltaRVsPt", "Foreground #DeltaR vs pair p_{T};#DeltaR;pair p_{T} [GeV]", 150, 0.0, 1.5, 300, 0.0, 30.0);

  h_diphotonBgDeltaRVsPt = new TH2F("h_diphotonBgDeltaRVsPt", "Background #DeltaR vs pair p_{T};#DeltaR;pair p_{T} [GeV]", 150, 0.0, 1.5, 300, 0.0, 30.0);

  h_diphotonFgDeltaRVsPtMassSideband = new TH2F("h_diphotonFgDeltaRVsPtMassSideband","Foreground #DeltaR vs pair p_{T} (mass sideband);#DeltaR;pair p_{T} [GeV]", 150, 0.0, 1.5, 300, 0.0, 30.0);

  h_diphotonFgMassPtDeltaR = new TH3F("h_diphotonFgMassPtDeltaR", "Foreground m_{#gamma#gamma} vs pair p_{T} vs #DeltaR;m_{#gamma#gamma} [GeV/c^{2}];pair p_{T} [GeV];#DeltaR", 100, 0.0, 1.0, 300, 0.0, 30.0, 150, 0.0, 1.5);

  h_diphotonBgMassPtDeltaR = new TH3F("h_diphotonBgMassPtDeltaR", "Background m_{#gamma#gamma} vs pair p_{T} vs #DeltaR;m_{#gamma#gamma} [GeV/c^{2}];pair p_{T} [GeV];#DeltaR", 100, 0.0, 1.0, 300, 0.0, 30.0, 150, 0.0, 1.5);
  h_diphotonBgMassPtDeltaR->Sumw2();

  //this weight is only for the mixed bg, weight is 1/(ncc * mcc * mixEvt)
  h_diphotonBgMassVsPtWeighted = new TH2F("h_diphotonBgMassVsPtWeighted", "Mixing-weighted background mass vs pair p_{T};m_{#gamma#gamma} [GeV/c^{2}];pair p_{T} [GeV]",100, 0.0, 1.0, 300, 0.0, 30.0);
  h_diphotonBgMassVsPtWeighted->Sumw2();

}

void Pi0EtaEfficiency::InitTruthHistos()
{

  h_truthZVertex = new TH1F("h_truthZVertex","Truth z-vertex;truth z-vertex [cm];Counts",600,-150,150);

  //h_truthEventCentrality = new TH1F("h_truthEventCentrality","Truth centrality;truth centrality;Counts",100,0,100);

  // Truth meson pT when BOTH decay photons fall within |eta|<photonEtaMax.
  h_truthMesonPtAccepted = new TH1F("h_truthMesonPtAccepted","Truth meson p_{T} (both photons accepted);p_{T}^{truth} [GeV/c];Counts",30,0,30);

  // Truth meson pT for all 2-gamma decays, regardless of photon acceptance.
  h_truthMesonPtAll = new TH1F("h_truthMesonPtAll","Truth meson p_{T} from #gamma#gamma decays;p_{T}^{truth} [GeV/c];Counts",30,0,30);

  // Meson pT vs pseudorapidity when both photons are within acceptance.
  h_truthMesonPtVsEtaAccepted = new TH2F("h_truthMesonPtVsEtaAccepted","Truth meson #eta vs p_{T} (both photons accepted);#eta^{truth};p_{T}^{truth} [GeV/c]",160, -2.0, 2.0, 30, 0, 30);

  // Meson pT vs pseudorapidity for all photon positions.
  h_truthMesonPtVsEtaAll = new TH2F("h_truthMesonPtVsEtaAll","Truth meson #eta vs p_{T} (all photon positions);#eta^{truth};p_{T}^{truth} [GeV/c]",160, -2.0, 2.0, 30, 0, 30);

  h_truthDecayMode = new TH1F("h_truthDecayMode","Truth decay-mode tally;decay mode;Counts",6,0,6);

  // (zvtx, meson eta, meson pT) for events where a decay photon is out of acceptance.
  h_truthLostPhotonZvtxEtaPt = new TH3F("h_truthLostPhotonZvtxEtaPt","Lost-photon kinematics;truth z-vertex [cm];#eta^{truth};p_{T}^{truth} [GeV/c]",300,-150,150,160,-2,2,30,0,30);

  h_truthMesonPtDist = new TH1F("h_truthMesonPtDist","Truth meson p_{T} distribution;p_{T} [GeV];Counts",30,0,30);

}


//____________________________________________________________________________..
int Pi0EtaEfficiency::process_event(PHCompositeNode *topNode)
{

  if(!topNode)
    {
      std::cout << "Pi0EtaEfficiency::process_event - topnode PHCompositeNode not valid! return -1" << std::endl;
      return -1;
    }

  if (ievent % 10 == 0)
    {
      std::cout << std::endl;
      std::cout << "Event number " << ievent <<" running." << std::endl;
      std::cout << "====================================" << std::endl;
    }
  
  
  ievent++;

  std::string towergeomnode = "TOWERGEOM_" + _caloname;

  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());
  if (!towergeom)
    {
      std::cout << PHWHERE << ": Could not find node " << towergeomnode << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  // ======= global vertex stuff =======  
  if(isDATA)
    {
      GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

      if(!vertexmap)
	{
	  std::cout << PHWHERE << "Could not find node GlobalVertexMap" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

      if(vertexmap)
	{
	  if (!vertexmap->empty())
	    {	      
	      GlobalVertex *vtx = vertexmap->begin()->second;
	      
	      if(!vertexmap->isValid() || !vtx)
		{
		  std::cout << "GlobalVertexMap Node is not empty, but not valid or iter failed" << PHWHERE << std::endl;
		  noVtx++;
		  return Fun4AllReturnCodes::ABORTEVENT;
		}
	      
	      vx = vtx->get_x();
	      vy = vtx->get_y();
	      vz = vtx->get_z();
	    }
	  else
	    {
	      noVtx++;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	}

      if(abs(vz) > myzvtx_cut)
	{      
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      
      _vertex[0] = vx;
      _vertex[1] = vy;
      _vertex[2] = vz;

      h_zVertex->Fill(vz);
    }

  if(isHIJING)
    {
      //mbd determined zvtx
      MbdVertexMap *mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode,"MbdVertexMap");
      
      bool isglbvtx=true;
      
      if(!mbdvtxmap)
	{
	  std::cout << "WARNING! Did not get MDBVERTEXMAP node " << PHWHERE << std::endl;
	  isglbvtx=false;
	  //return Fun4AllReturnCodes::ABORTEVENT;
	}
      
      if(mbdvtxmap->empty())
	{
	  std::cout << "WARNING! MDBVERTEXMAP is empty " << PHWHERE << std::endl;	  
	  isglbvtx=false;	  
	  //return Fun4AllReturnCodes::ABORTEVENT;
	}
     
      if(isglbvtx)
	{
	  MbdVertex *bvertex = nullptr;
	  
	  if (mbdvtxmap)
	    {
	      for (MbdVertexMap::ConstIter mbditer= mbdvtxmap->begin(); mbditer != mbdvtxmap->end(); ++mbditer)
		{
		  bvertex = mbditer->second;
		}
	      
	      if(!bvertex)
		{
		  std::cout << "could not find globalvtxmap iter :: set vtx to (-999,-999,-999)" << std::endl;
		}
	      else if(bvertex)
		{
		  vz = bvertex->get_z();
		  vy = bvertex->get_y();
		  vx = bvertex->get_x();
		}
	    }
	}
      
      
      //truth zvtx
      PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");  
      if(!truthinfo)
	{
	  std::cout << PHWHERE << "Could not find truthinfo node G4TruthInfo" << std::endl;
	}

      if (truthinfo)
	{
	  PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());

	  truth_vx = gvertex->get_x();
	  truth_vy = gvertex->get_y();
	  truth_vz = gvertex->get_z();
	}

      if(isglbvtx == false)
	{
	  if(abs(truth_vz) > myzvtx_cut)
	    {      
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  vx = truth_vx;
	  vy = truth_vy;
	  vz = truth_vz;
	  
	  _vertex[0] = vx;
	  _vertex[1] = vy;
	  _vertex[2] = vz;
	}
 
      
      if(isglbvtx == true)
	{
	  if(abs(vz) > myzvtx_cut)
	    {      
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	  
	  _vertex[0] = vx;
	  _vertex[1] = vy;
	  _vertex[2] = vz;	  
	}

      _truthZvtx[0] = truth_vx;
      _truthZvtx[1] = truth_vy;
      _truthZvtx[2] = truth_vz;
	
      h_zVertex->Fill(vz);
      h_truthZVertex->Fill(truth_vz);
    }
  
  // ====== GL1 Information ======
  if(isDATA)
    {
      Gl1Packet *gl1_packet = findNode::getClass<Gl1Packet>(topNode, "14001");
   
      if(!gl1_packet)
	{
	  std::cout << "EROR - did not find GL1 packet node" << std::endl;
	}
   
      if(gl1_packet)
	{
	  uint64_t gl1_scaledtriggervector = gl1_packet->lValue(0, "ScaledVector");
	  uint64_t gl1_livetriggervector = gl1_packet->lValue(0, "TriggerVector");

	  for (int i = 0; i < 64; i++)
	    {
	      ScaledTriggerBit[i] = ((gl1_scaledtriggervector & 0x1U) == 0x1U);
	      LiveTriggerBit[i] = ((gl1_livetriggervector & 0x1U) == 0x1U);
	      count_raw[i] = gl1_packet->lValue(i, 0);
	      count_live[i] = gl1_packet->lValue(i, 1);    
	      count_scaled[i] = gl1_packet->lValue(i, 2);
	      gl1_scaledtriggervector = (gl1_scaledtriggervector >> 1U) & 0xffffffffU;
	      gl1_livetriggervector = (gl1_livetriggervector >> 1U) & 0xffffffffU;
	    }
	}
    }

  // ======= MinBias stuff =======
  if(isDATA)
    {
      MinimumBiasInfov1 *minBiasInfo = findNode::getClass<MinimumBiasInfov1>(topNode,"MinimumBiasInfo");
  
      if(!minBiasInfo)
	{
	  std::cout << "Error -- did not get minBiasInfo node." << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT; 
	}
  
      bool isMinBias = minBiasInfo->isAuAuMinimumBias();
      if(!isMinBias)
	{
	  evt_isNotMB++;
	  //std::cout << "Not minbias" << std::endl;    
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      else
	{
	  evt_isMB++;
	}     
    }
  
  // ======= Centrality stuff (commented out -- not used for p+p) =======
  /*
  if(isDATA)
    {
      CentralityInfov2 *cent = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");
      if (!cent)
	{
	  std::cout << " ERROR -- can't find CentralityInfo node. " << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;

	}

      int centval = cent->get_centrality_bin(CentralityInfo::PROP::mbd_NS);


      _Centrality = centval;
      h_eventCentrality->Fill(centval);
    }
  */
  
 /* if(isHIJING)
    {
      // -------------------------------------------------------------------
      //  Event-level centrality bookkeeping.
      //
      //  PHYSICS (p+p, this analysis): Npart, Ncoll, impact parameter (bimp)
      //  and "centrality" are Glauber quantities that exist ONLY for
      //  nucleus-nucleus collisions. PYTHIA p+p has a single nucleon per side
      //  -> no impact-parameter geometry, no centrality; every minimum-bias
      //  p+p event is the same class. So in p+p mode we do NOT depend on the
      //  centrality node (it may be absent or hold meaningless values); we
      //  store neutral sentinels and put every event in ONE inclusive bin.
      //  The Au+Au branch (m_isPP == false) keeps the original behavior so
      //  this module can still be reused on HI samples.
      // -------------------------------------------------------------------
      EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
      if(event_header)
	{
	  _npart = event_header->get_intval("npart");
	  _ncoll = event_header->get_intval("ncoll");
	  _bimp = event_header->get_floatval("bimp");
	}
      else if(!m_isPP)
	{
	  // Au+Au needs the EventHeader (npart/ncoll/bimp); p+p does not.
	  std::cout << "Did not get EVENT_HEADER " << PHWHERE << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

      if(m_isPP)
	{
	  // No centrality for p+p: one inclusive class for every event.
	  _truthCentrality = 0;
	  _centImpactParam = _bimp;   // ~0 for p+p; carried only for the record
	  h_truthEventCentrality->Fill(_truthCentrality);
	}
      else
	{
	  // Au+Au: read the real centrality node. Guard against a missing node
	  // so we never dereference a null pointer (the original code did).
	  CentralityInfov1* cent_node = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
	  if (!cent_node)
	    {
	      std::cout << "Error can not find centrality node " << PHWHERE << std::endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	  _truthCentrality = cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
	  _centImpactParam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
	  h_truthEventCentrality->Fill(_truthCentrality);
	}
    }

  */
  
  

  // ======= loop over calibrated towers to fill isGood map =======
  if(isDATA)
    {
      std::string towernode = "TOWERINFO_CALIB_CEMC";

      TowerInfoContainer *tower_container = findNode::getClass<TowerInfoContainer>(topNode, towernode.c_str());

      if (!tower_container)
	{
	  std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}

      unsigned int nchannels = tower_container->size();

      TowerInfo *towerInfo;

      for(unsigned int channel = 0; channel < nchannels; channel++)
	{
	  unsigned int ieta = 999;
	  unsigned int iphi = 999;
	  //float energy = -1.0;

	  towerInfo = tower_container->get_tower_at_channel(channel);

	  unsigned int towerkey = tower_container->encode_key(channel);
   
	  ieta = tower_container->getTowerEtaBin(towerkey);
	  iphi = tower_container->getTowerPhiBin(towerkey);
	  
	  //energy = towerInfo->get_energy();
	  //ieta_iphi_energy->Fill(ieta,iphi,energy);
     
	  bool towerStatus = towerInfo->get_isGood();
	  bool isHot = towerInfo->get_isHot();
      
	  if(towerStatus == true)
	    {
	      h_towerGoodMap->Fill(ieta,iphi);
	    }

	  if(isHot == true)
	    {
	      h_towerHotMap->Fill(ieta,iphi);
	    }
	 

	}
    }
    

  // =========== Cluster Info ============
  // want this to execute always for either data or sims
  
  CLHEP::Hep3Vector vertex(vx, vy, vz);
  
  std::string clusnodename = "CLUSTERINFO_CEMC";
  RawClusterContainer *recal_clusters = findNode::getClass<RawClusterContainer>(topNode, clusnodename.c_str());
   if (!recal_clusters)
    {
      clusnodename = "CLUSTER_POS_COR_CEMC";
      recal_clusters = findNode::getClass<RawClusterContainer>(topNode, clusnodename.c_str());
    }

  if (!recal_clusters)
    {
      std::cout << PHWHERE << "ERROR! Could not find " << clusnodename << " node."  << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  RawClusterContainer::ConstRange clusterRange = recal_clusters->getClusters();
  RawClusterContainer::ConstIterator clusterIter;

  RawCluster *savedClusters[3000];

  int nSavedClusters = 0;

  // ---- Build the "reco photon" candidate list ----------------------------
  // Each EMCal cluster is the experimental proxy for one photon. We keep only
  // clusters that look photon-like by two standard sPHENIX EMCal cuts:
  //   * ecore > 0.6 GeV : core energy above noise (a real shower, not a blip)
  //   * chi2   < 6.0    : shower-shape consistent with an EM shower
  // The survivors (savedClusters[]) are what PASS 2 pairs into diphotons.
  for (clusterIter = clusterRange.first; clusterIter != clusterRange.second; ++clusterIter)
    {
      RawCluster *cluster = clusterIter->second;

      float clusterEcore = cluster->get_ecore();
      float clusterChi2 = cluster->get_chi2();

      if(clusterEcore < 0.6)
	continue;

      if(clusterChi2 > 6.0)
	continue;

      savedClusters[nSavedClusters++] = cluster;
    }

  _nClusters = nSavedClusters; //holds saved number of clusters


  TowerInfoContainer *emc_container = findNode::getClass<TowerInfoContainer>(topNode,"TOWERINFO_CALIB_CEMC");
  if(!emc_container)
    {
      std::cout << "ERROR! Could not get TOWERINFO_CALIB_CEMC node " << PHWHERE << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // looping on the saved clusters savedClusters[]
  for (int iCluster = 0; iCluster < nSavedClusters; iCluster++)
    {
      //energy of cluster
      CLHEP::Hep3Vector clusterEcoreVec = RawClusterUtility::GetECoreVec(*savedClusters[iCluster], vertex);

      //cluster attributes
      float clusterEnergy = clusterEcoreVec.mag();
      float clusterEta = clusterEcoreVec.pseudoRapidity();
      float clusterPt = clusterEcoreVec.perp();
      float clusterPhi = clusterEcoreVec.phi();

      //fill histo
      h_clusterEnergy->Fill(clusterEnergy);

      //apply to ntuple
      _clusterEnergies[iCluster] = clusterEnergy;
      _clusterPts[iCluster] = clusterPt;
      _clusterEtas[iCluster] = clusterEta;
      _clusterPhis[iCluster] = clusterPhi;

      _clusterChi2[iCluster] = savedClusters[iCluster]->get_chi2();
      _photonProb[iCluster] = savedClusters[iCluster]->get_prob();


      //vector to hold all the towers etas, phis, and energy in this cluster
      std::vector<int> toweretas;
      std::vector<int> towerphis;
      std::vector<float> towerenergies;

      TowerInfo *towinfo;      
    
      RawCluster::TowerConstRange towers = savedClusters[iCluster]->get_towers();
      RawCluster::TowerConstIterator toweriter;
    
      // loop over the towers of each cluster
      for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
	{

	  int iphi = RawTowerDefs::decode_index2(toweriter->first);
	  // index2 is phi in CYL
	  int ieta = RawTowerDefs::decode_index1(toweriter->first);  
	  // index1 is eta in CYL	    
	  unsigned int towerkey = iphi + (ieta << 16U);

	  unsigned int towerindex =  emc_container->decode_key(towerkey);
	  
	  towinfo = emc_container->get_tower_at_channel(towerindex);
	    
	  double towerenergy = towinfo->get_energy();

	  // put the eta, phi, energy into corresponding vectors
	  toweretas.push_back(ieta);
	  towerphis.push_back(iphi);
	  towerenergies.push_back(towerenergy);
	  	 			
	}
    
      //get the index with largest energy in energy vector 
      int maxTowerIndex = max_element(towerenergies.begin(), towerenergies.end()) - towerenergies.begin();
      
      //save corresponding eta,phi
      maxTowerEta = toweretas[maxTowerIndex];
      maxTowerPhi = towerphis[maxTowerIndex];
      
      //apply to ntuple
      _maxTowerEtas[iCluster] = maxTowerEta;
      _maxTowerPhis[iCluster] = maxTowerPhi;

      //fill histo w/ lead tower eta/phi
      //clusterMap->Fill(maxTowerEta,maxTowerPhi);

    }//cluster iteration


  _eventTree->Fill();

  if(isHIJING)
    {
      ProcessFillTruthParticle(topNode, vz);
    }

  
  
  return Fun4AllReturnCodes::EVENT_OK;
}


/*
void Pi0EtaEfficiency::ProcessFillTruthParticle(PHCompositeNode *topNode, float truth_vz)
{
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    std::cout << "No truth info node found — skipping." << std::endl;
    return;
  }

  acceptedEtaVec->Clear();
  acceptedPhoton1Vec->Clear();
  acceptedPhoton2Vec->Clear();
  rejectedEtaVec->Clear();
  rejectedPhoton1Vec->Clear();
  rejectedPhoton2Vec->Clear();

  int netas = 0;
  int rej_etas = 0;
  int all_etas = 0;

  PHG4TruthInfoContainer::Range primaries = truthinfo->GetPrimaryParticleRange();

  for (auto iter = primaries.first; iter != primaries.second; ++iter)
  {
    PHG4Particle* particle = iter->second;

    int pid = particle->get_pid();
    int track_id = particle->get_track_id();

    //  Select only embedded pi0s
    if (pid != PDGPID) continue; // pi0 PDG code
    if (truthinfo->isEmbeded(track_id) <= 0) continue;
    if (!truthinfo->is_primary(particle)) continue;
    if (particle->get_parent_id() != 0) continue;

    all_etas++;

    TLorentzVector pi0LV;
    pi0LV.SetPxPyPzE(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());

    double shiftedEta = GetShiftedEta(particle);
    h_truthMesonPtDist->Fill(pi0LV.Pt());

    PHG4TruthInfoContainer::Range secondaries = truthinfo->GetSecondaryParticleRange();
    if (secondaries.first == secondaries.second) continue;

    TLorentzVector photon1, photon2;
    bool photon1_filled = false;
    bool photon2_filled = false;

    int recoveredPhotons = 0;
    int lostPhotons = 0;

    for (auto iter2 = secondaries.first; iter2 != secondaries.second; ++iter2)
    {
      PHG4Particle* daughter = iter2->second;
      if (daughter->get_parent_id() != track_id) continue;
      if (daughter->get_pid() != 22) continue; // photon

      double daughterEta = GetShiftedEta(daughter);
      TLorentzVector gamma;
      gamma.SetPxPyPzE(daughter->get_px(), daughter->get_py(), daughter->get_pz(), daughter->get_e());

      if (fabs(daughterEta) < photonEtaMax)
      {
        recoveredPhotons++;
        if (!photon1_filled) { photon1 = gamma; photon1_filled = true; }
        else if (!photon2_filled) { photon2 = gamma; photon2_filled = true; }
      }
      else
      {
        h_truthLostPhotonZvtxEtaPt->Fill(truth_vz, shiftedEta, pi0LV.Pt());
        lostPhotons++;
        if (!photon1_filled) { photon1 = gamma; photon1_filled = true; }
        else if (!photon2_filled) { photon2 = gamma; photon2_filled = true; }
      }
    }

    int totalPhotons = recoveredPhotons + lostPhotons;
    if (totalPhotons != 2) continue; // only interested in 2γ decays

    if (recoveredPhotons == 2)
    {
      h_truthMesonPtAccepted->Fill(pi0LV.Pt());
      h_truthMesonPtVsEtaAccepted->Fill(shiftedEta, pi0LV.Pt());
      new ((*acceptedEtaVec)[netas]) TLorentzVector(pi0LV);
      new ((*acceptedPhoton1Vec)[netas]) TLorentzVector(photon1);
      new ((*acceptedPhoton2Vec)[netas]) TLorentzVector(photon2);
      _etaEmbedStatus[netas] = true;
      netas++;
    }
    else
    {
      new ((*rejectedEtaVec)[rej_etas]) TLorentzVector(pi0LV);
      new ((*rejectedPhoton1Vec)[rej_etas]) TLorentzVector(photon1);
      new ((*rejectedPhoton2Vec)[rej_etas]) TLorentzVector(photon2);
      rej_etas++;
    }

    h_truthMesonPtAll->Fill(pi0LV.Pt());
    h_truthMesonPtVsEtaAll->Fill(shiftedEta, pi0LV.Pt());
    h_truthDecayMode->Fill(0); // only tracking γγ decays in this setup
  }

  _nEtas = netas;
  _NrejectedEtas = rej_etas;
  _allEtas = all_etas;

  _truthEventTree->Fill();
}


*/



//____________________________________________________________________________..
void Pi0EtaEfficiency::ProcessFillTruthParticle(PHCompositeNode *topNode, float truth_vz)
{
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
  
  if(!truthinfo)
    {
      std::cout << "no truth info node... just skip the whole part.." << std::endl;
      return;
    }

  //need to clear TClonesArrays
  acceptedEtaVec->Clear();
  acceptedPhoton1Vec->Clear();
  acceptedPhoton2Vec->Clear();

  rejectedEtaVec->Clear();
  rejectedPhoton1Vec->Clear();
  rejectedPhoton2Vec->Clear();

  //these eta counters include hijing and embedding
  int netas = 0;   //number of etas that ARE within acceptance from 2photon channel
  int rej_etas = 0;//number of etas that are NOT within acceptance from 2photon channel
  int all_etas = 0; //all etas, regardless of decay mode (for a low level sanity check)
  
 
  // ===== Loop over primaries ===== //
  PHG4Particle *g4particle = nullptr;
  
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
   
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      g4particle = iter->second;

      int primary_id = g4particle->get_pid();

      int eta_track_id = g4particle->get_track_id();
      
      //only grab eta mesons
      if(primary_id != PDGPID)
	continue;

      if (! truthinfo->is_primary(g4particle) )
	continue;

      if(g4particle->get_parent_id() != 0)
	continue;

      // ---- SIGNAL ISOLATION (embed id) ----
      // Keep ONLY the embedded gun signal meson (embed id == m_signalEmbedId,
      // = 2 in embedded production). The pythia8_Detroit underlying event is
      // itself embed_id 1, so the old isEmbeded()>0 test let its soft mesons
      // contaminate the truth spectra and the efficiency denominator. Applying
      // the gate HERE guarantees every downstream truth histogram
      // (h_truthMesonPtDist / *All / *Accepted) and acceptedEtaVec is
      // signal-only -> the raw truth yield is FLAT (the gun is flat in pT).
      if (truthinfo->isEmbeded(eta_track_id) != m_signalEmbedId)
	continue;

      all_etas++;
         	
      double shiftedEta = GetShiftedEta(g4particle);
      
      TLorentzVector etaLV;
      etaLV.SetPxPyPzE(g4particle->get_px(),g4particle->get_py(),g4particle->get_pz(),g4particle->get_e());

      h_truthMesonPtDist->Fill(etaLV.Pt());

      // ===== Loop over secondaries ===== //
      noRange = 0;
      
      PHG4TruthInfoContainer::Range range2 = truthinfo->GetSecondaryParticleRange();      
      if(range2.first == range2.second)
	{
	  std::cout << "Secondary particle range is empty!" << std::endl;
	  
	  noRange++;

	  continue;
	}

      //resets after each eta
      recoveredPhotons = 0;
      not_recoveredPhotons = 0;
      recoveredPi0 = 0;
      recoveredPiPlus = 0;
      recoveredPiMinus = 0;

      TLorentzVector photon1, photon2;
      
      bool photon1_filled = false;
      bool photon2_filled = false;
	  
      for(PHG4TruthInfoContainer::ConstIterator iter2 = range2.first; iter2 != range2.second; ++iter2)
	{
	  	  
	  PHG4Particle *decay = iter2-> second;
	  
	  int parentid = decay->get_parent_id();

	  int truthpid = decay->get_pid();
	  
	  //make sure the decay daughter is associated with the track created by eta
	  if(eta_track_id != parentid)
	    {
	      continue;
	    }
	  
	  double daughterShiftedEta = GetShiftedEta(decay);
	  
	  TLorentzVector daughterLV;
	  daughterLV.SetPxPyPzE(decay->get_px(), decay->get_py(), decay->get_pz(), decay->get_e());
	  
	  //track daughters
	  if(truthpid == 111)//pi0 
	    {
	      recoveredPi0++;
	    }	  
	  else if(truthpid == 211)//pi+
	    {
	      recoveredPiPlus++;
	    }	  
	  else if(truthpid == (-211))//pi-
	    {
	      recoveredPiMinus++;
	    }	    
	  else if(truthpid == 22 && abs(daughterShiftedEta) > photonEtaMax)//lost photon case
	    {
	      h_truthLostPhotonZvtxEtaPt->Fill(truth_vz, shiftedEta, etaLV.Pt());

	      if (!photon1_filled)
		{
		  photon1 = daughterLV;
		  photon1_filled = true;
		}
	      else if (!photon2_filled)
		{
		  photon2 = daughterLV;
		  photon2_filled = true;
		}
	      
	      not_recoveredPhotons++;
	    }	  
	  else if(truthpid == 22 && abs(daughterShiftedEta) < photonEtaMax)//good photon case
	    {
	      if (!photon1_filled)
		{
		  photon1 = daughterLV;
		  photon1_filled = true;
		}
	      else if (!photon2_filled)
		{
		  photon2 = daughterLV;
		  photon2_filled = true;
		}
	      
	      recoveredPhotons++;	      
	    }
	  
	}// end daughter loop
     

      //only analyze eta->2gamma decay channel
      int primePhotons = not_recoveredPhotons + recoveredPhotons;
      
      if(primePhotons == 2)
	{
       
	  //should only execute when both photons fall within acceptance
	  if(recoveredPhotons == 2)
	    {
	      h_truthMesonPtVsEtaAccepted->Fill(shiftedEta, etaLV.Pt());

	      h_truthMesonPtAccepted->Fill(etaLV.Pt());

	      new ((*acceptedEtaVec)[netas])TLorentzVector(etaLV);
	      new ((*acceptedPhoton1Vec)[netas])TLorentzVector(photon1);
	      new ((*acceptedPhoton2Vec)[netas])TLorentzVector(photon2);

	      // Signal isolation is enforced at the top of the primary loop
	      // (isEmbeded == m_signalEmbedId), so every accepted meson here is
	      // guaranteed to be the embedded gun signal.
	      _etaEmbedStatus[netas] = true;

	      netas++;
	    }

	  //keep track of eta->2gamma when only 1 or both photons not within acceptance
	  if( (not_recoveredPhotons == 1 && recoveredPhotons == 1) || not_recoveredPhotons == 2)
	    {
	      std::cout << "eta rejected: either 1 or both photons not within acceptance" << std::endl;

	      new ((*rejectedEtaVec)[rej_etas])TLorentzVector(etaLV);
	      new ((*rejectedPhoton1Vec)[rej_etas])TLorentzVector(photon1);
	      new ((*rejectedPhoton2Vec)[rej_etas])TLorentzVector(photon2);
	      
	      rej_etas++;
	    }
		
	  //this will accept eta meson that decays to 2 photons and does not care if one is not within acceptance
	  if( (not_recoveredPhotons == 1 && recoveredPhotons == 1) || not_recoveredPhotons == 2 || recoveredPhotons == 2)
	    {
	      h_truthMesonPtAll->Fill(etaLV.Pt());
	    }

	}//end looking at only eta->2gamma channel

      //fill decay Modes histo
      if(primePhotons == 2)
	{
	  std::cout << "2gamma" << std::endl;
	  h_truthMesonPtVsEtaAll->Fill(shiftedEta, etaLV.Pt());

	  h_truthDecayMode->Fill(0);
	}
      else if(recoveredPi0 == 3)
	{
	  std::cout<<"3pi0"<<std::endl;
	  h_truthDecayMode->Fill(1);
	}
      else if(recoveredPi0 == 1 && recoveredPiPlus == 1  && recoveredPiMinus == 1)
	{  std::cout << "pi0pi+pi-"<<std::endl;
	  h_truthDecayMode->Fill(2);
	}
      else if(recoveredPiPlus == 1 && recoveredPiMinus == 1 && (recoveredPhotons == 1 || not_recoveredPhotons == 1) )
	{
	  std::cout<<"pi+pi-gamma"<<std::endl;
	  h_truthDecayMode->Fill(3);
	}
      else if(noRange == 1)
	{
	  std::cout << "no range in primary particles" << std::endl;
	  h_truthDecayMode->Fill(4);
	}
      else
	{
	  std::cout<<"other"<<std::endl;
	  h_truthDecayMode->Fill(5);
	}
	    
    }//end primary particle range loop

  _nEtas = netas;
  _NrejectedEtas = rej_etas;
  _allEtas = all_etas;
  
  _truthEventTree->Fill();
  
}//end ProcessFillTruthParticle




//____________________________________________________________________________..
double Pi0EtaEfficiency::GetShiftedEta(PHG4Particle *particle)
{
  float px = particle -> get_px();
  float py = particle -> get_py();
  float pz = particle -> get_pz();
  float p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

  return 0.5*log((p+pz)/(p-pz));

}

//____________________________________________________________________________..
int Pi0EtaEfficiency::End(PHCompositeNode *topNode)
{

  std::cout << "Fraction of events that are MB: " << (double)evt_isMB/ievent << std::endl;
  std::cout << "Fraction of events that are not MB: " << (double)evt_isNotMB/ievent << std::endl;
  std::cout << "Fraction of events w/ no vtx: " << (double)noVtx/ievent<< std::endl;

  if (anaOutFile)
    {
      anaOutFile->cd();
     
      anaOutFile->Write();
      
      anaOutFile->Close();
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//really just used for looping bc dont need to rewrite out ttree and no need to make flags to handle this
void Pi0EtaEfficiency::End()
{
  if (anaOutFile)
    {
      std::cout << "Writing output file..." << std::endl;
      anaOutFile->cd();
      
      anaOutFile->Write();
      
      anaOutFile->Close();     
    }
  
  else
    {
      std::cout<<"Error writing out file"<<std::endl;
    }

}
#if 0

void Pi0EtaEfficiency::Loop(TTree *eventTree, float leadPhotonPtCut, float subPhotonPtCut, float asymmetryCut, float deltaRCut, float clusterChi2Cut, float zVertexCut, int minClusterCut, int maxClusterCut, int firstEntry, int lastEntry)
{

  std::cout << "============================================" << std::endl;
  std::cout << "Running in Loop mode to make foreground" << std::endl;

  TTree *inputTree = eventTree;
  if(!inputTree)
    {
      std::cout << "Error assigning tree" << std::endl;
      exit(-1);
    }

  // Set Branches
  //inputTree->SetBranchAddress("_runnumber",&_runnumber);
  inputTree->SetBranchAddress("_nClusters", &_nClusters);
  inputTree->SetBranchAddress("_Centrality", &_Centrality);
  inputTree->SetBranchAddress("_vertex", _vertex);
  inputTree->SetBranchAddress("_towerEnergy", &_towerEnergy);
  inputTree->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  inputTree->SetBranchAddress("_clusterPts", _clusterPts);
  inputTree->SetBranchAddress("_clusterEtas", _clusterEtas);
  inputTree->SetBranchAddress("_clusterPhis", _clusterPhis);
  inputTree->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  inputTree->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  inputTree->SetBranchAddress("_clusterChi2", _clusterChi2);
  inputTree->SetBranchAddress("_photonProb",&_photonProb);
  //inputTree->SetBranchAddress("ScaledTriggerBit",ScaledTriggerBit);

  std::cout << "Total events in eventTree: " << inputTree->GetEntries() << std::endl;

  //keep count of number of used events in given centrality range for eventual scaling in invariant cross seciton
  int nUsedEvents = 0;

  for (int i = firstEntry; i <= lastEntry; i++) // for each event
    {

      inputTree->GetEntry(i);

      /*
      if(_runnumber <= 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				 ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				 ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	{
	  continue;
	}
      //while we could ommit this extra if statement below, due to the trigger indices being the exact same as above, they have different meaning so just want to              separate this out now
      if(_runnumber > 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	{
	  continue;
	}
      */
      
      if(fabs(_vertex[2]) > zVertexCut)
      continue;

      if(_nClusters < minClusterCut  || _nClusters > maxClusterCut)
	continue;

      //if(_nClusters < 2)
      //	continue;

      h_zVertex->Fill(_vertex[2]);

      h_eventCentrality->Fill(_nClusters);

      h_zVertexVsCentrality->Fill(_vertex[2],_nClusters);

      int zIdx = 0;
      int cent_ClusIdx = 0;

      if(makeBuff == true)
	{
	  zIdx = static_cast<int>(floor((_vertex[2] + zVertexCut)));
	  cent_ClusIdx = _nClusters - minClusterCut ;
	}

       nUsedEvents++;

      //fill scaled trigger bit histo
      /*
      for(int i = 0; i < 64; i++)
	{
	  bool val = ScaledTriggerBit[i];
	  if(val == true)
	    {
	      if(_runnumber <=54548)
		{
		  h_scaledTriggerBit->Fill(i);
		}
	      if(_runnumber > 54548)
		{
		  h2_scaledTriggerBit->Fill(i);
		}
	    }
	}
      */
      
      std::vector<TLorentzVector> selectedClusters;
      int nSelectedClusters = 0;
      
      for (int j = 0; j < _nClusters; j++) // for each cluster in event
	{

	  if(_maxTowerEtas[j] >= 92 || _maxTowerEtas[j] <= 3)
	    {
	      continue;
	    }
	  
	  if(_clusterChi2[j] > clusterChi2Cut)
	    continue;

	  float clusterPt  = _clusterPts[j];
	  float clusterEta = _clusterEtas[j];
	  float clusterPhi = _clusterPhis[j];
	  float clusterE   = _clusterEnergies[j];

	  TLorentzVector clusterLV;
	  clusterLV.SetPtEtaPhiE(clusterPt, clusterEta, clusterPhi, clusterE);

	  selectedClusters.push_back(clusterLV);

	  h_clusterEnergy->Fill(clusterLV.E());

	  nSelectedClusters++;

	}//end saving clusters in event

      if(nSelectedClusters < 2)
	continue;

      if(makeBuff == true)
	{
	  if(!selectedClusters.empty())
	    {
	      mixedEventStruct myStruct;
	      myStruct.clusters = selectedClusters;
	      myStruct.eventNumber = i;
      
	      eventBuffer[zIdx][cent_ClusIdx].push_back(myStruct);
	    }
	  else
	    continue;
	}
      
      TLorentzVector photon1, photon2;
		
      for (size_t iLead = 0; iLead < selectedClusters.size(); iLead++) //to iterate on saved clusters
	{

	  photon1 = selectedClusters[iLead];

	  if (photon1.Pt() < leadPhotonPtCut || photon1.Pt() < 0)
	    {
	      continue;
	    }

	  for (size_t iSub = iLead+1; iSub < selectedClusters.size(); iSub++) //to iterate clusters, remove double counting, make foreground
	    {

	      photon2 = selectedClusters[iSub];

	      if (photon2.Pt() < subPhotonPtCut || photon2.Pt() < 0)
		{
		  continue;
		}

	      double asymmetryAlpha = fabs(photon1.E() - photon2.E())/(photon1.E()+ photon2.E());

	      if (asymmetryAlpha > asymmetryCut)
		{
		  continue;
		}

	      if (photon1.DeltaR(photon2) > deltaRCut)
		{
		  continue;
		}

	      h_diphotonFgAsymmetry->Fill(asymmetryAlpha);
	      h_diphotonFgDeltaR->Fill(photon1.DeltaR(photon2));

	      TLorentzVector diphotonLV = photon1 + photon2;

	      float diphotonMass = diphotonLV.M();
	      float diphotonPt = diphotonLV.Perp();
	      float diphotonEta = diphotonLV.Eta();

	      if(diphotonMass < 0.01 || diphotonMass > 1.0)
		continue;


	      if(diphotonLV.Pt() > 1.0 && diphotonEta < 1.0)
		{

		  h_diphotonFgMass->Fill(diphotonMass);
		  h_diphotonFgMassVsPt->Fill(diphotonMass, diphotonPt);
		  //pt1_pt2_invmass_fg->Fill(photon1.Pt(),photon2.Pt(),diphotonMass);

		  h_diphotonFgDeltaRVsPt->Fill(photon1.DeltaR(photon2), diphotonLV.Pt());
		  h_diphotonFgMassPtDeltaR->Fill(diphotonLV.M(), diphotonLV.Pt(), photon1.DeltaR(photon2));

		}

	      //only fill if invmass is outside of pi0 and eta peak - these are 4sigma cuts
	      //pi0 mass(sigma) -- 0.1500(0.0195)
	      //eta mass(sigma) -- 0.614(0.064)
	      if( (diphotonMass < 0.072 || (diphotonMass > 0.228 && diphotonMass < 0.358) ||  diphotonMass > 0.87) )
		{
		  h_diphotonFgDeltaRVsPtMassSideband->Fill(photon1.DeltaR(photon2), diphotonLV.Pt());
		}

	    } // clus2

	} // clus1
      
    }//event loop
  
  h_usedEventCount->SetBinContent(1,nUsedEvents); //make change here

}//end LOOP

#endif

/*
void Pi0EtaEfficiency::Loop_background_event_mixing(TTree *eventTree, float leadPhotonPtCut, float subPhotonPtCut, float asymmetryCut, float deltaRCut, float clusterChi2Cut, float zVertexCut, int nMixedEvents, int minClusterCut, int maxClusterCut)
{

  
  double trck_not_enough_BGE = 0.0;

  std::cout << "=======================================================" << std::endl;
  std::cout << "running in Event Mixing to construct background" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      std::cout << "Error! Cannot get TTree. Exiting." << std::endl;
      exit(-1);
    }
  
  // Set Branches
  t1->SetBranchAddress("_runnumber",&_runnumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_Centrality", &_Centrality);
  t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_towerEnergy", &_towerEnergy);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);
  t1->SetBranchAddress("_photonProb",&_photonProb);
  t1->SetBranchAddress("ScaledTriggerBit",ScaledTriggerBit);


  for (int evt1 = 0; evt1 < t1->GetEntries(); evt1++) // iteration for first EVENT
    {
      if(evt1 % 100 == 0)
	std::cout << "Evt1 : " << evt1 << std::endl;
      
      t1->GetEntry(evt1);
      
      if(_runnumber <= 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				 ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				 ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	{
	  
	  continue;
	}
      //while we could ommit this extra if statement below, due to the trigger indices being the exact same as above, they have different meaning so just want to              separate this out now
      if(_runnumber > 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	{
	  continue;
	}
      
      
      double zvtx_1 = _vertex[2];
      
      if(abs(zvtx_1) > zvtx_cut)
	continue;

      int cent_1 = _Centrality;
      
      if(cent_1 < centlo || cent_1 > centhi)
	continue;

      
      //need to get actual trigger of event and compare to mixed event to ensure we dont mix triggers when making bkg
      int trigger_index1 = -999;
      
      for(int i = 0; i < 64; i++)
	{
	  if(ScaledTriggerBit[i] == 1)
	    trigger_index1 = i;
	}
      
      if(trigger_index1 == -999)
	{
	  std::cout << "Something is wrong.... trigger index is still -999." << std::endl;
	  continue;
	}
      

      TLorentzVector *savClusLV1 =  new TLorentzVector[1000];
      
      int clusIter = 0; //for indexing saved clusters that are good

      //each CLUSTER within first event
      for (int j1 = 0; j1 < _nClusters; j1++)
	{

	  //reject clusters near edges of detector
	  int max_ieta_1 = _maxTowerEtas[j1];
	      
	  if(max_ieta_1 >= 92 || max_ieta_1 <= 3)
	    {
	      continue;
	    }

	  if(_clusterChi2[j1] > chi2_cut)
	    continue;

	  float pt1, eta1, phi1, E1;

	  pt1 = _clusterPts[j1];
	  eta1 = _clusterEtas[j1];
	  phi1 = _clusterPhis[j1];
	  E1 = _clusterEnergies[j1];
	  
	  savClusLV1[j1].SetPtEtaPhiE(pt1, eta1, phi1, E1);

	  clusIter++;
	}

      
      //keep track of  # of events that we used to make bkg
      int passedEvtCut = 0;

      // get 2nd event to mix with
      for (int evt2 = 0; evt2 < t1->GetEntries(); evt2++)
	{
	  	  
	  if(evt1 == evt2)
	    continue; 

	  if(passedEvtCut == bckgnd_evnts)
	    {
	      break;		
	    }
		  
	  if(evt2 == (t1->GetEntries()-1) && passedEvtCut != bckgnd_evnts)
	    {
	      trck_not_enough_BGE++;
	    }

	  t1->GetEntry(evt2);

	  
	  if(_runnumber <= 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				     ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				     ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	    {
	      continue;
	    }
	  //while we could ommit this extra if statement below, due to the trigger indices being the exact same as above, they have different meaning so just want to              separate this out now
	  if(_runnumber > 54548 && (ScaledTriggerBit[10] !=1 && ScaledTriggerBit[16] !=1 && ScaledTriggerBit[17] !=1 &&
				    ScaledTriggerBit[18] !=1 && ScaledTriggerBit[19] !=1 && ScaledTriggerBit[24] !=1 &&
				    ScaledTriggerBit[25] !=1 && ScaledTriggerBit[26] !=1 && ScaledTriggerBit[27] !=1 ) )
	    {
	      continue;
	    }
	  	 
	  int trigger_index2 = -999;
	  
	  for(int i = 0; i < 64; i++)
	    {
	      if(ScaledTriggerBit[i] == 1)
		trigger_index2 = i;
	    }
	  
	  if(trigger_index2 == -999)
	    {
	      std::cout << "Something is wrong ... trigger index 2 still is -999" << std::endl;
	      continue;
	    }
	  
	  if(trigger_index1 != trigger_index2)
	    continue;
	  
	  
	  float zvtx_2 = _vertex[2];
	  
	  //if(abs(zvtx_2) > zvtx_cut)
	  //continue;
	  
	  if (fabs(zvtx_1 - zvtx_2) > 20) 
	    {
	      continue;
	    } 

	  
	  float cent_2 = _Centrality;
	  
	  //if(cent_2 < centlo ||cent_2 > centhi)
	  //continue;
	  
	  if (fabs(cent_1 - cent_2) > 8.0)
	    {
	      continue; 
	    }
	  
	  passedEvtCut++;

	  TLorentzVector *savClusLV2 =  new TLorentzVector[1000];

	  int nClusters2 = _nClusters;

	  //for indexing saved clusters that are good
	  int clusIter2 = 0; 
   
	  for (int j2 = 0; j2 < nClusters2; j2++) //loop CLUSTERS of second event
	    {
	
	      //reject clusters near edges of detector
	      int max_ieta_2 = _maxTowerEtas[j2];
	      
	      if(max_ieta_2 >= 92 || max_ieta_2 <= 3)
		{
		  continue;
		}

	      if(_clusterChi2[j2] > chi2_cut)
		continue;
	
	      float pt2, eta2, phi2, E2;
	      pt2 = _clusterPts[j2];
	      eta2 = _clusterEtas[j2];
	      phi2 = _clusterPhis[j2];
	      E2 = _clusterEnergies[j2];
	      
	      savClusLV2[j2].SetPtEtaPhiE(pt2, eta2, phi2, E2);

	      clusIter2++;

	    }		
	  
	  for (int iCs = 0; iCs < clusIter ; iCs++) // loop over first CLUSTER to generate bg
	    {
	      TLorentzVector pho1 = savClusLV1[iCs];
	
	      if (fabs(pho1.Pt()) < pt1_cut )	
		continue; 
	      
	      for (int jCs = iCs + 1; jCs < clusIter2 ; jCs++) // loop over second CLUSTER
		{

		  TLorentzVector pho2 = savClusLV2[jCs];

		  if (fabs(pho2.Pt()) < pt2_cut) 
		    continue; 

	
		  float mix_alpha = fabs((pho1.E() - pho2.E())) / (pho1.E()+ pho2.E());

		  float deltaR = pho1.DeltaR(pho2);
		    
		  if (mix_alpha > alpha_cut) 
		    continue;	
								
		  if (deltaR > delR_cut) 
		    continue;

		  asym_bg->Fill(mix_alpha);
		  deltaR_bg->Fill(deltaR);
		  
		  TLorentzVector etalv = pho1 + pho2;

		  float pairInvMass = etalv.M();		
		  float pair_pt = etalv.Perp();
		  float pair_eta = etalv.Eta();

		  if(pairInvMass < 0.01 || pairInvMass > 1.0)
		    continue;

		  if(fabs(etalv.Pt()) > 1.0  && pair_eta < 1.0)
		    {
		      pairInvMassTotal_bg->Fill(pairInvMass);
		      eta_bg_invmass_pt->Fill(pairInvMass, pair_pt);

		      pairInvMassPtdelR_bkg->Fill(etalv.M(),etalv.Pt(),deltaR);
		      delR_pairpT_bkg->Fill(deltaR, etalv.Pt());
		    }


		} // inner loop (cluster2, jCs)
	      
	    } // outer loop (cluster1, iCs)

	  delete [] savClusLV2;
	  savClusLV2 = nullptr;
	} //for loop for second event evt2

        delete [] savClusLV1;
	savClusLV1 = nullptr;
  
    }// iterating over all events evt1)

  float bkgEff = (trck_not_enough_BGE/t1->GetEntries())*100.0;

  bkg_eff->Fill(1,bkgEff);

  */
  //--------------------------------------------------------------------------------------
  /*
  
  TTree * t1 = intree;
  if (!intree)
    {
      std::cout << "Error! Cannot get TTree. Exiting." << std::endl;
      exit(-1);
    }
  
  // Set Branches
  //t1->SetBranchAddress("_runnumber",&_runnumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_Centrality", &_Centrality);
  t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_towerEnergy", &_towerEnergy);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  //t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);
  //t1->SetBranchAddress("_photonProb",&_photonProb);
  //t1->SetBranchAddress("ScaledTriggerBit",ScaledTriggerBit);
  
   // Use a deque to store clusters from previous events
  struct MixedEvent
  {
    float zvtx;
    int centrality;
    std::vector<TLorentzVector> clusters;
  };
  
    std::deque<MixedEvent> event_buffer;
    
    for (int evt_current = 0; evt_current < t1->GetEntries(); ++evt_current)
      {
            
        t1->GetEntry(evt_current);

        // -- Cuts for current event --
        if (abs(_vertex[2]) > zvtx_cut)
	  continue;
	
        if (_Centrality < centlo || _Centrality > centhi)
	  continue;

        // Save clusters for current event
        std::vector<TLorentzVector> current_clusters;
	
        for (int j1 = 0; j1 < _nClusters; ++j1)
	  {
            if (_maxTowerEtas[j1] >= 92 || _maxTowerEtas[j1] <= 3)
	      continue;
	    
            if (_clusterChi2[j1] > chi2_cut)
	      continue;

            float pt1 = _clusterPts[j1];
            float eta1 = _clusterEtas[j1];
            float phi1 = _clusterPhis[j1];
            float E1 = _clusterEnergies[j1];
            
            // Only add clusters that pass first cut
            if(pt1 > pt1_cut)
	      {
		TLorentzVector lv;
		lv.SetPtEtaPhiE(pt1, eta1, phi1, E1);
		current_clusters.push_back(lv);
	      }
	  }
        
        // -- Mixing with buffered events --
        for (const auto& buffered_event : event_buffer)
	  {
            // Apply cuts based on event properties
            if (fabs(_vertex[2] - buffered_event.zvtx) > 20)
	      continue;
	    
            if (fabs(_Centrality - buffered_event.centrality) > 8.0)
	      continue;
            
            for (const auto& pho1 : current_clusters)
	      {
                for (const auto& pho2 : buffered_event.clusters)
		  {
                    if (fabs(pho2.Pt()) < pt2_cut)
		      continue;
                    
                    float mix_alpha = fabs((pho1.E() - pho2.E())) / (pho1.E() + pho2.E());
                    float deltaR = pho1.DeltaR(pho2);
                    
                    if (mix_alpha > alpha_cut)
		      continue;
		    
                    if (deltaR > delR_cut)
		      continue;

                    asym_bg->Fill(mix_alpha);
                    deltaR_bg->Fill(deltaR);

                    TLorentzVector etalv = pho1 + pho2;
                    float pairInvMass = etalv.M();
                    float pair_eta = etalv.Eta();
                    
                    if(pairInvMass < 0.01 || pairInvMass > 1.0)
		      continue;
                    
                    if (fabs(etalv.Pt()) > 1.0 && pair_eta < 1.0)
		      {
                        pairInvMassTotal_bg->Fill(pairInvMass);
                        eta_bg_invmass_pt->Fill(pairInvMass, etalv.Pt());
                        pairInvMassPtdelR_bkg->Fill(etalv.M(), etalv.Pt(), deltaR);
                        delR_pairpT_bkg->Fill(deltaR, etalv.Pt());
		      }
		  } // inner cluster loop (buffered event)
	      } // inner cluster loop (current event)
	  } // event buffer loop

        // -- Manage the event buffer --
        // Add current event to the buffer
        MixedEvent new_event;
        new_event.zvtx = _vertex[2];
        new_event.centrality = _Centrality;
        new_event.clusters = current_clusters;
        event_buffer.push_back(new_event);

        // Remove the oldest event if buffer is full
        if (event_buffer.size() > (unsigned long int)bckgnd_evnts)
	  {
            event_buffer.pop_front();
	  }
	
      } // single event loop
  
  
  TTree * t1 = intree;
  if (!intree)
    {
      std::cout << "Error! Cannot get TTree. Exiting." << std::endl;
      exit(-1);
    }
  
  std::cout << "=======================================================" << std::endl;
  std::cout << "running in Event Mixing to construct background" << std::endl;

  // Set Branches
  //t1->SetBranchAddress("_runnumber",&_runnumber);
  //t1->SetBranchAddress("_nClusters", &_nClusters);
  //t1->SetBranchAddress("_Centrality", &_Centrality);
  //t1->SetBranchAddress("_vertex", _vertex);
  //t1->SetBranchAddress("_towerEnergy", &_towerEnergy);
  //t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  //t1->SetBranchAddress("_clusterPts", _clusterPts);
  //t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  //t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  //t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  //t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  //t1->SetBranchAddress("_clusterChi2", _clusterChi2);
  //t1->SetBranchAddress("_photonProb",&_photonProb);
  //t1->SetBranchAddress("ScaledTriggerBit",ScaledTriggerBit);

  for(size_t z = 0; z < eventBuffer.size();++z)
    {
      
      for(size_t c = 0; c < eventBuffer[z].size(); ++c)
	{
	  
	  std::vector<mixedEventStruct>& evts = eventBuffer[z][c];
	  
	  if(evts.empty() || evts.size() < 2)
	    {
	      continue;	    
	    }
       	  
	  //loop over all events in a given z/c bin
	  //evt1 is just the element number within evts, not the acutal numbered entry in the TTree
	  for(size_t evt1 = 0; evt1 < evts.size(); evt1++)
	    {

	      //copy over the evts, but this will be mixed randomly
	      std::vector<mixedEventStruct> mixed_evts = evts;
	      
	      std::shuffle(mixed_evts.begin(), mixed_evts.end(), std::mt19937{std::random_device{}()});

	      size_t nMix = std::min(mixed_evts.size() -1, (size_t)nMixedEvents);

	      // get clusters for current event
	      auto &current_clusters = evts[evt1].clusters;

	      //go to event 2
	      for(size_t evt2 = 0; evt2 < nMix; evt2++)
		{
		  if(evts[evt1].eventNumber == mixed_evts[evt2].eventNumber)
		    continue;
		  
		  auto &mix_clusters = mixed_evts[evt2].clusters;

		  if(mix_clusters.size() < 2)
		    continue;
		   
		  //loop over all clusters to make bkg
		  size_t ncc = current_clusters.size();
		  
		  for(size_t iCs = 0; iCs < ncc; iCs++) // loop over first CLUSTER 
		    {

		      if(current_clusters[iCs].Pt() < leadPhotonPtCut || current_clusters[iCs].Pt() < 0)
			continue;

		      size_t mcc = mix_clusters.size();
			
		      for(size_t jCs = 0; jCs < mcc; jCs++) 
			{
			  
			  if(mix_clusters[jCs].Pt() < subPhotonPtCut || mix_clusters[jCs].Pt() < 0)
			    continue;
			  
			  float mix_alpha = fabs(current_clusters[iCs].E() - mix_clusters[jCs].E()) / (current_clusters[iCs].E() + mix_clusters[jCs].E());

			  float deltaR = current_clusters[iCs].DeltaR(mix_clusters[jCs]);
		    
			  if (mix_alpha > asymmetryCut)
			    continue;

			  if (deltaR > deltaRCut)
			    continue;

			  h_diphotonBgAsymmetry->Fill(mix_alpha);
			  h_diphotonBgDeltaR->Fill(deltaR);

			  TLorentzVector diphotonLV = current_clusters[iCs] + mix_clusters[jCs];

			  float diphotonMass = diphotonLV.M();
			  float diphotonPt = diphotonLV.Perp();
			  float diphotonEta = diphotonLV.Eta();

			  if(diphotonMass < 0.01 || diphotonMass > 1.0)
			    continue;

			  if(diphotonLV.Pt() > 1.0  && diphotonEta < 1.0)
			    {
			      double weight = 1. / ( (double)ncc * (double)mcc * (double)nMixedEvents);

			      h_diphotonBgMass->Fill(diphotonMass);

			      h_diphotonBgMassVsPtWeighted->Fill(diphotonMass,diphotonPt,weight);

			      h_diphotonBgMassVsPt->Fill(diphotonMass,diphotonPt);

			      h_diphotonBgMassPtDeltaR->Fill(diphotonLV.M(),diphotonLV.Pt(),deltaR);

			      h_diphotonBgDeltaRVsPt->Fill(deltaR, diphotonLV.Pt());
			    }

			} //mix_cluster from evt2
	      
		    } //current_clusters from evt1
		  
		}//evt2 loop

	    }//evt1 loop
	  
	}//loop on cent bin
    
    }//loop on z bin

}//end event mixing f'n

*/


void Pi0EtaEfficiency::recoTruthMatch(eventTree &recoTree, truthEventTree &truthTree, double matchDeltaRCut, double asymmetryCut, double clusterChi2Cut, double diphotonDeltaRCut, double clusterEnergyCut, double leadPhotonPtCut, double subPhotonPtCut, double diphotonPtCut, TH1F *ptReweight)
{

/*	truth_cent_pt_mass3D->Write();
truth_cent_pt_mass3D_matched->Write();*/

  //make some histos
  h_recoMassVsRecoPt = new TH2F("h_recoMassVsRecoPt","Reco diphoton mass vs reco p_{T};reco meson p_{T} [GeV];reco M_{#gamma#gamma} [GeV]",30,0,30,50,0,1);

  h_truthMassVsTruthPt = new TH2F("h_truthMassVsTruthPt","Truth diphoton mass vs truth p_{T};truth meson p_{T} [GeV];truth M_{#gamma#gamma} [GeV]",30,0,30,50,0,1);

  h_recoMassVsRecoPtUnweighted = new TH2F("h_recoMassVsRecoPtUnweighted","Reco diphoton mass vs reco p_{T} (unweighted);reco meson p_{T} [GeV];reco M_{#gamma#gamma} [GeV]",30,0,30,50,0,1);

  h_truthMassVsTruthPtUnweighted = new TH2F("h_truthMassVsTruthPtUnweighted","Truth diphoton mass vs truth p_{T} (unweighted);truth meson p_{T} [GeV];truth M_{#gamma#gamma} [GeV]",30,0,30,50,0,1);

  h_truthDeltaR = new TH1F("h_truthDeltaR","Truth two-photon #DeltaR;#DeltaR_{truth};Counts",40,0,2);

  h_recoDeltaR = new TH1F("h_recoDeltaR","Reco two-cluster #DeltaR;#DeltaR_{reco};Counts",40,0,2);

  h_truthPhotonEnergy = new TH1F("h_truthPhotonEnergy","Truth photon energy;E_{#gamma}^{truth} [GeV];Counts",80,0,40);

  h_recoClusterEnergy = new TH1F("h_recoClusterEnergy","Reco cluster energy;E_{cluster}^{reco} [GeV];Counts",80,0,40);

 // h_truthPtVsCentrality = new TH2F("h_truthPtVsCentrality","Truth meson p_{T} vs centrality;truth meson p_{T} [GeV];centrality [%]",60,0,30,100,0,100);

  //h_recoPtVsCentrality = new TH2F("h_recoPtVsCentrality","Reco-matched meson p_{T} vs centrality;truth meson p_{T} [GeV];centrality [%]",60,0,30,100,0,100);

  //h_truthPtVsCentralityUnweighted = new TH2F("h_truthPtVsCentralityUnweighted","Truth meson p_{T} vs centrality (unweighted);truth meson p_{T} [GeV];centrality [%]",60,0,30,100,0,100);

//  h_recoPtVsCentralityUnweighted = new TH2F("h_recoPtVsCentralityUnweighted","Reco-matched meson p_{T} vs centrality (unweighted);truth meson p_{T} [GeV];centrality [%]",60,0,30,100,0,100);

  // ---- p+p efficiency vs pT (centrality-independent) ----
  // 1D denominator (all accepted truth mesons) and numerator (reco-matched +
  // reco-selected), same pT binning as the 2D histos. No centrality axis.
  h_truthPt = new TH1F("h_truthPt","Truth meson p_{T} (efficiency denominator);truth meson p_{T} [GeV];Counts",60,0,30);
  h_recoPt  = new TH1F("h_recoPt","Reco-matched meson p_{T} (efficiency numerator);truth meson p_{T} [GeV];Counts",60,0,30);
  h_truthPtUnweighted = new TH1F("h_truthPtUnweighted","Truth meson p_{T} (denominator, unweighted);truth meson p_{T} [GeV];Counts",60,0,30);
  h_recoPtUnweighted  = new TH1F("h_recoPtUnweighted","Reco-matched meson p_{T} (numerator, unweighted);truth meson p_{T} [GeV];Counts",60,0,30);
  h_truthPt->Sumw2();
  h_recoPt->Sumw2();

 // h_centralityFromImpactParam = new TH1F("h_centralityFromImpactParam","Centrality from impact parameter;centrality [%];Counts",100,0,100);

  h_truthPhotonAsymmetry = new TH1F("h_truthPhotonAsymmetry","Truth photon energy asymmetry;#alpha_{truth} = |E_{1}-E_{2}|/(E_{1}+E_{2});Counts",100,0,5);

  h_recoPhotonAsymmetry = new TH1F("h_recoPhotonAsymmetry","Reco cluster energy asymmetry;#alpha_{reco} = |E_{1}-E_{2}|/(E_{1}+E_{2});Counts",100,0,5);

  h_deltaRResidual = new TH1F("h_deltaRResidual","#DeltaR residual (reco - truth);#DeltaR_{reco} - #DeltaR_{truth};Counts",200,-1,1);
  h_energyResidual = new TH1F("h_energyResidual","Photon energy residual (reco - truth);E_{reco} - E_{truth} [GeV];Counts",200,-10,10);

  h_truthPtVsDeltaR = new TH2F("h_truthPtVsDeltaR","Truth meson p_{T} vs two-photon #DeltaR;truth meson p_{T} [GeV];#DeltaR_{truth}",120, 0, 30, 1500, 0, 1.5);

  h_truthPhotonEtaPhi = new TH2F("h_truthPhotonEtaPhi","Truth photon #eta-#phi;#eta_{#gamma}^{truth};#phi_{#gamma}^{truth}",96,-1.2,1.2, 256, -1*M_PI, M_PI);

  h_recoClusterEtaPhi = new TH2F("h_recoClusterEtaPhi","Reco cluster #eta-#phi;#eta_{cluster}^{reco};#phi_{cluster}^{reco}",96,-1.2,1.2, 256, -1*M_PI, M_PI);


 //Truth 3D histo filling
// --- 3D truth: (truth pT, centrality, truth mass)
//TH3D* truthPtCentralityMass3D = new TH3D("h_truthPtCentralityMass","Truth meson p_{T} vs centrality vs mass;truth meson p_{T} [GeV];centrality [%];truth M_{#gamma#gamma} [GeV]",
//                                     60, 0, 30,     // truth pT
//                                     100, 0, 100,   // centrality
//                                     50, 0, 1.0);   // truth mass


// --- 3D truth for reco-matched & selected candidates
//TH3D* truthPtCentralityMass3DMatched = new TH3D("h_truthPtCentralityMassMatched","Truth meson p_{T} vs centrality vs mass (reco-matched);truth meson p_{T} [GeV];centrality [%];truth M_{#gamma#gamma} [GeV]",
//                                             60, 0, 30,
//                                             100, 0, 100,
//                                             50, 0, 1.0);
 


  //-------------- end histo making -----------------//
  
  //loop is per-event
  for(int i = 0; i < truthTree.fChain->GetEntries(); i++)
    {
      
      truthTree.GetEntry(i);
      recoTree.GetEntry(i);

      //h_centralityFromImpactParam->Fill(getCentralityBin(truthTree._bimp));
      
      //keep track of which clusters are matched with truth photon
      std::set<int> usedClusters;
      
      //loop on etas within acceptance
      for(int k = 0; k < truthTree._nEtas; k++)
	{
    	  
	  TLorentzVector truthEta  = *(TLorentzVector *)truthTree.acceptedEtaVec->At(k);	  
	  TLorentzVector truthPho1 = *(TLorentzVector *)truthTree.acceptedPhoton1Vec->At(k);
	  TLorentzVector truthPho2 = *(TLorentzVector *)truthTree.acceptedPhoton2Vec->At(k);

	  // ==== get weights ===
	  // Weighted by default (MC pT reweight). setUseWeights(false) makes
	  // every fill unit-weight so the unweighted efficiency can be produced
	  // with the same binary (decision #6). The "weights" hist is still
	  // required to be non-null when weighting is on.
	  int weight_bin = ptReweight->FindBin(truthEta.Pt());
	  const double WEIGHT = m_useWeights ? ptReweight->GetBinContent(weight_bin) : 1.0;
	  //===

	  float truthPho1E = truthPho1.E();
	  float truthPho2E = truthPho2.E();
	  float truth_alpha = fabs((truthPho1E - truthPho2E)) / (truthPho1E + truthPho2E);
	  float truthEtaPt = truthEta.Pt();
	  float truthEtaM = truthEta.M();
	  float f_truthDelR = truthPho1.DeltaR(truthPho2);  

	  
	  h_truthPhotonAsymmetry->Fill(truth_alpha);

	  h_truthMassVsTruthPt->Fill(truthEtaPt,truthEtaM, WEIGHT);

	  h_truthMassVsTruthPtUnweighted->Fill(truthEtaPt, truthEtaM);

	  h_truthDeltaR->Fill(f_truthDelR);

	  h_truthPhotonEnergy->Fill(truthPho1E);

	  h_truthPhotonEnergy->Fill(truthPho2E);

	  h_truthPtVsDeltaR->Fill(truthEtaPt, f_truthDelR);

	  h_truthPhotonEtaPhi->Fill(truthPho1.Eta(), truthPho1.Phi());
	  h_truthPhotonEtaPhi->Fill(truthPho2.Eta(), truthPho2.Phi());

	  // "centrality" bin from the truth impact parameter. In p+p this is a
	  // FIXED inclusive bin (getCentralityBin returns 1), so the y-axis of the
	  // efficiency histos below is a single bin == the whole p+p sample. The
	  // efficiency is thus inclusive vs pT. (Histogram NAMES keep "Centrality"
	  // for backward compatibility with the extraction macros.)
	 // int converted_bimp = getCentralityBin(truthTree._bimp);

	 // h_truthPtVsCentrality->Fill(truthEta.Pt(), converted_bimp, WEIGHT);

	 // h_truthPtVsCentralityUnweighted->Fill(truthEta.Pt(), converted_bimp);

	  // p+p efficiency DENOMINATOR vs pT (no centrality)
	  h_truthPt->Fill(truthEta.Pt(), WEIGHT);
	  h_truthPtUnweighted->Fill(truthEta.Pt());

//	  truthPtCentralityMass3D->Fill(truthEtaPt, converted_bimp, truthEtaM, WEIGHT);

   /*
	  
	  if(truthPho1.Pt() < leadPhotonPtCut)
	    continue;

	  if(truthPho2.Pt() < subPhotonPtCut)
	    continue;

	  if (truth_alpha > asymmetryCut)
	    continue;

	 // if(f_truthDelR > diphotonDeltaRCut)
	 //   continue;

	  if(truthEtaM < 0.01 || truthEtaM > 1.0)
	    continue;

	  if (truthEtaPt < diphotonPtCut)
	    continue;

	  if(fabs(truthEta.Eta()) > 1)
	    continue;
	 */ 
	  
	  // ---- TRUTH-MATCHING (this is what makes the efficiency combinatorial-free) ----
	  // For THIS truth meson we already have its two decay photons (truthPho1,
	  // truthPho2). We find the reco cluster closest in (eta,phi) -> DeltaR to
	  // each truth photon, requiring DeltaR < matchDeltaRCut. "best" = smallest
	  // DeltaR. usedClusters[] prevents a cluster from matching two photons, and
	  // photon2 may not reuse photon1's cluster. If BOTH photons find a match the
	  // meson counts in the NUMERATOR -> it was successfully reconstructed.
	  int bestClus1 = -1;   // index of reco cluster matched to truth photon 1
	  int bestClus2 = -1;   // index of reco cluster matched to truth photon 2
	  double bestDR1 = 999; // smallest DeltaR found so far for photon 1
	  double bestDR2 = 999; // smallest DeltaR found so far for photon 2

	  //loop on reco clusters to get best match on truth photon1
	  for(int j = 0; j < recoTree._nClusters; j++)
	    {
	      //make TLV
	      TLorentzVector recoClus;
	      recoClus.SetPtEtaPhiE(recoTree._clusterPts[j], recoTree._clusterEtas[j], recoTree._clusterPhis[j], recoTree._clusterEnergies[j]);
	      
	      h_recoClusterEtaPhi->Fill(recoClus.Eta(), recoClus.Phi());

	      if (usedClusters.count(j))
		continue;

	      if(recoTree._clusterChi2[j] > clusterChi2Cut)
		continue;

	      if(recoTree._clusterPts[j] < leadPhotonPtCut)
		continue;

	      double dR1 = recoClus.DeltaR(truthPho1);
	      if (dR1 < bestDR1 && dR1 < matchDeltaRCut)
		{
		  bestDR1 = dR1;
		  bestClus1 = j;
		}
	     
	    }
	  
	  
	  //loop on reco clusters to get best match for truth photon2
	  for(int j = 0; j < recoTree._nClusters; j++)
	    {
	      
	      if (usedClusters.count(j))
		continue;
	      
	      if (j == bestClus1)
		continue;

	      if(recoTree._clusterChi2[j] > clusterChi2Cut)
		continue;

	      if(recoTree._clusterPts[j] < subPhotonPtCut)
		continue;

	      //make TLV
	      TLorentzVector recoClus;
	      recoClus.SetPtEtaPhiE(recoTree._clusterPts[j], recoTree._clusterEtas[j], recoTree._clusterPhis[j], recoTree._clusterEnergies[j]);

	      double dR2 = recoClus.DeltaR(truthPho2);
	      if (dR2 < bestDR2 && dR2 < matchDeltaRCut)
		{
		  bestDR2 = dR2;
		  bestClus2 = j;
		}
	    }

		  
	  //make diphoton
	  if (bestClus1 != -1 && bestClus2 != -1)
	    {
	      if(bestClus1 == bestClus2)
		continue;
	      
	      TLorentzVector recoPho1,recoPho2, recoEta;
	      
	      recoPho1.SetPtEtaPhiE(recoTree._clusterPts[bestClus1], recoTree._clusterEtas[bestClus1], recoTree._clusterPhis[bestClus1], recoTree._clusterEnergies[bestClus1]);
	      recoPho2.SetPtEtaPhiE(recoTree._clusterPts[bestClus2], recoTree._clusterEtas[bestClus2], recoTree._clusterPhis[bestClus2], recoTree._clusterEnergies[bestClus2]);
	      recoEta = recoPho1 + recoPho2;

	      double reco_alpha = fabs((recoPho1.E() - recoPho2.E())) / (recoPho1.E() + recoPho2.E());
	      double delR_reco = recoPho1.DeltaR(recoPho2);
	      double delR_res = delR_reco -  f_truthDelR;
	      double p1_res = recoPho1.E() - truthPho1E;
	      double p2_res = recoPho2.E() - truthPho2E;

	      //Fill some histos
	      h_energyResidual->Fill(p1_res);

	      h_energyResidual->Fill(p2_res);

	      h_deltaRResidual->Fill(delR_res);

	      h_recoDeltaR->Fill(delR_reco);

	      h_recoPhotonAsymmetry->Fill(reco_alpha,WEIGHT);

	      if(reco_alpha > asymmetryCut)
		continue;

	       if(delR_reco > diphotonDeltaRCut)
		continue;

	      // Species-driven reco diphoton mass acceptance (set via setMeson()).
	      if(recoEta.M() < m_massWinLo || recoEta.M() > m_massWinHi)
		continue;

	      if(recoEta.Pt() < diphotonPtCut)
		continue;

	      if(fabs(recoEta.Eta()) > 1)
		continue;


	      // Fill histograms
	      h_recoMassVsRecoPt->Fill(recoEta.Pt(), recoEta.M(),WEIGHT);

	      h_recoMassVsRecoPtUnweighted->Fill(recoEta.Pt(), recoEta.M());

       	     // h_recoPtVsCentrality->Fill(truthEta.Pt(),converted_bimp,WEIGHT);

	      //h_recoPtVsCentralityUnweighted->Fill(truthEta.Pt(),converted_bimp);

	      // p+p efficiency NUMERATOR vs pT (no centrality)
	      h_recoPt->Fill(truthEta.Pt(), WEIGHT);
	      h_recoPtUnweighted->Fill(truthEta.Pt());

	      h_recoClusterEnergy->Fill(recoPho1.E());

	      h_recoClusterEnergy->Fill(recoPho2.E());


               // Fill truth 3D but only for reco-matched + reco-selected candidates
              // truthPtCentralityMass3DMatched->Fill(truthEtaPt, converted_bimp, truthEtaM, WEIGHT);

	   	    		  	       
	      // Mark clusters as used
	      usedClusters.insert(bestClus1);
	      usedClusters.insert(bestClus2);
	    }
	  
	}//loop on number of etas
      
    }//loop on total events/entries


  //--------------- do some more stuff (mass resolution, efficiencies) ------------------------------




  double _ptCenters, _resolutionVals, _resolutionErrs;

  TTree *t = new TTree("t","mytree");
  t->Branch("_ptCenters", &_ptCenters, "_ptCenters/D");
  t->Branch("_resolutionVals", &_resolutionVals, "_resolutionVals/D");
  t->Branch("_resolutionErrs", &_resolutionErrs, "_resolutionErrs/D");


  std::vector<double> yields;
  std::vector<double> yieldsErr;

  double ptVals[9] = {1,2,3,3.3,4,5,6,7,8};
  double widths[8] = {0.5,0.5,0.15,0.35,0.5,0.5,0.5,0.5};


  for (int i = 0; i < 8; ++i)
  {
    int ptMin = ptVals[i];
    int ptMax = ptVals[i+1];
    double ptAvg = (ptMin+ptMax)/2.0;
    
    int binLow  = h_recoMassVsRecoPtUnweighted->GetXaxis()->FindBin(ptMin);
    int binHigh = h_recoMassVsRecoPtUnweighted->GetXaxis()->FindBin(ptMax);

    std::string projName = Form("mass_proj_pt_%d_%d", ptMin, ptMax);
    TH1D* h_proj = h_recoMassVsRecoPtUnweighted->ProjectionY(projName.c_str(), binLow, binHigh);

    // Species-driven fit range/seed (set via setMeson()): gaus peak on a
    // pol3 background. pi0 ~ [0.065,0.45] mean 0.135; eta ~ [0.35,0.75] mean 0.548.
    double lowFit = m_fitLo;
    double hiFit = m_fitHi;
    TF1* fitFcn = new TF1("fitFcn", " gaus(0) + pol3(3) ", lowFit, hiFit);


    fitFcn->SetParameters(1000, m_fitMeanSeed, 0.012);
    fitFcn->SetParLimits(0, 0, 1e9);
    fitFcn->SetParLimits(1, m_fitMeanLo, m_fitMeanHi);
    fitFcn->SetParLimits(2, 0.005, 0.075);

    fitFcn->SetLineColor(kRed);

    h_proj->Fit(fitFcn, "RQ");

    TF1* bestFn = (TF1*)fitFcn->Clone("bestFn");

    double mean = 0;
    double sigma = 0;
    double mean_err = 0;
    double sigma_err = 0;

    mean = bestFn->GetParameter(1);
    sigma = bestFn->GetParameter(2);
    
    mean_err = bestFn->GetParError(1);
    sigma_err = bestFn->GetParError(2);
      

    if (mean > 0.10 && mean < 0.17 && sigma > 0)
    {
      _ptCenters = ptAvg; 
      _resolutionVals = sigma / mean;

      double rel_err = sqrt(pow(sigma_err / sigma, 2) + pow(mean_err / mean, 2));
      _resolutionErrs = _resolutionVals * rel_err;

      t->Fill();

      yields.push_back( bestFn->Integral(mean - 3*sigma, mean + 3*sigma) / h_proj->GetBinWidth(1) );
      yieldsErr.push_back( sqrt(bestFn->Integral(mean - 3*sigma, mean + 3*sigma) / h_proj->GetBinWidth(1) ) );
    }
    
    h_proj->Write();
    
  }//end for loop on pt bins for mass peak fits

  // Plotting
  int nPoints = t->GetEntries();
  std::vector<double> x, y, ex, ey;
  
  for (int i = 0; i < nPoints; ++i)
  {
    t->GetEntry(i);
    x.push_back(_ptCenters);
    y.push_back(_resolutionVals);
    ex.push_back(widths[i]);
    ey.push_back(_resolutionErrs);

  }

  TGraphErrors *mass_res = new TGraphErrors(nPoints, &x[0], &y[0], &ex[0], &ey[0]);
  mass_res->SetTitle("Mass Resolution vs pT; #pi^{0}_{reco} p_{T} [GeV]; #sigma_{M} / M");
  mass_res->SetMarkerStyle(20);
  mass_res->SetMarkerColor(kBlue);
  mass_res->SetLineColor(kBlue);

  TGraphErrors *eta_yields = new TGraphErrors(nPoints, &x[0], &yields[0], &ex[0], &yieldsErr[0]);
  eta_yields->SetTitle("; p_{T} [GeV];E #frac{1}{N_{evt}} #frac{d^{3}N}{dp^{3}}");
  eta_yields->SetMarkerStyle(20);
  eta_yields->SetMarkerColor(kBlack);
  eta_yields->GetYaxis()->SetRangeUser(1e-11,1e3);
  

  t->Write();
  mass_res->Write("resolutionGraph");
  eta_yields->Write("etaYields");

  

  

  //---------------- efficiency vs pT (p+p, centrality-independent) ---------------------
  // Primary p+p efficiency: divide the 1D numerator by the 1D denominator
  // directly. No centrality projection, no dependence on the getCentralityBin
  // hack. "B" gives binomial errors appropriate for an efficiency (unweighted);
  // the weighted variant uses a plain divide (Sumw2 propagates the weights).

  h_truthPt->Write();
  h_recoPt->Write();
  h_truthPtUnweighted->Write();
  h_recoPtUnweighted->Write();

  h_efficiencyVsPt = (TH1F *)h_recoPt->Clone("h_efficiencyVsPt");
  h_efficiencyVsPt->SetDirectory(0);
  h_efficiencyVsPt->SetTitle("Efficiency vs p_{T} (p+p);truth meson p_{T} [GeV];Efficiency");
  h_efficiencyVsPt->SetYTitle("Efficiency");
  h_efficiencyVsPt->Divide(h_recoPt, h_truthPt);
  h_efficiencyVsPt->GetYaxis()->SetRangeUser(0,1);
  h_efficiencyVsPt->Write();

  h_efficiencyVsPtUnweighted = (TH1F *)h_recoPtUnweighted->Clone("h_efficiencyVsPtUnweighted");
  h_efficiencyVsPtUnweighted->SetDirectory(0);
  h_efficiencyVsPtUnweighted->SetTitle("Efficiency vs p_{T} (p+p, unweighted);truth meson p_{T} [GeV];Efficiency");
  h_efficiencyVsPtUnweighted->SetYTitle("Efficiency");
  h_efficiencyVsPtUnweighted->Divide(h_recoPtUnweighted, h_truthPtUnweighted, 1.0, 1.0, "B");
  h_efficiencyVsPtUnweighted->GetYaxis()->SetRangeUser(0,1);
  h_efficiencyVsPtUnweighted->Write();


  /*
  //---------------- efficiencies (legacy centrality-projected, kept for Au+Au) ---------------------

  TH2F *numerator = (TH2F *)h_recoPtVsCentrality->Clone("numerator");
  TH2F *denom = (TH2F *)h_truthPtVsCentrality->Clone("denom");


  TH2F *noW_numerator = (TH2F *)h_recoPtVsCentralityUnweighted->Clone("noW_numerator");
  TH2F *noW_denom = (TH2F *)h_truthPtVsCentralityUnweighted->Clone("noW_denom");

  numerator->SetDirectory(0);
  denom->SetDirectory(0);
  noW_numerator->SetDirectory(0);
  noW_denom->SetDirectory(0);

  int centBin[4] = {0,30,60,92};

  for(int i = 0; i < 3; i++)
    {
      if(i==1)continue; //will skip the 92/0 bin combo....
      
      int lbin = denom->GetYaxis()->FindBin(centBin[i]);
      int hbin = denom->GetYaxis()->FindBin(centBin[i+1]) -1;

      //with weight
      TH1D *numProj = numerator->ProjectionX(Form("Ncent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);
      TH1D *denProj = denom    ->ProjectionX(Form("Dcent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);

      numProj->SetDirectory(0);
      denProj->SetDirectory(0);

      //no weights
      TH1D *noW_numProj = noW_numerator->ProjectionX(Form("noW_Ncent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);
      TH1D *noW_denProj = noW_denom    ->ProjectionX(Form("noW_Dcent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);

      noW_numProj->SetDirectory(0);
      noW_denProj->SetDirectory(0);

      TH1D *eff = (TH1D *)numProj->Clone(Form("eff_%d_%d",centBin[i],centBin[i+1]));
      eff->SetXTitle("#eta_{truth} p_{T} [GeV]");
      eff->SetYTitle("Efficiency");
      eff->SetTitle(Form("%d-%d%% Centrality",centBin[i],centBin[i+1]));
      eff->Divide(denProj);
      eff->GetYaxis()->SetRangeUser(0,1);
      eff->Write();
      eff=nullptr;

      TH1D *noW_eff = (TH1D *)noW_numProj->Clone(Form("noW_eff_%d_%d",centBin[i],centBin[i+1]));
      noW_eff->SetXTitle("#eta_{truth} p_{T} [GeV]");
      noW_eff->SetYTitle("Efficiency");
      noW_eff->SetTitle(Form("%d-%d%% Centrality",centBin[i],centBin[i+1]));
      noW_eff->Divide(noW_denProj);
      noW_eff->GetYaxis()->SetRangeUser(0,1);
      noW_eff->Write();
      noW_eff=nullptr;
    }

*/

/*
 
  //---------------- efficiencies ---------------------

  TH2F *numerator = (TH2F *)reco_cent_pt->Clone("numerator");
  TH2F *denom = (TH2F *)truth_cent_pt->Clone("denom");

 
  TH2F *noW_numerator = (TH2F *)noW_reco_cent_pt->Clone("noW_numerator");
  TH2F *noW_denom = (TH2F *)noW_truth_cent_pt->Clone("noW_denom");

  numerator->SetDirectory(0);
  denom->SetDirectory(0);
  noW_numerator->SetDirectory(0);
  noW_denom->SetDirectory(0);

  int centBin[7] = {0,20,40,60,92,0,92};

  for(int i = 0; i < 6; i++)
    {
      if(i==4)continue; //will skip the 92/0 bin combo....
      
      int lbin = denom->GetYaxis()->FindBin(centBin[i]);
      int hbin = denom->GetYaxis()->FindBin(centBin[i+1]);

      //with weight
      TH1D *numProj = numerator->ProjectionX(Form("Ncent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);
      TH1D *denProj = denom    ->ProjectionX(Form("Dcent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);

      numProj->SetDirectory(0);
      denProj->SetDirectory(0);

      //no weights
      TH1D *noW_numProj = noW_numerator->ProjectionX(Form("noW_Ncent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);
      TH1D *noW_denProj = noW_denom    ->ProjectionX(Form("noW_Dcent_%d_%d",centBin[i],centBin[i+1]),lbin,hbin);

      noW_numProj->SetDirectory(0);
      noW_denProj->SetDirectory(0);

      TH1D *eff = (TH1D *)numProj->Clone(Form("eff_%d_%d",centBin[i],centBin[i+1]));
      eff->SetXTitle("#pi^{0}_{truth} p_{T} [GeV]");
      eff->SetYTitle("Efficiency");
      eff->SetTitle(Form("%d-%d%% Centrality",centBin[i],centBin[i+1]));
      eff->Divide(denProj);
      eff->GetYaxis()->SetRangeUser(0,1);
      eff->Write();
      eff=nullptr;

      TH1D *noW_eff = (TH1D *)noW_numProj->Clone(Form("noW_eff_%d_%d",centBin[i],centBin[i+1]));
      noW_eff->SetXTitle("#pi^{0}_{truth} p_{T} [GeV]");
      noW_eff->SetYTitle("Efficiency");
      noW_eff->SetTitle(Form("%d-%d%% Centrality",centBin[i],centBin[i+1]));
      noW_eff->Divide(noW_denProj);
      noW_eff->GetYaxis()->SetRangeUser(0,1);
      noW_eff->Write();
      noW_eff=nullptr;
    }

*/








  
}//end recoTruthMatch




/*

int Pi0EtaEfficiency::getCentralityBin(double b_imp)
{
  // PHYSICS: this converts an Au+Au Glauber IMPACT PARAMETER b (in fm) into a
  // 1..100 centrality percentile via fixed b-edges (small b = central). It is
  // meaningful ONLY for nucleus-nucleus collisions.
  //
  // p+p (PYTHIA): there is no impact-parameter geometry, so we return a SINGLE
  // inclusive bin (1) for every event. The (pT vs centrality) efficiency
  // histograms therefore become 1 bin wide in centrality == an inclusive p+p
  // efficiency vs pT. (Returning 1 also matches the old numeric behavior, where
  // p+p's bimp ~ 0 already mapped to the first bin -- so results are unchanged.)
  if (m_isPP) return 1;

  // Edges are fixed -> build the lookup table ONCE (static), not on every call
  // (this function runs per-meson, per-event).
  static const std::vector<double> b_edges = {1.58,2.23,2.73,3.15,3.52,3.85,4.15,4.43,4.69,4.94,5.18,5.4,5.62,5.83,6.03,6.22,6.4,6.59,6.76,6.94,7.11,7.27,7.43,7.58,7.73,7.88,8.03,8.17,8.31,8.45,8.59,8.73,8.86,8.99,9.12,9.25,9.37,9.49,9.61,9.73,9.85,9.97,10.08,10.2,10.31,10.42,10.53,10.64,10.75,10.85,10.96,11.07,11.17,11.27,11.37,11.47,11.57,11.67,11.77,11.87,11.97,12.06,12.16,12.25,12.34,12.44,12.53,12.62,12.71,12.8,12.89,12.98,13.07,13.16,13.25,13.34,13.43,13.51,13.6,13.69,13.78,13.87,13.96,14.05,14.15,14.24,14.34,14.44,14.55,14.66,14.78,14.9,15.04,15.19,15.36,15.56,15.8,16.13,16.64,19.99};

  auto it = std::lower_bound(b_edges.begin(), b_edges.end(), b_imp);
  
  if (it != b_edges.end())
    {
      return std::distance(b_edges.begin(), it) + 1;  // Centrality bin 1–100
    }
  else
    {
      std::cout << "WARNING! Impact parameter outside of range! " << PHWHERE << std::endl;
      return -1;  // b_val outside range
    }
 
}

*/
void Pi0EtaEfficiency::setBufferSize(float z, int clo, int chi)
{
  std::cout << "Making event buffer ... " << std::endl;
  
  int nZ = (2 * static_cast<int>(z)) + 1;
  int nCent = (chi - clo) + 1;

  std::cout << "Size of buffer (z,cent): " << nZ << ", " << nCent << std::endl;
  
  eventBuffer.resize(nZ);
  for (int i = 0; i < nZ; ++i)
    eventBuffer[i].resize(nCent);
}








