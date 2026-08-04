#include "CaloCalibEmc_eta.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>




// Tower includes
#include <calobase/RawTowerGeom.h>
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



//Trigger
#include <calotrigger/TriggerRunInfov1.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <calotrigger/MinimumBiasInfo.h>



//GL1
#include <ffarawobjects/Gl1Packet.h>


// phool
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>


// MBD
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>


#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>


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

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TVectorF.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TList.h>

#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

#include <algorithm>  // for max, max_element
#include <cmath>      // for abs
#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>   // for vector
#include <TArray.h>

#include <centrality/CentralityInfo.h>    // for CentralityInfo, CentralityI...
#include <centrality/CentralityInfov1.h>
#include <centrality/CentralityInfov2.h>


#include <phool/getClass.h>
#include <phool/phool.h>

#include "TRandom3.h"
#include <string>   // for string
#include <random>
#include <fstream>

#include <algorithm>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <TRandom3.h>
#include <ctime>
#include <random>
#include <TMath.h>

#include "TVector3.h"
#include <TLorentzVector.h>
#include <TVector3.h>

// G4
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

//using namespace std;

std::vector<int> bad_towers;

//____________________________________________________________________________..
CaloCalibEmc_eta::CaloCalibEmc_eta(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  ,m_ievent(0)
  , m_Filename(filename)
  , cal_output(0)
  , _caloname("CEMC")
  //, fitp1_eta_phi2d(0)
  //, pairInvMassTotal(0)
  ,_eventTree(0)
  ,_eventNumber(-1)
  ,_nClusters(-1)
  ,maxTowerEta(-1)
  ,maxTowerPhi(-1)
  ,alphaCut(-1.0)   
  ,f_temp(0)  
{
	// different type of fore-ground we are making
  pairInvMassTotal = NULL;
  pairInvMassTotalBkgd = NULL;

	// different types of background we are making
	pairInvMassPtEta = NULL;
	pairInvMassPtEtaBkgd = NULL;
	
	//m_cent_nclus_cut = 1000;

}

//____________________________________________________________________________..
int CaloCalibEmc_eta::InitRun(PHCompositeNode *topNode)
{
  std::cout << "CaloCalibEmc_eta::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_ievent = 0;

  cal_output = new TFile(m_Filename.c_str(), "RECREATE");

  pairInvMassTotal = new TH1F("pairInvMassTotal", "Pair_Mass_Histo", 300, 0.0, 3.0);
  pairInvMassTotal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  pairInvMassTotal->SetYTitle("Counts");

	pairInvMassTotalBkgd = new TH1F("pairInvMassTotalBkgd", "Pair_Mass_Histo_Bkgd", 300, 0.0, 3.0);

  pairInvMassPtEta = new TH3F("pairInvMassPtEta", "Pair_Mass_Histo_PtEta", 300, 0.0, 3.0, 200, 0.0, 20.0, 22, -1.1, 1.1);
  pairInvMassPtEtaBkgd = new TH3F("pairInvMassPtEtaBkgd", "Pair_Mass_Histo_PtEta_Bkgd", 300, 0.0, 3.0, 200, 0.0, 20.0, 22, -1.1, 1.1);

  pairInvMassPtdelR = new TH3F("pairInvMassPtdelR", "Pair_Mass_Histo_PtdelR", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);
  pairInvMassPtdelRBkgd = new TH3F("pairInvMassPtdelRBkgd", "Pair_Mass_Histo_PtdelR_Bkgd", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);

	pairpTDelPhiDelEta = new TH3F("pairpTDelPhiDelEta", "pairpT_delphi_deleta", 200, 0.0, 20.0, 100, -3.2, 3.2, 110, -1.1, 1.1);
	pairpTDelPhiDelEtaBkgd = new TH3F("pairpTDelPhiDelEtaBkgd", "pairpT_delphi_deleta", 200, 0.0, 20.0, 100, -3.2, 3.2, 110, -1.1, 1.1);
		
	DelR_pairpT_f = new TH2F("DelR_pairpT_f", "DelR_pairpT", 300, 0.0, 1.5, 200, 0.0, 20.0);
	DelR_pairpT_b = new TH2F("DelR_pairpT_b", "DelR_pairpT", 300, 0.0, 1.5, 200, 0.0, 20.0);
  
	h_etaphi_clus = new TH2F("h_etaphi_clus", "etaphi_clus", 140, -1.2, 1.2, 64, -1 * M_PI, M_PI);  

	//for different run # range that have different trigger names for certain indices
        //h2_scaledTriggerBit = new TH1I("h2_scaledTriggerBit","",63,0,63);
       // h2_scaledTriggerBit->GetXaxis()->SetNdivisions(64);
       // h2_scaledTriggerBit->GetXaxis()->SetLabelSize(0.024);


     //  bkg_eff = new TH1F("bkg_eff","",1,0,1);
     //  bkg_eff->SetXTitle("");
    //   bkg_eff->SetYTitle("Efficiency");

    //   centrality_count = new TH1F("centrality_count","",1,0,1);


 // deltaR_fg = new TH1F("deltaR_fg","#Delta R (p_{T}^{1}, p_{T}^{2}) Foreground",2000.0, -10.0, 10.0);
  //deltaR_fg->SetXTitle("#Delta R");
 // deltaR_fg->SetYTitle("Counts");

 // deltaR_bg = new TH1F("deltaR_bg","#Delta R (p_{T}^{1}, p_{T}^{2}) Background",2000.0, -10.0, 10.0);
 // deltaR_bg->SetXTitle("#Delta R");
 // deltaR_bg->SetYTitle("Counts");

  //asym_fg = new TH1F("asym_fg","Asymmetry (Foreground)", 110.0, 0.0, 1.1);
 // asym_fg->SetXTitle("Assymetry Value");
 // asym_fg->SetYTitle("Counts");

 // asym_bg = new TH1F("asym_bg","Asymmetry (Background)", 110.0, 0.0, 1.1);
  //asym_bg->SetXTitle("Assymetry Value");
  //asym_bg->SetYTitle("Counts");

  //ieta_iphi_notGood = new TH2F("ieta_iphi_notGood","",96.0, -0.5, 95.5, 256.0, -0.5, 255.5);
 // ieta_iphi_notGood->SetXTitle("#eta bin");
 // ieta_iphi_notGood->SetYTitle("#phi bin");
 // ieta_iphi_notGood->SetZTitle("Flagged as 'notGood'");

 // ieta_iphi_energy = new TH3F("ieta_iphi_energy","",96, -0.5, 95.5, 256, -0.5, 255.5, 20.0, 0.0, 10.0);
 // ieta_iphi_energy->SetXTitle("#eta bin");
 // ieta_iphi_energy->SetYTitle("#phi bin");
//  ieta_iphi_energy->SetZTitle("Energy [GeV]");


  // eta_bg_invmass_pt = new TH2F("eta_bg_invmass_pt", "",3000, 0.0, 3.0, 290, 1.0, 30.0);
 // eta_bg_invmass_pt->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
 // eta_bg_invmass_pt->SetYTitle("pair p_{T} [GeV]");

  //clusterMap = new TH2F("clusterMap","Cluster Position", 96, -0.5, 95.5, 256, -0.5, 255.5);
  //clusterMap->SetXTitle("Cluster #eta");
  //clusterMap->SetYTitle("Cluster #phi");

  //clusterEnergy = new TH1F("clusterEnergy","",160,0,40);
 // clusterEnergy->SetXTitle("Cluster Energy [GeV]");
 // clusterEnergy->SetYTitle("Counts / (0.250 MeV)");

  //diPhotonMap = new TH2F("diPhotonMap","diPhoton Map",  96, -1.2, 1.2, 256, -3.2, 3.2);
  //diPhotonMap->SetXTitle("Reco #eta");
  //diPhotonMap->SetYTitle("Reco #phi");

  //leadPhoton = new TH3F("leadPhoton","Leading Photon Map",96, -0.5, 95.5, 256, -0.5, 255.5, 50, 0.0, 25.0);
  //leadPhoton->SetXTitle("#eta Bin");
  //leadPhoton->SetYTitle("#phi Bin");
  //leadPhoton->SetZTitle("p_{T}");

 // z_vertex = new TH1F("z_vertex","z_vertex",3000,-150,150);
//  z_vertex->SetXTitle("z_vertex [cm]");
 // z_vertex->SetYTitle("Counts");

  //evtCentrality = new TH1F("evtCentrality","",100,0,100);
 // evtCentrality->SetXTitle("Centrality Class");
  //evtCentrality->SetYTitle("Counts")

	if (topNode != 0)
  {
    // TTree declare
    _eventTree = new TTree("_eventTree", "An event level info Tree");
    // TTree branches
   // _eventTree->Branch("_runnumber",&_runnumber);
    _eventTree->Branch("_eventNumber", &_eventNumber, "_eventNumber/I");
    _eventTree->Branch("_nClusters", &_nClusters, "_nClusters/I");
    _eventTree->Branch("_nCentrality", &_nCentrality, "_nCentrality/F");
    _eventTree->Branch("_vertex", _vertex, "_vertex[3]/F"); // save co-ordinates of global-vertex (x, y and z co-ordinates)
    //_eventTree->Branch("_nTowers", &_nTowers, "_nTowers[_nClusters]/I");
    //_eventTree->Branch("_clusterIDs", _clusterIDs, "_clusterIDs[_nClusters]/I");
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
    // histlumicount = new TH1F("histlumicount",";;",1,0,1);
     //for mc single particle
 // _eventTree->Branch("truth_vertex",truth_vertex,"truth_vertex[3]/F");


  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_eta::process_event(PHCompositeNode *topNode)
{
	// list of all cuts we are going to use within process event
	float pt1_cut = 0.0;
	float pt2_cut = 0.0;
	int nCluster_cut = 10;
	float alpha_cut = 0.8;
	float delR_cut = 1.1;

  if (m_ievent % 500 == 0)
  {
    std::cout << std::endl;
    std::cout << "Beginning of the event " << m_ievent << std::endl;
    std::cout << "====================================" << std::endl;
  }

  _eventNumber = m_ievent++;

  // create a cluster object
  std::string clusnodename = "CLUSTERINFO_" + _caloname; // two options (pos-corr node and raw cluster
  //std::string clusnodename = "CLUSTER_" + _caloname; // Cluster Cemc

  RawClusterContainer *recal_clusters = findNode::getClass<RawClusterContainer>(topNode,clusnodename.c_str());

  if (!recal_clusters)
		{
			std::cout << PHWHERE << clusnodename << " node is missing" << std::endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	
	// get towers and tower-level information
	TowerInfoContainer *_towerinfos = NULL;
		
	std::string towernode = "TOWERINFO_CALIB_" + _caloname;
	//std::string towernode = "TOWERS_" + _caloname;
	_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towernode.c_str());
			
	if (!_towerinfos) 
		{
			std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}

  // create a tower geometry object
	// geometry object is same in tower info and raw tower
  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());

  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnode << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get Vertex
  float vx = 0;
  float vy = 0;
  float vz = 0;
		
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
	
/*	
  if ((!vertexmap) || (vertexmap->empty()))
  {
    std::cout << "Vertexmap either does not exist or is empty" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
	}
	

  if (vertexmap)
  {
    if (!vertexmap->empty())
    {
      GlobalVertex *vtx = (vertexmap->begin()->second);
      vx = vtx->get_x();
      vy = vtx->get_y();
      vz = vtx->get_z();
    }
  }
	
	if (abs(vz) < 30.0) {
	// saving vertex in ntuple
	_vertex[0] = vx;
	_vertex[1] = vy;
	_vertex[2] = vz;
	}

*/ 

float noVtx = 0;

if (vertexmap)
    {
      if (!vertexmap->empty())
	{
	  GlobalVertex *vtx = vertexmap->begin()->second;
	  if(!vertexmap->isValid() || !vtx)
	    {
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

  if(abs(vz) > 150)
    {      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _vertex[0] = vx;
  _vertex[1] = vy;
  _vertex[2] = vz;


/*
  // Get Vertex
  float truth_vx = 0;
  float truth_vy = 0;
  float truth_vz = 0;




  //======= truth container =======
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");  
  if(!truthinfo)
    {
      std::cout << PHWHERE << "Could not find truthinfo node G4TruthInfo" << std::endl;
    }

  
  // ======= truth vertex stuff ========
  if (truthinfo)
    {
      PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      truth_vx = gvertex->get_x();
      truth_vy = gvertex->get_y();
      truth_vz = gvertex->get_z();
    }

  
  if(abs(truth_vz) > 30.0)
    return Fun4AllReturnCodes::ABORTEVENT;

  truth_vertex[0] = truth_vx; 
  truth_vertex[1] = truth_vy; 
  truth_vertex[2] = truth_vz; 



  //z_vertex->Fill(vz);
 
*/
 

 

 // ====== GL1 Information ======
  Gl1Packet *gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
   
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

    

  
  // ======= MinBias stuff =======
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
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  else if(isMinBias)
    {
      evt_isMB++;
    }   

 	
  //get centrality node
  CentralityInfov2 *cent = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");
  if (!cent)
    {
      std::cout << " Centrality Node is missing." << std::endl;
       return Fun4AllReturnCodes::ABORTEVENT;

    }

  //get centrality of event
  float cent_val = cent->get_centrality_bin(CentralityInfo::PROP::mbd_NS);
  _nCentrality = cent_val;
  //evtCentrality->Fill(centval); 
	
  
	CLHEP::Hep3Vector vertex(vx, vy, vz); // creating 3-vector

  // ------------------------------

  // loop over the clusters
  RawClusterContainer::ConstRange t_rbegin_end = recal_clusters->getClusters();
  RawClusterContainer::ConstIterator t_rclusiter;

  RawCluster *savCs[10000];  // savingClusters that has 1 GeV or more
  int iCs = 0;

  // saving the clusters
  for (t_rclusiter = t_rbegin_end.first; t_rclusiter != t_rbegin_end.second; ++t_rclusiter)
  {
    RawCluster *t_recalcluster = t_rclusiter->second;

    //float cluse = t_recalcluster->get_ecore();
    //if ((t_recalcluster->get_ecore() > 0.6) || (t_recalcluster->get_chi2() < 4.0))  savCs[iCs++] = t_recalcluster; // only taking cluster with energy (e-core) more than 0.6 GeV
    if ((abs(vz) < 20.0) && (t_recalcluster->get_ecore() > 0.6) && (t_recalcluster->get_chi2() < 4.0))  savCs[iCs++] = t_recalcluster; // only taking cluster with energy (e-core) more than 0.6 GeV
  }

  _nClusters = iCs;

	//if (_nClusters > m_cent_nclus_cut)
	//	return Fun4AllReturnCodes::EVENT_OK;

  // looping on the saved clusters savCs
  // outer loop (we want to do pair of the loops to get invariant mass distribution)
  for (int jCs = 0; jCs < _nClusters; jCs++)
  {
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*savCs[jCs], vertex);



    //_clusterIDs[jCs] = savCs[jCs]->get_id();

		_clusterChi2[jCs] = savCs[jCs]->get_chi2();
		 _photonProb[jCs] = savCs[jCs]->get_prob();



    //vector to hold all the towers etas, phis, and energy in this cluster
    std::vector<int> toweretas;
    std::vector<int> towerphis;
    std::vector<float> towerenergies;
	
    TowerInfo * towinfo;
		
    // loop over the towers from the outer loop cluster
    // and find the max tower location and save the 
    // histogram on that max tower location with this 
    // invariant mass
    
    RawCluster::TowerConstRange towers = savCs[jCs]->get_towers();
    RawCluster::TowerConstIterator toweriter;
    
    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
      {
				//RawTowerGeom *tower_geom_instance = geom->get_tower_geometry(toweriter->first);
				//unsigned int towerindex = _towerinfos->decode_key(toweriter->first);
				int iphi = RawTowerDefs::decode_index2(toweriter->first);  // index2 is phi in CYL
				int ieta = RawTowerDefs::decode_index1(toweriter->first);  // index1 is eta in CYL	    
				unsigned int towerkey = iphi + (ieta << 16U);
				unsigned int towerindex =  _towerinfos->decode_key(towerkey);
	    
				towinfo = _towerinfos->get_tower_at_channel(towerindex);
	    
				double towerenergy = towinfo->get_energy();

			  // put the eta, phi, energy into corresponding vectors
				toweretas.push_back(ieta);
				towerphis.push_back(iphi);
				towerenergies.push_back(towerenergy);
       }
    
    // cout << endl;
    // cout << "Cluster energy: " << tt_clus_energy << endl;
    // cout << "Total number of towers (getNTowers()): " << savCs[jCs]->getNTowers() << endl;
    // cout << "Total number of towers size(toweretas): " << toweretas.size() << endl;
    // float maxTowerEnergy = *max_element(towerenergies.begin(), towerenergies.end());
    // cout << "The maxTowerEnergy: " << maxTowerEnergy << endl;

    int maxTowerIndex = max_element(towerenergies.begin(), towerenergies.end()) - towerenergies.begin();
    maxTowerEta = toweretas[maxTowerIndex];
    maxTowerPhi = towerphis[maxTowerIndex];
		
    _maxTowerEtas[jCs] = maxTowerEta;
    _maxTowerPhis[jCs] = maxTowerPhi;
   
	 //std::cout << E_vec_cluster.mag() << std::endl; 
		
		// extracting few values and setting to Ttree
		// most important area to convert ADC to GeV ("coversion_factor" is used for conversion)
    float tt_clus_energy = E_vec_cluster.mag();
    float tt_clus_eta = E_vec_cluster.pseudoRapidity();
    float tt_clus_pt = E_vec_cluster.perp();
    float tt_clus_phi = E_vec_cluster.getPhi();

    _clusterEnergies[jCs] = tt_clus_energy;
    _clusterPts[jCs] = tt_clus_pt;
    _clusterEtas[jCs] = tt_clus_eta;
    _clusterPhis[jCs] = tt_clus_phi;

    
		if ((tt_clus_pt < pt1_cut) || (_nClusters < nCluster_cut)) {continue;} //first photon cut and centrality cut (based on number of clusters)
				
		// another loop to go into the saved cluster
		for (int kCs = 0; kCs < _nClusters; kCs++)
			{
				if (jCs == kCs) {continue;}

				CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*savCs[kCs], vertex);
	      
				float tt2_clus_energy = E_vec_cluster2.mag();
				float tt2_clus_pt = E_vec_cluster2.perp();
							
				if (tt2_clus_pt < pt2_cut) {continue;}  // second photon cut
				
				// applying alpha cut (energy assymetry cut)
				alphaCut = abs(tt_clus_energy - tt2_clus_energy) / (tt_clus_energy + tt_clus_energy);
				
				if (alphaCut > alpha_cut) {continue;}
						
				float tt2_clus_eta = E_vec_cluster2.pseudoRapidity();
				float tt2_clus_phi = E_vec_cluster2.getPhi();
		      
				TLorentzVector pho1, pho2, mesonlv;
		      
				pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
				pho2.SetPtEtaPhiE(tt2_clus_pt, tt2_clus_eta, tt2_clus_phi, tt2_clus_energy);
		      
				if (pho1.DeltaR(pho2) > delR_cut) {continue;}
											
				mesonlv = pho1 + pho2;
				float pairInvMass = mesonlv.M();

				//std::cout << pairInvMass << std::endl;
		      

				if (fabs(mesonlv.M()) > 0.02)
					{ 
						pairInvMassTotal->Fill(pairInvMass);
					 }

			 } // closing for second loop cluster
		} // closing for first loop cluster
							
							
  _eventTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_eta::End(PHCompositeNode * topNode)
{


  //std::cout << "Total events: " << m_ievent << std::endl;
 // std::cout << "Fraction of events that are MB: " << (double)evt_isMB/m_ievent << std::endl;
 // std::cout << "Fraction of events that are not MB: " << (double)evt_isNotMB/m_ievent << std::endl;
 // std::cout << "Fraction of events w/ no vtx: " << (double)noVtx/m_ievent<< std::endl;




  if (topNode == 0 && f_temp)
    {
      cal_output->Close();
      f_temp->Close();
      delete f_temp;
      delete cal_output;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  
  cal_output->cd();
	//_eventTree->Write();
  cal_output->Write();
  cal_output->Close();
   // histlumicount->SetBinContent(1,evt_isMB);
   // histlumicount->Write();
  delete cal_output;

  return Fun4AllReturnCodes::EVENT_OK;
}


//______________________________________________________________________________..
// get the inv mass histogram for pi0/eta meson after collecting information from Ttree
// for time being we can consider this as a foreground
void CaloCalibEmc_eta::Loop(int nevts, TString _filename, TTree * intree)
{
  std::cout << "running in loop mode to get invariant mass plot" << std::endl;
	
	// list of parameters used in the loop
	float pt1_cut = 0.6;
	float pt2_cut = 0.5;
	float alpha_cut = 0.6;
	int max_nCluster_cut = 650;
	int min_nCluster_cut = 100;
	float delR_cut = 1.1;
	
	float chi2_cut = 4.0; // chi2 cut (this removes clusters that do not represents electron/photon)
	float zver_cut = 20.0;

  TTree * t1 = intree;
  if (!intree)
  {
		TFile *f = new TFile(_filename);
    t1 = (TTree *) f->Get("_eventTree");
  }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
	//t1->SetBranchAddress("_nCentrality", &_nCentrality);
  //t1->SetBranchAddress("_clusterIDs", _clusterIDs);
	t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
	
	t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
	t1->SetBranchAddress("_clusterChi2", _clusterChi2);


  // pre-loop to save all the clusters LorentzVector
  //TLorentzVector *savClusLV[10000];

	// making vector that has TLorentz vectors
	std::vector<TLorentzVector> savClusLV1;
	std::vector<TLorentzVector> savClusLV;
	std::vector<int> cluster_coll;

  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if ((nevts < 0) || (nEntries < nevts))
    nevts2 = nEntries;
	
  for (int i = 0; i < nevts2; i++) // for each event
  {
		// clear the vector after every events
		cluster_coll.clear();

    // load the ith instance of the TTree
    t1->GetEntry(i);

    if (i % 1000 == 0) {std::cout << "evt no " << i << std::endl;}

		//std::cout << "Number of Clusters = " << nClusters << std::endl;
    
		if ((_nClusters > max_nCluster_cut) || (_nClusters <  min_nCluster_cut)) continue; // ignoring higher cluster events (more central events)
		
		float vertex_z = _vertex[2]; // co-ordinates are 0, 1, 2
		if (abs(vertex_z) > zver_cut) continue;
		

		savClusLV1.clear();
		for (int j = 0; j < _nClusters; j++) // for each cluster in event
    {
			if (_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
			if (_clusterChi2[j] > chi2_cut) continue; 

		  TLorentzVector LV_Clus;

		  LV_Clus.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
		  savClusLV1.push_back(LV_Clus);
		 }

		//////////////////////////////////////////////////////////////////////////////////////////////////
    for (size_t a = 0; a < savClusLV1.size(); a++) // outer loop to iterate clusters
    {
			TLorentzVector& pho_a = savClusLV1[a];

			//if (fabs(pho_a.Pt()) < pt1_cut) continue; //second P_T cut
      
			for (size_t b = a+1; b < savClusLV1.size(); b++) // inner loop to iterate clusters, remove double counting
      {

        TLorentzVector& pho_b = savClusLV1[b];
				TLorentzVector mes_ab;

				if (pho_a.DeltaR(pho_b) > 0.2) continue;

				//alphaCut = fabs((pho_a.E() - pho_b.E())/(pho_a.E()+ pho_b.E()));

				//if (alphaCut < 0.25) continue;

				
				/*
				if (fabs(pho_b.Pt()) < pt2_cut) continue; //second P_T cut
							
				alphaCut = fabs((pho_a.E() - pho_b.E())/(pho_a.E()+ pho_b.E()));

				if (alphaCut > alpha_cut) continue;
				
				//if (pho1->Eta()/pho2->Eta() < 0) continue;
				
				if (pho_a.DeltaR(pho_b) > delR_cut) continue;
				if (pho_a.DeltaR(pho_b) > 0.35) continue;
				//if (pho_a.DeltaR(pho_b) < 0.05) continue;
				*/
				
        //mes_ab = pho_a + pho_b;

					//if ((mes_lv.M() < 0.02) || (mes_lv.Pt() < 1.0)) continue;
				//if((mes_ab.M() > 0.05) && (mes_ab.M() < 0.25))
					//std::cout << _delR << std::endl;
				cluster_coll.push_back(a);
				cluster_coll.push_back(b);
			 
			 }
		 }


		/////////
		// now get update from cluster_coll to savClusLV
	  // Sort the vector
    std::sort(cluster_coll.begin(), cluster_coll.end());

    // Remove duplicates
    auto last = std::unique(cluster_coll.begin(), cluster_coll.end());

    // Erase the 'extra' elements after unique to resize the vector
    cluster_coll.erase(last, cluster_coll.end());

		savClusLV.clear();
		// update to to the savClusLV
		for (size_t k = 0; k < cluster_coll.size(); k++) // for each modified one
    {
			savClusLV.push_back(savClusLV1[cluster_coll[k]]);
		 }

		/////////////////////////////////////////////////////////

    for (size_t iCs = 0; iCs < savClusLV.size(); iCs++) // outer loop to iterate clusters
    {
			TLorentzVector& pho1 = savClusLV[iCs];

      
			if (fabs(pho1.Pt()) < pt1_cut)	continue; // first P_T cut

			h_etaphi_clus->Fill(pho1.Eta(), pho1.Phi());

      for (size_t jCs = iCs+1; jCs < savClusLV.size(); jCs++) // inner loop to iterate clusters, remove double counting
      {

        TLorentzVector& pho2 = savClusLV[jCs];
				
				if (fabs(pho2.Pt()) < pt2_cut) continue; //second P_T cut
							
				alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));

				if (alphaCut > alpha_cut) continue;
				
				//if (pho1->Eta()/pho2->Eta() < 0) continue;
				
				if (pho1.DeltaR(pho2) > delR_cut) continue;
				
				// getting deleta and delphi for two photons
				float _deleta = pho1.Eta() - pho2.Eta(); 
				float _delphi = pho1.Phi() - pho2.Phi();

				if (_delphi > M_PI) {_delphi -= 2 * M_PI;}       
				if (_delphi < -1 * M_PI) {_delphi += 2 * M_PI;}	

				
				// getting the 4-vector for eta-meson
				TLorentzVector mes_lv;

        mes_lv = pho1 + pho2;
								
				pairInvMassTotal->Fill(mes_lv.M());
					
				pairInvMassPtEta->Fill(mes_lv.M(), mes_lv.Pt(), mes_lv.Eta());
				pairpTDelPhiDelEta->Fill(mes_lv.Pt(), _delphi, _deleta);
				
				TVector3 _pho1 = pho1.Vect();
				TVector3 _pho2 = pho2.Vect();

				DelR_pairpT_f->Fill(pho1.DeltaR(pho2), mes_lv.Pt());			
	  		 
       } // inner lopp
     } // outer loop
   		
	}
}


//______________________________________________________________________________..
// getting background (using inv mass histogram) for pi0/eta meson using Ttree
// this code is inspired (mostly based on) from Loop function that creates fore-ground

void CaloCalibEmc_eta::Loop_background_event_mixing(int nevts, TString _filename, TTree * intree)
{
  std::cout << "creating ""Combinatorial Background"" based on Event Mixing. This will create foreground as well." << std::endl;

	// list of parameters used in the loop
	float pt1_cut = 1.0; // first photon cut
	float pt2_cut = 0.5; // second photon cut
	float alpha_cut = 0.6; // energy assymetry cut
	float delR_cut = 1.1; // cone angle cut
	//float chi2_cut = 4.0; // cluster cut (this removes clusters that do not represents electron/photon)
	
	//int max_nCluster_cut = 650; // this will control the most central event (remove events with large number of clusters)
	//int min_nCluster_cut = 120; // min cluster cut to ignore peripheral event
		
	float zver_cut = 30.0; // ignore any event if their z-pos of global vertex is outside of this value (+-)
	float diff_zver_cut = 5.0; // ignore any event if their diff of z-pos of global vertex is outside of this value (+-)
	int bckgnd_evnts = 50; //number of events to go higher (and lower) from current event number (required for background reconstruction)

	
	// different centrality cuts (diff in cent, max cent and min cent in order)
	int cent_range = 10;
	int max_cent_cut = 100;
	int min_cent_cut = 60;
	float photon_prob = 0.05;  // another form of clusters cut of electron/hadron
	

	// making vectors that has TLorentz vectors
	std::vector<TLorentzVector> savClusLV1;
	std::vector<TLorentzVector> savClusLV2;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTree * t1 = intree;
  if (!intree)
  {
	 TFile *f = new TFile(_filename);
   t1 = (TTree *) f->Get("_eventTree");
  }

	// Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_nCentrality", &_nCentrality);
  //t1->SetBranchAddress("_clusterIDs", _clusterIDs);
	t1->SetBranchAddress("_vertex", _vertex); // contains co-ordinates of global-vertex	
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);
  t1->SetBranchAddress("_photonProb",&_photonProb);


  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts) nevts2 = nEntries;
	

  for (int i1 = 0; i1 < nevts2; i1++) // iterating for each/every event(s)
  {
    // load the i1_th instance of the TTree
    t1->GetEntry(i1);

    if (i1 % 1000 == 0) {std::cout << "event nunber = " << i1 << std::endl;}
		
		// ncluster cuts (this can always be treated similar to centrality cut)
		int nClusters1 = _nClusters;
/*		if ((nClusters1 > max_nCluster_cut) || (nClusters1 < min_nCluster_cut)) continue; // will adjust this to look either peripheral or central event
*/
		float vertex_z1 = _vertex[2]; // co-ordinates are 0, 1, 2 for "x" "y" and "z" respectively
		if (abs(vertex_z1) > zver_cut) continue;

		
    float Cent1 = _nCentrality; // centraity of the event
		if ((Cent1 < min_cent_cut) || (Cent1 > max_cent_cut)) continue; // centrality cut
		

		// saving all cluster in vector of TLorentz Vector event #1
		savClusLV1.clear();
		for (int i = 0; i < nClusters1; i++) // for each cluster in event
    {
			if (_clusterPts[i] < std::min(pt1_cut, pt2_cut)) continue;
			//if (_clusterChi2[i] > chi2_cut) continue; 
			if(_photonProb[i] < photon_prob) continue;


		  TLorentzVector LV_Clus;

		  LV_Clus.SetPtEtaPhiE(_clusterPts[i], _clusterEtas[i], _clusterPhis[i], _clusterEnergies[i]);
		  savClusLV1.push_back(LV_Clus);
		 }

		for (int i2 = i1 - bckgnd_evnts; i2 < i1 + bckgnd_evnts; i2++) // we take some events up and down each wrt to event of our concern
			{
				if (i2 < 0) continue; // discard negative event number (which does not exist at all)
				if (i2 > nevts2) continue; // when max event number is reached, we exit our iteration
	
				t1->GetEntry(i2);

				int nClusters2 = _nClusters;
			/*	if ((nClusters2 > max_nCluster_cut) || (nClusters2 < min_nCluster_cut)) continue;
                         */
				float vertex_z2 = _vertex[2]; // co-ordinates are 0, 1, 2
				if (abs(vertex_z2) > zver_cut) continue;

				if (abs(vertex_z1 - vertex_z2) > diff_zver_cut) continue; // ignore those events if their z-vertex is far from 5 units
				
				
				float Cent2 = _nCentrality; // centraity of the event
				if ((Cent2 < min_cent_cut) || (Cent2 > max_cent_cut)) continue; // centrality cut
				
				if (abs(Cent1-Cent2) > cent_range) continue; // we use two events with almost same centrality
				

				// saving all cluster in vector of TLorentz Vector event #2
				savClusLV2.clear();
				for (int j = 0; j < nClusters2; j++) // for each cluster in event
				 {
					if (_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
					//if (_clusterChi2[j] > chi2_cut) continue; 
					if(_photonProb[j] < photon_prob) continue;

					TLorentzVector LV_Clus2;

					LV_Clus2.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
					savClusLV2.push_back(LV_Clus2);
				 }

			
				for (size_t iCs = 0; iCs < savClusLV1.size(); iCs++) // outer loop to iterate clusters
					{
						TLorentzVector& pho1 = savClusLV1[iCs];

						if (fabs(pho1.Pt()) < pt1_cut)	continue; // first P_T cut
						
						for (size_t jCs = iCs+1; jCs < savClusLV2.size(); jCs++) // inner loop to iterate clusters
						 {
								TLorentzVector& pho2 = savClusLV2[jCs];

								if (fabs(pho2.Pt()) < pt2_cut) continue; //second P_T cut
								
								alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
								if (alphaCut > alpha_cut) continue;	
								
								//if ((pho1.DeltaR(pho2) > delR_cut) || (pho1.DeltaR(pho2) < 0.05))  continue;
								if (pho1.DeltaR(pho2) > delR_cut)  continue;

								TLorentzVector mesonlv = pho1 + pho2;

								//if ((mesonlv.M() < 0.05) || (mesonlv.Pt() < 1.0)) continue;

								double _delphi = pho1.Phi() - pho2.Phi();
								double _deleta = pho1.Eta() - pho2.Eta();
								double _pairpT = mesonlv.Pt();
								double _delR = pho1.DeltaR(pho2);

								if (_delphi > M_PI) {_delphi -= 2 * M_PI;}       
								if (_delphi < -1 * M_PI) {_delphi += 2 * M_PI;}
								
								if (i1 == i2) // histogram for foreground
								 {

									pairInvMassTotal->Fill(mesonlv.M());
									pairInvMassPtEta->Fill(mesonlv.M(), mesonlv.Pt(), mesonlv.Eta());
									pairInvMassPtdelR->Fill(mesonlv.M(), mesonlv.Pt(), _delR);

									pairpTDelPhiDelEta->Fill(mesonlv.Pt(), _delphi, _deleta); // for foreground
									DelR_pairpT_f->Fill(_delR, _pairpT);

		 						 }

								else if (i1 != i2) // histogram for background
								 {

								 DelR_pairpT_b->Fill(_delR, _pairpT);

								 pairInvMassTotalBkgd->Fill(mesonlv.M());
						 		 pairInvMassPtEtaBkgd->Fill(mesonlv.M(), mesonlv.Pt(), mesonlv.Eta());
						 		 pairInvMassPtdelRBkgd->Fill(mesonlv.M(), mesonlv.Pt(), _delR);
								 pairpTDelPhiDelEtaBkgd->Fill(_pairpT, _delphi, _deleta);	// for background

								}
							else {std::cout << "Something is not good for event mixing function." << std::endl;} 
									
						} // inner loop (cluster)
				 } // outer loop (cluster)
		 } // second-set of events (i2) (when i1 == i2, it makes foreground)
	}// iterating over all events (first set :: i1)
}
