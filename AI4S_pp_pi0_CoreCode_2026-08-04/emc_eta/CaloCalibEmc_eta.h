// Tell emacs that this is a C++ source
//  -*- C++ -*-
#ifndef CALOCALIBEMC_ETA_H
#define CALOCALIBEMC_ETA_H

#include <fun4all/SubsysReco.h>
#include <string>

class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TH1;
class TNtuple;
class TTree;
class TString;
class TCanvas;

class CaloCalibEmc_eta : public SubsysReco
{
 public:
  CaloCalibEmc_eta(const std::string &name = "CaloCalibEmc_eta",
                   const std::string &fnm = "outJF");

  ~CaloCalibEmc_eta() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void Loop(int nevts, TString _filename, TTree *intree);
  void Loop_background_event_mixing(int nevts, TString _filename, TTree *intree);

  // Setters to support external configuration
 void SetScaledowns(int scaledowns[])
  {
    for (int i = 0; i < 64; ++i)
    {
    m_scaledowns[i] = scaledowns[i];
   }
  }

  void setRunNumber(int runnum)
  {
    _runnumber = runnum;
  }

 private:
  int m_ievent = 0;
  std::string m_Filename;
  TFile *cal_output = nullptr;
  std::string _caloname = "CEMC";
  bool m_UseTowerInfo;

  // Histograms
  TH1F *pairInvMassTotal = nullptr;
  TH1F *pairInvMassTotalBkgd = nullptr;
  TH3F *pairInvMassPtEta = nullptr;
  TH3F *pairInvMassPtEtaBkgd = nullptr;
  TH3F *pairpTDelPhiDelEta = nullptr;
  TH3F *pairpTDelPhiDelEtaBkgd = nullptr;
  TH3F *pairInvMassPtdelR = nullptr;
  TH3F *pairInvMassPtdelRBkgd = nullptr;
  TH2F *DelR_pairpT_f = nullptr;
  TH2F *DelR_pairpT_b = nullptr;
  TH2F *h_etaphi_clus = nullptr;

  // Tree and event-level variables
  TTree *_eventTree = nullptr;
  int _eventNumber = -1;
  int _nClusters = -1;
  float _nCentrality = -1.0;
  float _vertex[3] = {0};
  float _clusterIDs[10000] = {0};
  float _clusterEnergies[10000] = {0};
  float _clusterPts[10000] = {0};
  float _clusterEtas[10000] = {0};
  float _clusterPhis[10000] = {0};
  float _clusterChi2[10000] = {0};
  int maxTowerEta = -1;
  int maxTowerPhi = -1;
  int _maxTowerEtas[10000] = {0};
  int _maxTowerPhis[10000] = {0};

  float alphaCut = -1.;
  int m_cent_nclus_cut;
  TFile *f_temp = nullptr;

  //Added variables
  int _runnumber = -1;
  float _photonProb[10000] = {0};
  //float _photonProb[MAX_SIZE] = {0};
  bool ScaledTriggerBit[64] = {false};
  bool LiveTriggerBit[64] = {false};
  int m_scaledowns[64] = {0};
 uint64_t count_raw[64] = {0};
  uint64_t count_live[64] = {0};
  uint64_t count_scaled[64] = {0};
  int evt_isNotMB = 0;
  int evt_isMB = 0;
   //single particle mc for ttree
  float truth_vertex[3] = {0.0};
  //float truth_vz=-999;
  //float truth_vx=-999;
  //float truth_vy=-999;
};

#endif  // CALOCALIBEMC_ETA_H

