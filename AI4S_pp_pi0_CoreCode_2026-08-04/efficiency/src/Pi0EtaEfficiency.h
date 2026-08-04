// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PI0ETAEFFICIENCY_H
#define PI0ETAEFFICIENCY_H
//#pragma once

#include <fun4all/SubsysReco.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>
#include <TClonesArray.h>


class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TH3I;
class TH3D;
class TF1;
class TH1;
class TNtuple;
class TTree;
class TString;
class TCanvas;
class PHG4Particle;
class eventTree;
class truthEventTree;


class Pi0EtaEfficiency : public SubsysReco
{
 public:

  Pi0EtaEfficiency(const std::string &name = "Pi0EtaEfficiency", const char *outfile = "out_dst_calo.root");

  ~Pi0EtaEfficiency() override;

  int Init(PHCompositeNode* topNode) override;

  void InitOutputFile();
  
  void InitTree();

  void InitTruthTree();

  void InitHistos();

  void InitTruthHistos();

  int process_event(PHCompositeNode *topNode) override;

  void ProcessFillTruthParticle(PHCompositeNode *topNode, float vz);

  double GetShiftedEta(PHG4Particle *particle);

  int End(PHCompositeNode *topNode) override;

  void End();

 /* void Loop(TTree *eventTree, float leadPhotonPtCut, float subPhotonPtCut, float asymmetryCut, float deltaRCut, float clusterChi2Cut, float zVertexCut, int minClusterCut, int maxClusterCut, int firstEntry, int lastEntry);

  void Loop_background_event_mixing(TTree *eventTree, float leadPhotonPtCut, float subPhotonPtCut, float asymmetryCut, float deltaRCut, float clusterChi2Cut, float zVertexCut, int nMixedEvents, int minClusterCut, int maxClusterCut); */
	
  void SetScaledowns(int scaledowns[])
  {
    for(int i=0; i<64;i++)
      {
	m_scaledowns[i] = scaledowns[i];
      }
  };

  void setRunNumber(int runnum)
  {
    _runnumber = runnum;
  };

  void setFlag(bool data, bool hijing)
  {
    isDATA = data;
    isHIJING = hijing;
  };

  // ---- Collision system: p+p (PYTHIA) vs nucleus-nucleus (HIJING/Au+Au) ----
  // PHYSICS: centrality, impact parameter (bimp), Npart and Ncoll are Glauber
  // quantities that only exist for NUCLEUS-NUCLEUS collisions. PYTHIA p+p has
  // one nucleon per side -> NO impact-parameter geometry and NO centrality;
  // every minimum-bias p+p event is the same "class". When p+p mode is on we
  // (a) never dereference the (possibly absent/meaningless) centrality node and
  // (b) assign every event to ONE inclusive centrality bin, so the efficiency
  // is reported inclusively vs pT. Default = true (this is the p+p space);
  // call setIsPP(false) to restore the original Au+Au centrality behavior.
  //void setIsPP(bool v) { m_isPP = v; };

  void set_zvtx_cut(int z)
  {
    myzvtx_cut = z;
  };

  // ---- Runtime meson selection (Run-24 p+p pi0/eta efficiency) ----
  // Retargets the truth-particle gate (PDGPID), the PDG mass used for
  // truth bookkeeping, the reco diphoton mass acceptance window, and the
  // invariant-mass fit range/seed. Call ONCE from the Fun4All macro before
  // run() so pi0 (111) and eta (221) produce independent outputs (decision
  // #4/#7: 2-gamma channel only, species kept separate). Defaults below
  // preserve the original pi0 behavior if setMeson() is never called.
  void setMeson(int pdgid)
  {
    if (pdgid == 221)  // eta -> gamma gamma
    {
      PDGPID       = 221;
      PDG_MASS     = 0.547862f;  // PDG eta mass (GeV)
      data_eta_mass= 0.560f;     // expected reco eta peak (tune to data/MC)
      m_massWinLo  = 0.35f;      // reco diphoton mass acceptance window
      m_massWinHi  = 0.75f;
      m_fitLo      = 0.35;       // gaus+pol3 fit range for eta peak
      m_fitHi      = 0.75;
      m_fitMeanLo  = 0.50;       // par1 (mean) limits, widened for eta
      m_fitMeanHi  = 0.60;
      m_fitMeanSeed= 0.548;
    }
    else  // default / pi0 -> gamma gamma (pdgid 111)
    {
      PDGPID       = 111;
      PDG_MASS     = 0.1349768f;
      data_eta_mass= 0.155f;
      m_massWinLo  = 0.05f;
      m_massWinHi  = 0.95f;
      m_fitLo      = 0.065;
      m_fitHi      = 0.45;
      m_fitMeanLo  = 0.11;
      m_fitMeanHi  = 0.185;
      m_fitMeanSeed= 0.135;
    }
  };

  // Toggle the MC reweight (MC_PHENIX_reweight.root "weights") ON/OFF
  // (decision #6: weighted by default, single flag to run unweighted).
  void setUseWeights(bool u) { m_useWeights = u; };

  // ---- Signal isolation by embed id (CRITICAL for p+p PYTHIA embedding) ----
  // The efficiency truth denominator must contain ONLY the embedded gun signal
  // mesons, never the mesons of the underlying event. sPHENIX G4_Input.C assigns
  // embed ids as: pileup/background PYTHIA8 = 1, and the FIRST embedded signal
  // generator = 2 (when Input::EMBED is on; it is 1 when the gun runs
  // standalone with no background). The original Au+Au code isolated the signal
  // with isEmbeded()>0, which is WRONG here because the pythia8_Detroit
  // background is itself embed_id 1 (>0) -> its soft mesons flooded the truth
  // histograms and the denominator. We now require isEmbeded(track)==signal id.
  //   setSignalEmbedId(2) : embedded production (default)  -> gun signal only
  //   setSignalEmbedId(1) : standalone gun (USE_EMBED=false smoke test)
  void setSignalEmbedId(int id) { m_signalEmbedId = id; };

  void recoTruthMatch(eventTree &recoTree, truthEventTree &truthTree, double matchDeltaRCut, double asymmetryCut, double clusterChi2Cut, double diphotonDeltaRCut, double clusterEnergyCut, double leadPhotonPtCut, double subPhotonPtCut, double diphotonPtCut, TH1F *ptReweight);

  //int getCentralityBin(double);

  void setBufferSize(float z, int clo, int chi);

  void makeBuffer(bool mb)
  {
    makeBuff = mb;
  }

  struct mixedEventStruct
  {
    int eventNumber;
    std::vector<TLorentzVector> clusters;
  };
    
 private:

  std::vector<std::vector<std::vector<mixedEventStruct>>> eventBuffer;

  
  int myzvtx_cut = 150;  
  int evt_isNotMB = 0;
  int evt_isMB = 0;
  int noVtx = 0; 
  int ievent = 0;
  //int centlo=0;
 // int centhi=100;
  int min_nCluster_cut=0; 
  int max_nCluster_cut= 1000; 

  
  
  bool isDATA = true;
  bool isHIJING = false;
  bool makeBuff = false;
  //bool m_isPP = true;   // PYTHIA p+p => no centrality (see setIsPP / getCentralityBin)

  float photonEtaMax = 1.0;

  std::string _caloname = "CEMC";
  std::string outFileName;
  TFile *anaOutFile = nullptr;

  float vx = -1;
  float vy = -1;
  float vz = -1;
 
  //Histo list
  // Diphoton (gamma-gamma) invariant-mass histograms.
  // "fg" = same-event foreground (real pairs); "bg" = mixed-event background.
  TH2F *h_diphotonFgMassVsPt = nullptr;        // foreground: m_gg vs pair pT
  TH2F *h_diphotonBgMassVsPt = nullptr;        // mixed-event background: m_gg vs pair pT
  TH2F *h_diphotonBgMassVsPtWeighted = nullptr;// mixed-event background, mixing-weighted

  TH1F *h_diphotonFgMass = nullptr;            // foreground m_gg (1D)
  TH1F *h_diphotonBgMass = nullptr;            // mixed-event background m_gg (1D)


 // TH1F *h_eventCentrality = nullptr;
  TH1F *h_clusterEnergy = nullptr;
  TH1F *h_zVertex = nullptr;

  TH1F *h_diphotonFgDeltaR = nullptr;          // foreground photon-pair deltaR
  TH1F *h_diphotonBgDeltaR = nullptr;          // background photon-pair deltaR

  TH1F *h_diphotonFgAsymmetry = nullptr;       // foreground energy asymmetry alpha
  TH1F *h_diphotonBgAsymmetry = nullptr;       // background energy asymmetry alpha

  TH2F *h_towerGoodMap = nullptr;              // ieta-iphi map of isGood towers

  TH2F *h_diphotonFgDeltaRVsPt = nullptr;      // foreground deltaR vs pair pT
  TH2F *h_diphotonFgDeltaRVsPtMassSideband = nullptr; // same, mass-sideband selection
  TH2F *h_diphotonBgDeltaRVsPt = nullptr;      // background deltaR vs pair pT

  TH3F *h_diphotonFgMassPtDeltaR = nullptr;    // foreground m_gg vs pair pT vs deltaR
  TH3F *h_diphotonBgMassPtDeltaR = nullptr;    // background m_gg vs pair pT vs deltaR


  TH1F *h_usedEventCount = nullptr;            // number of events used (for cross-section scaling)
  TH1 *h_scaledTriggerBit = nullptr;
  TH1 *h2_scaledTriggerBit = nullptr;
  //TH2F *h_zVertexVsCentrality = nullptr;
  TH2F *h_towerHotMap = nullptr;               // ieta-iphi map of isHot towers

  //TH1F * bkg_eff = nullptr;
  //TH3F *ieta_iphi_energy = nullptr;
  //TH2F * diPhotonMap = nullptr;
  //TH2F *clusterMap = nullptr;

  static const int MAX_SIZE = 3000;
  
  //TTree info
  TTree *_eventTree = nullptr;     
  int _nClusters = 0;
  int maxTowerEta = -1;
  int maxTowerPhi = -1;
  int _runnumber = -1; 
  //float _Centrality = -1.0;
  float _vertex[3] = {0.0};
  float _clusterEnergies[MAX_SIZE] = {0};
  float _clusterPts[MAX_SIZE] = {0};
  float _clusterEtas[MAX_SIZE] = {0};
  float _clusterPhis[MAX_SIZE] = {0};
  float _clusterChi2[MAX_SIZE] = {0};
  float _photonProb[MAX_SIZE] = {0};
  float _towerEnergy[MAX_SIZE] = {0};
  int _maxTowerEtas[MAX_SIZE] = {0};
  int _maxTowerPhis[MAX_SIZE] = {0};
  bool ScaledTriggerBit[64];
  bool LiveTriggerBit[64];
  long long int count_raw[64];
  long long int count_live[64];
  long long int count_scaled[64];
  int m_scaledowns[64];


  // ===== All truth related stuff ===== //
  
  int PDGPID = 111;//neutral meson
  float PDG_MASS = 0.1349768;   // PDG π⁰ mass in GeV
  float data_eta_mass = 0.155;  // your expected reconstructed π⁰ mass (adjust if needed)

  // Species-driven tunables (set by setMeson(); pi0 defaults shown).
  // m_useWeights toggles the MC reweight; see setUseWeights().
  bool  m_useWeights = true;
  int   m_signalEmbedId = 2;    // embed id of the gun signal (see setSignalEmbedId)
  float m_massWinLo  = 0.05f;   // reco diphoton mass acceptance window (GeV)
  float m_massWinHi  = 0.95f;
  double m_fitLo     = 0.065;   // gaus(0)+pol3(3) fit range (GeV)
  double m_fitHi     = 0.45;
  double m_fitMeanLo = 0.11;    // par1 (peak mean) limits + seed
  double m_fitMeanHi = 0.185;
  double m_fitMeanSeed = 0.135;
  
      
  float truth_vx;
  float truth_vy;
  float truth_vz;
  
  int recoveredPhotons;
  int not_recoveredPhotons;
  int recoveredPi0;
  int recoveredPiPlus;
  int recoveredPiMinus;
  int noRange;

  static const int INIT = 0;
  
  TTree *_truthEventTree = nullptr;

  bool _etaEmbedStatus[50000000];

  TClonesArray *acceptedEtaVec = nullptr;
  TClonesArray *acceptedPhoton1Vec = nullptr;
  TClonesArray *acceptedPhoton2Vec = nullptr;

  TClonesArray *rejectedEtaVec = nullptr;
  TClonesArray *rejectedPhoton1Vec = nullptr;
  TClonesArray *rejectedPhoton2Vec = nullptr;

  /*
  TClonesArray *acceptedPi0Vec = nullptr;
  TClonesArray *acceptedPhoton1VecPi0 = nullptr;
  TClonesArray *acceptedPhoton2VecPi0 = nullptr;

  TClonesArray *rejectedPi0Vec = nullptr;
  TClonesArray *rejectedPhoton1VecPi0 = nullptr;
  TClonesArray *rejectedPhoton2VecPi0 = nullptr;
  */
  

  int _NrejectedEtas = 0;
  int _nEtas = 0;
  int _allEtas = 0;
  
 // int _npart = INIT;              //from event_header
 // int _ncoll = INIT;              //from event_header
  //loat _bimp = INIT;             //from event_header
  
  //int _truthCentrality;          
 // float _centImpactParam;
  float _truthZvtx[3] = {INIT};
  

   //Truth histo list
  TH1F *h_truthZVertex = nullptr;
  //TH1F *h_truthEventCentrality = nullptr;

  // Truth meson pT spectra. "accepted" = both decay photons inside |eta|<photonEtaMax.
  TH1F *h_truthMesonPtAccepted = nullptr; // truth meson pT when BOTH photons are in acceptance
  TH1F *h_truthMesonPtAll = nullptr;      // truth meson pT for all 2-gamma decays (any photon position)
  TH2F *h_truthMesonPtVsEtaAccepted = nullptr; // meson pT vs pseudorapidity, both photons accepted
  TH1F *h_truthDecayMode = nullptr;       // decay-mode tally
  TH3F *h_truthLostPhotonZvtxEtaPt = nullptr;  // (zvtx, eta, meson pT) when a photon is out of acceptance
  TH2F *h_truthMesonPtVsEtaAll = nullptr; // meson pT vs pseudorapidity, all photon positions

  TH1F *h_truthMesonPtDist = nullptr;     // raw truth meson pT distribution

  // truth-reco matching histos
  TH2F *h_recoMassVsRecoPt = nullptr;       // reco m_gg vs reco meson pT (weighted)
  TH2F *h_truthMassVsTruthPt = nullptr;     // truth m_gg vs truth meson pT (weighted)
  TH2F *h_recoMassVsRecoPtUnweighted = nullptr;
  TH2F *h_truthMassVsTruthPtUnweighted = nullptr;
  TH1F *h_truthDeltaR = nullptr;            // truth photon-pair deltaR
  TH1F *h_recoDeltaR = nullptr;             // reco photon-pair deltaR
  TH1F *h_truthPhotonEnergy = nullptr;
  TH1F *h_recoClusterEnergy = nullptr;
 // TH2F *h_truthPtVsCentrality = nullptr;       // efficiency DENOMINATOR (weighted)
 // TH2F *h_recoPtVsCentrality = nullptr;        // efficiency NUMERATOR (weighted)
  //TH2F *h_truthPtVsCentralityUnweighted = nullptr;
  //TH2F *h_recoPtVsCentralityUnweighted = nullptr;

  // ---- p+p efficiency vs pT (centrality-independent) ----
  // PYTHIA p+p has no centrality, so the efficiency is a plain 1D function of
  // truth meson pT: eff(pT) = N_reco(pT) / N_truth(pT). These 1D denominators
  // (truth) and numerators (reco-matched) are filled directly, with NO
  // centrality axis, and divided in End of recoTruthMatch to give the final
  // h_efficiencyVsPt. Weighted + unweighted variants are produced.
  TH1F *h_truthPt = nullptr;                    // efficiency DENOMINATOR vs pT (weighted)
  TH1F *h_recoPt = nullptr;                     // efficiency NUMERATOR   vs pT (weighted)
  TH1F *h_truthPtUnweighted = nullptr;          // DENOMINATOR vs pT (unweighted)
  TH1F *h_recoPtUnweighted = nullptr;           // NUMERATOR   vs pT (unweighted)
  TH1F *h_efficiencyVsPt = nullptr;             // eff(pT) = reco/truth (weighted)
  TH1F *h_efficiencyVsPtUnweighted = nullptr;   // eff(pT) = reco/truth (unweighted)
 // TH1F *h_centralityFromImpactParam = nullptr;
  TH1F *h_truthPhotonAsymmetry = nullptr;
  TH1F *h_recoPhotonAsymmetry = nullptr;
  TH1F *h_deltaRResidual = nullptr;            // reco - truth deltaR
  TH1F *h_energyResidual = nullptr;            // reco cluster E - truth photon E
  TH2F *h_truthPtVsDeltaR = nullptr;
  TH2F *h_truthPhotonEtaPhi = nullptr;
  TH2F *h_recoClusterEtaPhi = nullptr;
 // TH3D *h_truthPtCentralityMass = nullptr;
 // TH3D *h_truthPtCentralityMassMatched = nullptr;

  
  
};

#endif // PI0ETAEFFICIENCY_H
