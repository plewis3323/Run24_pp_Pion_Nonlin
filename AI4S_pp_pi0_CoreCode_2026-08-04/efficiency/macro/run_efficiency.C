// run_efficiency.C
// Runs Pi0EtaEfficiency::recoTruthMatch on the merged pi0 production file
// and writes the efficiency histograms into runEff/.
//
// Usage:
//   root -b -q 'run_efficiency.C()'
//   root -b -q 'run_efficiency.C("myfile.root","eff_out.root",221)'   // eta

R__LOAD_LIBRARY(/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/install/lib/libPi0EtaEfficiency.so)

#include "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/src/eventTree.h"
#include "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/src/truthEventTree.h"
#include "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/src/Pi0EtaEfficiency.h"

// ---- hardcoded paths ------------------------------------------------------
const TString BASE_DIR  = "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia";
// hadd'ed pi0 production merge (10000 condor jobs)
const TString PI0_MERGE = BASE_DIR + "/output/pi0_Prod/Hadd_Out/eff_pi0_run24pp_Merge1.root";
const TString ETA_MERGE = BASE_DIR + "/output/Eta_Prod/Hadd_Out/eff_eta_run24pp_Merge1.root";
// ---------------------------------------------------------------------------

void run_efficiency(const char* input   = "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/output/Eta_Prod/Hadd_Out/eff_eta_run24pp_Merge1.root",     // "" = use hardcoded merge for the species
                    const char* outfile = "Eta_eff_out.root",
                    int    pdgid       = 221,     // 111 pi0, 221 eta
                    double matchDeltaR = 0.3,
                    double asym        = 0.6,
                    double chi2        = 5.0,
                    double diphotonDR  = 1.1,
                    double leadPt      = 0.8,
                    double subPt       = 0.6,
                    double diphotonPt  = 2.0,
                    bool   useWeights  = true)    // Hagedorn pT reweight
{
  TString species = (pdgid == 221) ? "eta" : "pi0";

  // input file: caller's choice, otherwise the hardcoded production merge
  TString inFile = input;
  if (inFile.Length() == 0) inFile = (pdgid == 221) ? ETA_MERGE : PI0_MERGE;
  cout << "input:  " << inFile << endl;

  // everything gets written into runEff/
  TString resultDir = BASE_DIR + "/runEff";
  gSystem->mkdir(resultDir, kTRUE);
  TString outPath = resultDir + "/" + gSystem->BaseName(outfile);
  cout << "output: " << outPath << endl;

  // chain the trees
  TChain *ce = new TChain("_eventTree");
  TChain *ct = new TChain("_truthEventTree");
  ce->Add(inFile);
  ct->Add(inFile);
  cout << "entries: event=" << ce->GetEntries()
       << " truth=" << ct->GetEntries() << endl;
  if (ce->GetEntries() == 0) { cout << "no events, abort" << endl; return; }

  eventTree      recoTree(ce);
  truthEventTree truthTree(ct);

  Pi0EtaEfficiency *ana = new Pi0EtaEfficiency("species", outPath.Data());
  ana->setMeson(pdgid);

  // Hagedorn pT weights (w(pT) = Hagedorn/flat on the gun range [1,15] GeV)
  TH1F *ptWeights = nullptr;
  if (useWeights) {
    TString wPath = BASE_DIR + "/Main_Reweight_Plan/Analysis/weights/MC_HAGEDORN_reweight_" + species + ".root";
    TFile *fw = TFile::Open(wPath, "READ");
    if (fw && !fw->IsZombie() && fw->Get("weights")) {
      ptWeights = (TH1F*)fw->Get("weights")->Clone("ptWeights");
      ptWeights->SetDirectory(0);
      cout << "weights: " << wPath << endl;
      fw->Close();
    } else {
      cout << "WARNING: no weight file (" << wPath << "), running unweighted" << endl;
      useWeights = false;
    }
  }
  ana->setUseWeights(useWeights);

  TFile *fout = new TFile(outPath, "RECREATE");
  fout->cd();
  if (!ptWeights) {  // recoTruthMatch needs a non-null hist even when unweighted
    ptWeights = new TH1F("ptWeights", "flat weights", 60, 0, 30);
    for (int b = 1; b <= 60; b++) ptWeights->SetBinContent(b, 1.0);
    ptWeights->SetDirectory(0);
  }

  ana->recoTruthMatch(recoTree, truthTree, matchDeltaR, asym, chi2, diphotonDR,
                      0.0, leadPt, subPt, diphotonPt, ptWeights);

  fout->cd();
  fout->Write();

  // quick summary
  TH1F *num = (TH1F*)fout->Get("h_recoPt");
  TH1F *den = (TH1F*)fout->Get("h_truthPt");
  if (num && den && den->Integral() > 0)
    cout << "integrated efficiency = " << num->Integral()/den->Integral()
         << "  (num=" << num->Integral() << " den=" << den->Integral() << ")" << endl;

  TH1F *eff = (TH1F*)fout->Get(useWeights ? "h_efficiencyVsPt" : "h_efficiencyVsPtUnweighted");
  if (eff) {
    cout << "eff vs pT:" << endl;
    for (int b = 1; b <= eff->GetNbinsX(); b++)
      if (eff->GetBinContent(b) > 0)
        printf("  pT ~ %5.2f GeV : %.4f +/- %.4f\n",
               eff->GetBinCenter(b), eff->GetBinContent(b), eff->GetBinError(b));
  }

  fout->Close();
  cout << "wrote " << outPath << endl;
}
