#ifndef MACRO_FUN4ALLG4PYTHIAEFF_C
#define MACRO_FUN4ALLG4PYTHIAEFF_C

// =====================================================================
//  Fun4All_G4_pythia_eff.C
//  Run-24 p+p  pi0 / eta  two-photon reconstruction-EFFICIENCY driver.
//
//  This is the p+p analogue of Spring_Pi0_eff/Fun4All_G4_eta_embed.C and
//  follows that reference's structure directly. A single pi0/eta is fired
//  by the SimpleEventGenerator, transported through Geant4, and OVERLAID
//  (embedded) onto a pre-simulated MB pythia8_Detroit background G4Hits
//  DST. The embedded signal gets embed_id >= 2 while the background keeps
//  embed_id 0, so Pi0EtaEfficiency's isEmbeded()>0 cut isolates the signal.
//
//  One background DST is embedded per job, indexed by fileNum (exactly as
//  the Au+Au reference indexes its sHijing DST by fileNum).
//
//  USE_EMBED = true  (default) : the real efficiency configuration --
//                                embed the gun into one background DST.
//  USE_EMBED = false           : standalone smoke test (no background);
//                                the gun makes its own vertex so the whole
//                                chain can be exercised with NO input file.
// =====================================================================

#include <GlobalVariables.C>

// The analysis module header (forked p+p efficiency copy):
#include "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/src/Pi0EtaEfficiency.h"

#include <DisplayOn.C>
#include "/sphenix/user/plewis3323/Fall_2025/pi0_Eff/macro/G4Setup_sPHENIX.C"
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>
#include <Trkr_QA.C>
#include <Trkr_Diagnostics.C>

#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <Calo_Calib.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcalo_reco.so)
// Load the forked p+p efficiency module from its dedicated install prefix:
R__LOAD_LIBRARY(/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/install/lib/libPi0EtaEfficiency.so)

// Directory holding the pre-simulated MB pythia8_Detroit background G4Hits
// DSTs (run 28) the signal is embedded into. This is the physical lustre
// production location resolved from the FileCatalog LFNs (the SIM_DST_SPACE
// g4hits.list holds only bare LFNs, not paths, so it can't be fileopen'd
// directly). One segment per job, indexed by fileNum.
static const std::string BKG_DIR =
    "/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal/g4hits/run0028/detroit/";

// ---------------------------------------------------------------------
//  Arguments:
//    fileNum    : zero-padded index; selects the background DST to embed
//                 into AND names the output file.
//    pdgid      : 111 = pi0, 221 = eta   (threads to setMeson()).
//    useWeights : true = MC pT reweight ON, false = unweighted.
//    nEvents    : number of events to process (<0 build-only; 0 = all).
//    USE_EMBED  : true = embed the gun into one background DST (real run);
//                 false = standalone gun, no background (smoke test).
//    outDir     : directory for the analysis output ROOT file.
// ---------------------------------------------------------------------
int Fun4All_G4_pythia_eff(int fileNum = 0,
                          int pdgid = 111,
                          bool useWeights = false,
                          int nEvents = 300,
                          bool USE_EMBED = true,
                          const std::string &outDir = "/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/output/pi0_Prod/Condor_Out")
{
  // zero-pad the file index for tidy, sortable filenames
  std::string fileNumber = std::to_string(fileNum);
  fileNumber = std::string(6 - fileNumber.length(), '0') + fileNumber;

  std::string species = (pdgid == 221) ? "eta" : "pi0";

  // Background DST this job embeds into (one file per job, indexed by fileNum).
  const std::string embed_input_file =
      "G4Hits_pythia8_Detroit-0000000028-" + fileNumber + ".root";

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  Enable::CDB = true;

  // ---- Conditions database ----
  // MDC2 is the simulation tag the pythia8_Detroit background was produced
  // with (same tag family as the Au+Au reference). The run/segment come from
  // the background DST filename so the EMCal calibration matches the sample.
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");

  int runnumber = 0;
  int segment = 0;
  if (USE_EMBED)
  {
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(embed_input_file);
    runnumber = runseg.first;
    segment = runseg.second;
  }
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  if (runnumber != 0)
  {
    rc->set_IntFlag("RUNNUMBER", runnumber);
    Fun4AllSyncManager *syncman = se->getSyncManager();
    syncman->SegmentNumber(segment);
  }
  else
  {
    // MC standalone: set RUNNUMBER explicitly so subsystems that query it
    // (e.g. TPC SetDefaultParameters) do not warn about an unknown flag.
    rc->set_IntFlag("RUNNUMBER", runnumber);
  }

  // ---------------------------------------------------------------
  //  INPUT: SimpleEventGenerator gun, optionally embedded into a
  //         pre-simulated background DST.
  // ---------------------------------------------------------------
  Input::VERBOSITY = 0;

  // Input::EMBED must be set BEFORE InputInit(): InputInit() bumps the signal
  // EmbedId to 2 and enables vertex reuse when EMBED is on.
  Input::EMBED = USE_EMBED;

  if (USE_EMBED)
  {
    // One background DST per job (fileopen path, Au+Au reference style).
    // REPEAT=false: play the background once -- it terminates the event loop,
    // so nEvents=0 safely means "process every event in this background DST".
    INPUTEMBED::filename[0] = BKG_DIR + embed_input_file;
    INPUTEMBED::REPEAT = false;
    std::cout << "[Fun4All_G4_pythia_eff] EMBED ON -> background DST: "
              << BKG_DIR + embed_input_file << std::endl;
  }

  Input::SIMPLE = true;
  Input::HEPMC = false;

  InputInit();

  INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(species, 3);
  if (Input::EMBED)
  {
    // EMBEDDING: place the signal meson at the BACKGROUND event's primary
    // vertex (reuse it) so signal and background share one z. Without this the
    // reco would be nonsense.
    INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
  }
  else
  {
    // STANDALONE smoke test: no background, so the gun makes its own vertex.
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(
        PHG4SimpleEventGenerator::Gaus,
        PHG4SimpleEventGenerator::Gaus,
        PHG4SimpleEventGenerator::Gaus);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);
  }
  INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.1, 1.1);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(1.0, 15.);

  InputRegister();

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  //======================
  // What to run (detectors needed for an EMCal cluster efficiency study)
  //======================
  Enable::MBD = true;
  Enable::MBDRECO = Enable::MBD && true;
  Enable::PIPE = true;
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;

  Enable::CEMC = true;                      // the calorimeter that matters
  Enable::CEMC_ABSORBER = true;
  Enable::CEMC_CELL = Enable::CEMC && true;
  Enable::CEMC_TOWER = Enable::CEMC_CELL && true;
  Enable::CEMC_CLUSTER = Enable::CEMC_TOWER && true;

  Enable::HCALIN = true;
  Enable::MAGNET = true;
  Enable::HCALOUT = true;
  Enable::EPD = true;
  Enable::BEAMLINE = true;
  Enable::ZDC = true;
  Enable::PLUGDOOR = true;
  Enable::PLUGDOOR_ABSORBER = true;

  // Centrality is an Au+Au (Glauber) concept; kept enabled only to keep the
  // reco chain identical to the reference. The analysis runs in p+p mode
  // (setIsPP(true) below) and does not depend on the centrality node.
  Enable::CENTRALITY = true;

  Enable::BLACKHOLE = true;
  Enable::BLACKHOLE_SAVEHITS = false;

  G4Init();
  if (!Input::READHITS)
  {
    G4Setup();
  }

  // --- Calorimeter cell -> tower -> cluster reconstruction ---
  if ((Enable::MBD && Enable::MBDRECO) || Enable::MBDFAKE) Mbd_Reco();
  if (Enable::CEMC_CELL) CEMC_Cells();
  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  // --- Global vertex + centrality ---
  if (Enable::GLOBAL_RECO) Global_Reco();
  else if (Enable::GLOBAL_FASTSIM) Global_FastSim();
  if (Enable::CENTRALITY) Centrality();

  InputManagers();

  // nEvents < 0 : build everything but process no events (dry run).
  if (nEvents < 0) return 0;
  // A standalone gun (no embed) has no file to exhaust, so nEvents=0 would run
  // forever. With embedding the (REPEAT=false) background DST terminates the
  // loop, so nEvents=0 is fine there and means "all events in the DST".
  if (nEvents == 0 && !USE_EMBED)
  {
    std::cout << "[Fun4All_G4_pythia_eff] nEvents=0 with a standalone gun would "
                 "run forever; refusing. Pass a positive nEvents." << std::endl;
    return 0;
  }

  // ---------------------------------------------------------------
  //  ANALYSIS MODULE
  // ---------------------------------------------------------------
  std::string outFile = outDir + "/eff_" + species + "_run24pp_" + fileNumber + ".root";

  Pi0EtaEfficiency *ana = new Pi0EtaEfficiency(species.c_str(), outFile.c_str());
  // setFlag(isDATA, isHIJING): MC truth path requires isHIJING=1 so the truth
  // tree + truth histos (the efficiency DENOMINATOR) get filled.
  ana->setFlag(0, 1);
  //ana->setIsPP(true);             // PYTHIA p+p: no centrality (single inclusive bin)
  ana->set_zvtx_cut(30);          // |z_vtx| < 30 cm, standard fiducial cut
  ana->setMeson(pdgid);           // 111 pi0 / 221 eta -> retargets PDGID, mass window, fit
  ana->setUseWeights(useWeights); // MC pT reweight toggle
  // Isolate ONLY the embedded gun signal in the truth denominator. G4_Input.C
  // gives the gun embed_id 2 when embedding into a background (Input::EMBED),
  // but embed_id 1 when it runs standalone. The pythia8_Detroit background is
  // itself embed_id 1, so we MUST match the gun's exact id (not isEmbeded()>0)
  // or the soft background mesons flood the flat-signal truth spectrum.
  ana->setSignalEmbedId(USE_EMBED ? 2 : 1);
  se->registerSubsystem(ana);

  std::cout << "[Fun4All_G4_pythia_eff] species=" << species
            << " pdgid=" << pdgid
            << " useWeights=" << useWeights
            << " USE_EMBED=" << USE_EMBED
            << " out=" << outFile << std::endl;

  se->run(nEvents);
  se->End();

  std::cout << "All done" << std::endl;
  // NOTE: intentionally NOT calling `delete se`. End() has already run End on
  // every subsystem and written+closed the output ROOT file, so there is no
  // work left to do. `delete se` walks the module/node destructor chain, which
  // segfaults at teardown (some module destructors touch already-freed nodes) --
  // a benign crash AFTER all work is done, but it yields a nonzero exit code
  // that Condor/monitoring flags as a job failure. Terminate immediately with
  // _exit(0) (mode=kFALSE) so no destructor/atexit chain runs at all.
  gSystem->Exit(0, kFALSE);
  return 0;
}
#endif
