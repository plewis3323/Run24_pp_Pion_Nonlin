// Combined Fun4All Macro (Code 1 + Code 2, adapted for Condor usage)
#ifndef MACRO_FUN4ALLG4SLOPECAL_C
#define MACRO_FUN4ALLG4SLOPECAL_C

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <GlobalVariables.C>
#include <Calo_Calib.C>

// Fun4All framework (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllServer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllInputManager.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllDstInputManager.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllNoSyncDstInputManager.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllUtils.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllRunNodeInputManager.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllHistoManager.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/Fun4AllReturnCodes.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/fun4all/SubsysReco.h"

// Event
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/ffaobjects/EventHeaderv1.h"

// Trigger (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/trigger/TriggerRunInfov1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/trigger/TriggerAnalyzer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/trigger/MinimumBiasInfov1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/trigger/MinimumBiasInfo.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/trigger/MinimumBiasClassifier.h"

// Global vertex (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/globalvertex/GlobalVertexMap.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/globalvertex/GlobalVertexMapv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/globalvertex/MbdVertex.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/globalvertex/MbdVertexMapv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/globalvertex/GlobalVertex.h"

// Tower includes (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawTower.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawTowerContainer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawTowerGeom.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawTowerGeomContainer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawTowerGeomContainer_Cylinderv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoContainerv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfov1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoContainerSimv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoSimv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoContainerv2.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfov2.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoContainerv3.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfov3.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoContainerv4.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfov4.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/TowerInfoDefs.h"

// MBD (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdOut.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdPmtContainer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdPmtContainerV1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdPmtSimContainerV1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdPmtHit.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/mbd/MbdGeom.h"

// Cluster includes (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawCluster.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawClusterv1.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawClusterContainer.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/CaloBase/RawClusterUtility.h"

// phool (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/phool/recoConsts.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/phool/getClass.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/framework/phool/PHCompositeNode.h"

// Centrality MB (absolute paths)
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/centrality/CentralityInfo.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/centrality/CentralityInfov2.h"
#include "/sphenix/user/plewis3323/MainMastersWork/src/offline/packages/centrality/CentralityReco.h"

// My Analysis
#include "/sphenix/user/plewis3323/MainMastersWork/src/calo_emc_pi0_tbt/CaloCalibEmc_eta.h"


// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalibCaloEmc_eta.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libmbd.so) 
R__LOAD_LIBRARY(libffamodules.so) 
R__LOAD_LIBRARY(libg4vertex.so) 
R__LOAD_LIBRARY(libglobalvertex.so) 
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libcentrality.so) 
R__LOAD_LIBRARY(libcalotrigger.so) 

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");

  if (db)
  {
    printf("Server info: %s\n", db->ServerInfo());
  }
  else
  {
    printf("bad\n");
  }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];


  for (int is = 0; is < 64; is++)
  {
    sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
    printf("%s \n" , sql);

    res = db->Query(sql);

    int nrows = res->GetRowCount();

    int nfields = res->GetFieldCount();
    for (int i = 0; i < nrows; i++) {
      row = res->Next();
      for (int j = 0; j < nfields; j++) {
        scaledowns[is] = stoi(row->GetField(j));
      }
      delete row;
    }

    delete res;
  }
  delete db;
}


void Fun4All_G4_Pi0_Tbt_3(const std::string &inputRootFile, const int nEvents = 1000000000)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(inputRootFile);
  int runnumber = runseg.first;

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  Fun4AllInputManager *dstInput = new Fun4AllDstInputManager("DST");
  dstInput->AddFile(inputRootFile);
  se->registerInputManager(dstInput);

  MinimumBiasClassifier *mb = new MinimumBiasClassifier();
  se->registerSubsystem(mb);

  CentralityReco *cent = new CentralityReco();
  se->registerSubsystem(cent);

  int m_scaledowns[64];
  get_scaledowns(runnumber, m_scaledowns);

  std::string baseFilename = inputRootFile.substr(inputRootFile.find_last_of("/") + 1);
  baseFilename = baseFilename.substr(0, baseFilename.find(".root"));
  std::string outputfile = "output_" + baseFilename + "_gbcemc_eval.root";

  CaloCalibEmc_eta *eval = new CaloCalibEmc_eta("CEMC_CALIB_PI0", outputfile);
  eval->SetScaledowns(m_scaledowns);
  eval->setRunNumber(runnumber);
  se->registerSubsystem(eval);

  se->run(nEvents);
  se->End();

  delete se;
  delete mb;
  delete cent;
  delete eval;
  delete dstInput;

  gSystem->Exit(0);
}

#endif  // MACRO_FUN4ALLG4SLOPECAL_C



