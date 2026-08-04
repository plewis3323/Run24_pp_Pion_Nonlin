# AI4S Core Code Subset — p+p pi0 Cross-Section Measurement (sPHENIX)

Prepared 2026-08-04 by plewis3323. This is a **small, email-sized subset** of the
full AI4S Open Science package: the report plus the handful of source files that
carry the measurement. No ROOT files, PDFs, build artifacts, or bulk data.

Read `report/AI4S_Report.md` first — it gives the physics goal, the measurement
methodology as a numbered chain, and the per-step status of what is built.
Section references below point into that report.

## Contents

```
report/AI4S_Report.md          Physics goal, methodology, inventory, gap analysis.

efficiency/                    Reconstruction efficiency x acceptance (method step 5).
  src/Pi0EtaEfficiency.cc/.h     The SubsysReco: cluster pairing, dR truth matching,
                                 efficiency numerator/denominator with per-meson pT
                                 weights.
  macro/Fun4All_G4_pythia_eff.C  Pass 1 — particle gun -> Geant4 -> EMCal clustering,
                                 embedded into one MB pythia8-Detroit background DST
                                 per Condor job.
  macro/run_efficiency.C         Pass 2 — truth matching -> eff(pT), weights on by
                                 default.
  condor/pipeline.sh             Orchestrates the 10,000-job SUBMIT / WAIT / HADD
                                 production.
  condor/run.sh                  Per-job driver.
  condor/CondorRun.sh            Job wrapper for the Pass-1 embedding.
  condor/CondorHadd.sh           Merge stage.
  condor/merge_local_staged.sh   Local staged-output merge helper.

production/                    Ntuple production campaign (method steps 1-2).
  Fun4All_G4_Pi0_Tbt_3.C         The production macro: TriggerAnalyzer,
                                 MinimumBiasClassifier, GlobalVertexMap, Calo_Calib.
  condor/prepare_condor_files.sh Good-run selection and job-file generation.
  condor/grouping_A.sh           File-list batching.
  condor/run_job_treeAna.job     HTCondor submit description.
  condor/run_sh_treeAna.sh       Per-job execution script.
  condor/hadd_batches.py         Size-capped output merging.
  condor/good_runlist.txt        The good-run list used.

emc_eta/                       Ntuple + invariant-mass machinery (method steps 2-3).
  CaloCalibEmc_eta.cc/.h         The core SubsysReco: books the _eventTree ntuple and
                                 the M_gammagamma vs pT TH3s; contains Loop() and
                                 Loop_background_event_mixing() (mixed-event
                                 combinatorial background).
```

## Where these came from

Paths below are relative to the full package's `code/` directory, which mirrors
the on-disk layout at BNL SDCC.

| Here | Full package |
|---|---|
| `efficiency/` | `P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/` (macro `condor/` lives under `macro/`) |
| `production/` | `Run24_AuAu_New_2024_p007_DST_Calo/D_Calo_run2auau/` |
| `emc_eta/` | `home_area/MainMastersWork/src/calo_emc_pi0_tbt/` |

## What is deliberately not here

- `report/MANIFEST.txt` — it lists the full 1815-file package, so it would only
  describe files absent from this subset.
- The Hagedorn reweighting chain (`Main_Reweight_Plan/Analysis/step1-3`), the
  EMCal smearing and pileup systematics forks (`pi0eta_eff_pythia_SM`,
  `..._SM_Pileup`), the downstream Au+Au yield-extraction macros (`code/macro/`),
  and all efficiency/calibration ROOT outputs. These are described in the report
  and live in the full package.
- Build system files (`Makefile.am`, `configure.ac`, `autogen.sh`) for the
  efficiency module — the two source files here will not build standalone
  without them.

Full package: `AI4S_pp_pi0_CrossSection_Package_2026-08-04.zip`
(41 MB, 1815 files, MD5 `4f4d730d6769593fc2c3dc227f3d741e`) at
`/sphenix/tg/tg01/bulk/plewis3323/AI4S_OpenScience_Package/` on SDCC.
