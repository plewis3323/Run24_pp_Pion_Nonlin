# AI4S Open Science Report — p+p π⁰ Cross-Section Measurement at sPHENIX

**Author:** P. Lewis (plewis3323), sPHENIX Collaboration graduate student
**Date of snapshot:** 2026-08-04
**Purpose of this package:** Provide the AI4S open-source research platform with (1) a description of the physics goal — measuring the neutral-pion invariant cross section in p+p collisions at √s = 200 GeV with sPHENIX Run-24 data — (2) the standard methodology for that measurement, (3) an honest inventory of what has already been built and where each analysis stage stands, and (4) the source code and scripts themselves (in `code/`, mirroring the on-disk layout at BNL SDCC under `/sphenix/user/plewis3323` and `/sphenix/tg/tg01/bulk/plewis3323`).

---

## 1. Physics goal

Measure the π⁰ invariant differential cross section in p+p collisions at √s = 200 GeV using the sPHENIX electromagnetic calorimeter (EMCal), via the π⁰ → γγ decay channel (BR = 98.8%):

```
E d³σ/dp³ = (1/𝓛_int) · 1/(2π p_T Δp_T Δy) · N_raw(p_T) / (ε_reco(p_T) · A · ε_trig(p_T) · BR)
```

This p+p baseline serves two purposes: a QCD measurement in its own right (comparison to NLO pQCD and to PHENIX results), and the denominator for nuclear modification factors (R_AA, R_CP) in the author's parallel Au+Au and O+O program. The measurement was preceded by a completed Au+Au Run-24 π⁰ R_CP analysis (Masters thesis work), whose full machinery — invariant-mass extraction, mixed-event background, efficiency embedding, yield fitting — is being adapted to p+p.

## 2. How the measurement is done (methodology)

The standard chain, each step mapped in §3 to what already exists:

1. **Dataset definition.** Run-24 p+p good-run list; `DST_CALO` / `DST_CALOFITTING` / `DST_TRIGGER` data files enumerated in file lists from the sPHENIX file catalog. Triggers: MBD minimum-bias and EMCal photon triggers, with scaledown bookkeeping per run.
2. **Ntuple production (Fun4All).** A `SubsysReco` module runs over cluster DSTs: applies calorimeter calibrations (CDB), tower status masks, position-corrected clusters, z-vertex cut (|z| < 30 cm typically), and writes an `_eventTree` ntuple of cluster kinematics (E, p_T, η, φ, χ², max-tower coordinates) plus trigger/vertex info. Run on HTCondor, one job per DST segment; outputs `hadd`-merged in size-capped batches.
3. **Diphoton invariant-mass spectra.** Pair all cluster pairs passing photon-ID cuts (cluster χ², minimum cluster energies/p_T — asymmetric cuts, e.g. p_T1 > 1.1, p_T2 > 0.6 GeV; energy asymmetry α < 0.6). Fill M_γγ vs p_T (TH3 `pairInvMassPtEta`). Combinatorial background estimated by **mixed-event** technique (pairing clusters across events with similar vertex/multiplicity), normalized in sidebands away from the π⁰ (and η) peaks.
4. **Raw yield extraction.** In each p_T bin, subtract the normalized mixed-event background and fit the π⁰ peak (Gaussian or Crystal Ball + polynomial residual background) over ~0.12–0.16 GeV window; cross-check with sideband counting. Yields per p_T bin with fit-model systematics.
5. **Reconstruction efficiency × acceptance.** Simulate π⁰ (single-particle gun and/or PYTHIA-8 Detroit tune) through the full Geant4 sPHENIX detector, embed into real or simulated background events, reconstruct with the identical cut chain, and truth-match reco diphotons to generated π⁰ (ΔR matching, e.g. ΔR < 1.1 with asymmetry/χ² cuts). ε(p_T) = N_matched-reco / N_truth. Because a flat-p_T gun does not have the physical spectrum, **reweight truth p_T to a Hagedorn/power-law fit** of measured (PHENIX) π⁰ cross sections so bin migration and cut efficiencies are evaluated on a realistic spectrum.
6. **Energy-scale and resolution systematics.** Data/MC EMCal response differences handled by (a) tower-by-tower π⁰-mass calibration (`pi0EtaByEta` iterative fitting), (b) an MC smearing procedure (multiplicative energy smear + additive η/φ position smear tuned so MC π⁰ mass/width match data — following B. Seidlitz's method), and (c) nonlinearity studies. Pileup (double interactions in p+p) treated as an additional efficiency systematic.
7. **Trigger efficiency and luminosity.** ε_trig(p_T) for the photon trigger measured relative to minimum bias; integrated luminosity from sampled MB events × σ_MBD (the MBD cross section from vernier scans), with per-run scaledown corrections.
8. **Corrections and final assembly.** Bin-shift/unfolding correction for the steeply falling spectrum, per-p_T-bin combination of yields/efficiency/luminosity into the invariant cross section, statistical + systematic uncertainty propagation (yield-extraction variation, energy scale, efficiency/reweighting, trigger, luminosity).

## 3. What is already built — inventory and status

### 3.1 CURRENT FRONTIER (July 2026): p+p π⁰/η efficiency in PYTHIA — `code/P+P_NonLin_Run24_Space/`

The most recent and most directly cross-section-relevant work (files edited up to **2026-07-28**):

- **`Analysis_space/pi0eta_eff_pythia/`** — full production package for π⁰ and η reconstruction efficiency in Run-24 p+p conditions. `src/Pi0EtaEfficiency.cc/.h` (91 KB module) does cluster pairing, ΔR truth-matching, and efficiency numerator/denominator with per-meson p_T weights. `macro/Fun4All_G4_pythia_eff.C` is Pass 1 (particle gun → Geant4 → EMCal clustering, embedded into one MB pythia8-Detroit background DST per Condor job); `macro/run_efficiency.C` is Pass 2 (truth matching → ε(p_T), weights on by default). `condor/pipeline.sh` orchestrates a 10,000-job SUBMIT/WAIT/HADD production. Final efficiency histograms are included (`runEff/pi0_eff_out.root`, `Eta_eff_out.root`, `eff_pi0_corrected.root`).
- **`Main_Reweight_Plan/Analysis/`** — the Hagedorn reweighting chain, the newest physics code: `step1_invyield.C` (truth invariant yield with the 1/(2π p_T Δp_T Δy N_evt) Jacobian — this is literally the cross-section formula being exercised on MC truth), `step2_hagedorn_fit.C` (Hagedorn fit A[exp(−a p_T − b p_T²) + p_T/p₀]^−n to digitized PHENIX π⁰/η cross sections, included under `Weighting_Report/PHENIX_Weighting_Data/`), `step3_make_weights.C` (produces `weights/MC_HAGEDORN_reweight_{pi0,eta}.root`, included).
- **`Analysis_space/pi0eta_eff_pythia_SM/`** — fork adding EMCal smearing (`src/EmcalSmear.h`, re-implementing B. Seidlitz's `smearPhot()`; the smearing TF1 file `function_compare_wide.root` is included) with unit tests and docs.
- **`Analysis_space/pi0eta_eff_pythia_SM_Pileup/`** — fork adding a double-interaction pileup model (`src/PileupSim.h`: z_reco ≈ (z1+z2)/2, η shift at the EMCal radius R = 93.5 cm) with unit tests and docs.
- **`Blairs_Codes_files/`** — reference subset of B. Seidlitz's calibration/smearing code (`pi0EtaByEta.cc` with `smearPhot()`, `figMaker/mass_plot.C` which generates the smearing functions, `ttreeReader` analysis macros).
- **`Notes/`, `Claude_Outputs/`** — forensic logs of the collaborator code chain, smearing equation sheets, and LaTeX reports (weighting reports 1–4, cross-section plan `Chunk1_InvYield_CrossSection_Plan.tex`).

**Status:** efficiency pipeline runs end-to-end for π⁰ and η; Hagedorn reweighting implemented and validated on truth; smearing and pileup systematics frameworks written with unit tests. This covers methodology steps 5 and 6.

### 3.2 Data-side machinery proven on Au+Au Run-24 (directly portable to p+p)

- **`code/Copy_Sources/`, `code/Eta_Copy/`, `code/home_area/MainMastersWork/src/calo_emc_pi0_tbt/`** — `CaloCalibEmc_eta.cc/.h`: the core `SubsysReco` that books the `_eventTree` ntuple and the M_γγ vs p_T TH3s, with `Loop()` and `Loop_background_event_mixing()` (methodology steps 2–3). The `Eta_Copy` variant differs by exactly one line: raw (`CLUSTERINFO_`) vs position-corrected (`CLUSTERINFO_POS_COR_`) cluster node — a physics-relevant switch. `JB_Fun4all.txt` (J.B.'s macro) is the best in-package example of `TriggerAnalyzer` + `MinimumBiasClassifier` handling — directly relevant to p+p trigger normalization (step 7).
- **`code/Run24_AuAu_New_2024_p007_DST_Calo/`** — the production Condor campaign template: good-run selection (`prepare_condor_files.sh`, 54 runs), list batching (`grouping_A.sh`), the most complete Fun4All production macro (`Condor_process/Fun4All_G4_Pi0_Tbt_3.C` with TriggerAnalyzer, MinimumBiasClassifier, GlobalVertexMap, Calo_Calib), size-capped merging (`hadd_batches.py`), and a DST node-tree dump (`exmaple_debug_stage/Run24_auau_DST_node_info.txt`) invaluable for finding node names. This is the template for the p+p production pass (steps 1–2).
- **`code/macro/`** — downstream Au+Au yield chain: TH3 → p_T slices (`Pt_extraction.C`), mixed-event background normalization with sidebands 0.25–0.45 / 1.5–3.0 / 0.0–0.09 GeV (`Background_Sub_Code.C`), signal fits (`Sig_extraction.C`), sideband-counting cross-check (`Sideband-Y-Method.C`), R_CP assembly. Numeric yield/fit-parameter text outputs included (steps 3–4, proven on Au+Au).
- **`code/macro2/`** — earlier tree-production stage; contains the only executed p+p processing so far: `Fun4All_G4_Pi0_Tbt.C` ran over p+p simulation DSTs (`DST_Lists/DST_pp.list`, `DST_Data_pp_run24_Calo*.list`, `DST_Data_pp_run24-Trigger.list` — all included).
- **`code/Final_tuples{,2}/`** — event-count bookkeeping (`CountingTree.C`, `number_R24_AuAu_2024_build.txt`, `R24_pho_prob_Nevents.txt`) — normalization records for the Au+Au ntuple productions.

### 3.3 Efficiency machinery, generation 1–2 (Au+Au embedding) — completed

- **`code/home_area/Fall_2025/pi0_Eff/`** — the central efficiency codebase: `src/etaMesonAnalysis24.cc` (80 KB; reco+truth trees, `recoTruthMatch()` with ΔR/asymmetry/χ²/diphoton-p_T cuts and **p_T-dependent weights from a TH1**, centrality binning, mixed-event background weighting) and `src_agent/etaMesonAnalysis25.cc` (next generation). Embedding macros (`Fun4All_G4_eta_embed.C`, `G4Setup_sPHENIX.C`), Condor harnesses, `MC_PHENIX_reweight.root` (the PHENIX-based p_T reweighting histogram — predecessor of the Hagedorn chain), and the efficiency-curve production macros with numeric results.
- **`code/Spring_Pi0_eff/`** — Spring-2026 Au+Au efficiency campaign (10,000-job embedding into sHijing): `Pi0_Eff_Extract.C` (efficiency at fixed p_T points by three methods incl. Crystal-Ball + pol2 fits), occupancy-corrected R_CP (`RCP24_occ_yield.C`), final efficiency ROOT files included.
- **`code/Agents_Eff_Directory/`** — an agent-assisted efficiency iteration (April 2026) with a complete run log (`AGENT_LOG.md`) and final efficiency plots (PDFs included). Efficiency rises ~0.30 (central) → ~0.78 (peripheral) as expected.
- **`code/Efficiencies/`** — the MC input provenance: file lists for the single-π⁰, single-η, and PYTHIA `DST_CALO_CLUSTER` + `G4Hits` productions used.
- **`code/home_area/MainMastersWork/Efficiencies/pi0Efficiency_src/`** — standalone single-particle π⁰ efficiency module (5 η bins, truth matching), plus the PPG11 reference analysis was studied (not redistributed here; see §5).

### 3.4 Calibration work (feeds the energy-scale systematic)

- **`code/home_area/Spring_2026/calo_emc_pi0_tbt/`** — current fork of the official `pi0EtaByEta` tower-by-tower calibration module with a centrality port (`CENTRALITY_PORT_CHANGELOG.md`), real-data O+O macros (`Fun4All_DATA_EMCal.C`, 31,719-job Condor campaign, local-calib file creation from CDB), fit macros, and procedure/results markdown logs.
- **`code/home_area/macro_OO/`** — ten iterations of O+O tower-by-tower π⁰ calibration (`doFitAndCalibUpdate.C` iteration driver; the per-iteration calibration ROOT files are included, ~125 KB each).
- **`code/home_area/CDB_TTree/`** — CDB (Conditions Database) utilities for writing/reading EMCal and HCal calibration `CDBTTree`s.

### 3.5 Final plotting / thesis assembly

- **`code/home_area/RCP_Masters_work/`** — the git-tracked curated snapshot (github.com/plewis3323/RCP_Masters_work) of the Masters R_CP analysis, including `EFF_Developer_INFO.txt` (23 KB developer notes — best prose documentation in the package), final styled R_CP ROOT files with systematics, sPHENIX plot style, and the full fit/yield chain with numeric outputs.
- **`code/home_area/MainMastersWork/macro/RESULTS/`** — thesis figure macros (Ch. 3/4 plots, PPG11 efficiency fit comparisons with covariance-based error propagation, background-model studies).

### 3.6 Roadmap documents

- **`code/Summer_26_Space/Summer26_Research_List.txt`** — the author's stated research plan: build/verify the data framework for p+p and O+O (Au+Au as tester), single-particle p+p reconstruction-efficiency corrections for η/π⁰, O+O η/π⁰ R_AA and R_CP, calibration tasks.
- **`code/Agents_Eff_Directory/Eff_analysis.md`**, **`code/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/{README.md,CHANGELOG.md,docs/}`** — task specs, run guides, and changelogs.

## 4. Gap analysis — what remains for the p+p cross section

Mapping §2 against §3:

| Step | Status |
|---|---|
| 1. p+p dataset definition | **Partial.** Run-24 p+p DST lists exist (`macro2/DST_Lists/`, ~1,891 `dst_calo_run2pp-*` lists in `MainMastersWork/other_DST/run_24_DST/` — 3 samples included); no curated p+p good-run list yet (`macro/filelist/Run24_pp/` is empty). |
| 2. p+p data ntuple production | **Not yet run at scale.** All machinery proven on Au+Au Run-24 (`Run24_AuAu_New_2024_p007_DST_Calo` campaign); needs pointing at p+p DSTs with p+p trigger/vertex handling. |
| 3. M_γγ spectra + mixed-event background | **Built and proven on Au+Au**; mixed-event normalization may need retuning for low-multiplicity p+p. |
| 4. Raw yield extraction | **Built and proven on Au+Au** (Gaussian and CB+pol2 fits, sideband cross-check). |
| 5. Efficiency × acceptance (p+p) | **Essentially complete** — the July 2026 PYTHIA pipeline with Hagedorn reweighting; final ε(p_T) histograms in the package. |
| 6. Energy scale / smearing / pileup systematics | **Frameworks written** (SM and SM_Pileup forks, unit-tested); need to be run as full systematic variations. |
| 7. Trigger efficiency + luminosity | **Not started** beyond scaledown hooks (`SetScaledowns()` in `CaloCalibEmc_eta`) and the J.B. TriggerAnalyzer example. Needs photon-trigger turn-on curves and σ_MBD-based luminosity. |
| 8. Cross-section assembly + systematics | **Formula exercised on MC truth** (`step1_invyield.C`); a written plan exists (`Chunk1_InvYield_CrossSection_Plan.tex`). Data-side assembly not yet possible until steps 1–2, 7 complete. |

**Suggested next actions:** (a) curate the Run-24 p+p good-run list and clone the `Run24_AuAu_New_2024_p007_DST_Calo` Condor campaign for p+p; (b) measure photon-trigger efficiency vs MB; (c) obtain σ_MBD and per-run scaledowns for the luminosity normalization; (d) run the data-side M_γγ chain, then combine with the existing ε(p_T) via the already-written invariant-yield code.

## 5. Package notes, caveats, and reproducibility

- **Layout:** `code/` mirrors the on-disk layout. Top-level entries come from `/sphenix/tg/tg01/bulk/plewis3323/`; `code/home_area/` mirrors `/sphenix/user/plewis3323/`. Only source code, scripts, docs, small results (< ~2 MB ROOT/text), and provenance file lists were copied — no bulk ntuples (the on-disk footprint is ~4 TB; this package is ~100 MB).
- **Absolute paths:** nearly every macro/script hardcodes `/sphenix/user/plewis3323/...` or `/sphenix/tg/tg01/bulk/plewis3323/...` (e.g. `R__LOAD_LIBRARY`, `BASE_DIR` in `run_efficiency.C`, filelist paths). Reproduction elsewhere requires rewriting these or a BASE_DIR override.
- **Environment:** everything assumes the sPHENIX software stack (`/opt/sphenix/core/bin/sphenix_setup.sh`) with a local install prefix `MYINSTALL=/sphenix/user/plewis3323/install`. The C++ modules build with autotools (`autogen.sh`/`configure`/`make install`) against Fun4All/coresoftware. Headers pinned in `code/home_area/install_include/`; shared libraries must be rebuilt from the included sources. `code/home_area/bash_profile.txt` documents the login environment.
- **Third-party code:** `Blairs_Codes_files/` (B. Seidlitz), `JB_*` files (J.B.), and the `pi0EtaByEta`/PPG11-derived code are sPHENIX-internal collaborator/official code included as reference lineage — treat as sPHENIX collaboration material, not the author's original work. The full PPG11 SkimMaker clone and coresoftware clones were deliberately excluded.
- **Known wrinkles found during packaging:** `reco_truth.job` files request `600000MB` memory (typo); `macro_OO/Fun4All_EMCal.job` still has `Initialdir = /sphenix/user/sregmi/...` (inherited from the code's origin); `fit_eta_slice_macro.C` first line is corrupted (`ginclude`); `run_job_treeAna.job` in `macro2` has a username typo in its Error path; `Sijan_directory_copy/` on disk is a byte-identical duplicate of `macro2` (excluded); the nested `MainMastersWork/RCP_Masters_work` (687 GB) is superseded by the top-level git repo (the version included here).
- **Data access:** the underlying DSTs and ntuples live on sPHENIX storage at BNL SDCC and are enumerated by the included `.list` files; they are not redistributable here.
