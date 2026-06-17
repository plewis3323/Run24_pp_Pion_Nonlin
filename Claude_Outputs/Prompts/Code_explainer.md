

## Pre-Context: 

- You are an AI agent physics helper and software developer that is specialized in HEP physics (LHC and RHIC. etc) 
  Your job it read code help not neccarily teach it to me line by line but in terms of the entire system and the flow of 
  how it's built. And put the output in either latex pdf, txt, or other readable format.  



## Main Task: 

- Your main task is to read the code framework for smearing/nonlinearity matching for MC to data comparison for specifically 
  pions inside directory path: /sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/. 
  The relevant code lives in these locations (read these specifically):
    - /sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Blairs_Codes_files/calo_emc_pi0_tbt/
      (CaloCalibEmc_Pi0.cc / CaloCalibEmc_Pi0.h, pi0EtaByEta.cc / pi0EtaByEta.h, and macro/Fun4All_G4_Pi0_Tbt.C)
    - /sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/ (macros and src)
    - Supporting context: /sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Notes/
      (SMEARING_CODE_CHAIN_LOG.txt, SMEARING_EQUATION_SHEET.txt, findings_log.txt)
  What I want you to do is to teach the system of what the code does, how it built and in easy explantion of how to understand it.  
  Focus on Smearing  and/or resolution lines of importance. 


## Final Deliverables: 

- Have Final Output be in Path: 
/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Claude_Outputs/Code_explainer
- If a LaTeX engine (pdflatex) is available, produce a PDF. If not, output .tex + .txt and note that the PDF was not compiled.




