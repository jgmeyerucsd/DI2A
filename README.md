# **NOTE** this analysis pipeline was made obsolete with the completion of CsoDIAq:
CsoDIAq preprint: https://www.biorxiv.org/content/10.1101/2021.05.12.443833v1
CsoDIAq code: https://github.com/CCranney/CsoDIAq

# Data analysis code and results from the original Direct Infusion Shotgun Proteome Analysis (DISPA) paper
Python and R code for peptide identification and quantification from direct infusion shotgun proteome analysis.

### Publications:
preprint on chemrxiv: https://doi.org/10.26434/chemrxiv.12367376.v1
Final publication in Nature Methods: https://www.nature.com/articles/s41592-020-00999-z


## Organization
The code is split into folders for "Python" and "R". 
The data needed to do targeted and untargeted quantificaiton is under the "data" folder. 
Outputs from R and Python scripts are in the "outputs" folder. 

R is used to:
1. compute peptide-level FDR  
2. match peptides to proteins in the FASTA file using the biostrings package. Only peptides that match a single protein entry are kept. 
3. generate images of the raw, projected, and library spectra. 

Python is used to:
1. fix the MGF spectral library created by msconvert.exe conversion of the ms2 spectral library
2. do the quantification by accessing the raw data using the pyteomics package.

## Python file descriptions
1. **fixMGFlib.ipynb** - This file contains a few lines of code to add lines of text containing the peptide sequence to the MGF file for use with MSPLIT-DIA
2. **DI2A_make_quant_targets_for_R.ipynb** - code to process MSPLIT-DIA results and add the heavy and light fragment ion masses to use for quantification. Output is ready for protein-level FDR filtering with R to make the target lists for building targeted instrument methods (make_targetedMS2_lists_for_instrument_method.R) and to make the quantitative list for python targeted Quant (make_quantlist_withMass_for_python.R).
3. **DI2A_untargeted_peptide_quantification.ipynb** - Performs untargeted peptide quantification
4. **DI2A_targeted_protein_quantification.ipynb** - Performs targeted peptide quantification using the target list from "make_quantlist_withMass_for_python.R" 

## R file descriptions
1. **getSpectra.R** - code to generate figures of raw, projected, and library spectra
2. **pepFDR.R** - function to compute the peptide-level FDR from MSPLIT-DIA output
3. **getPepFDRloop.R** - code to loop through all the MSPLIT-DIA results in a directory and determine the number of peptides identified at a given FDR
4. **make_quantlist_withMass_for_python.R** - code to read the python-processed MSPLIT-DIA results (processed with "DI2A_make_quant_targets_file_for_R.ipynb"), do protein FDR filtration, and generate a file containing the quantitative targets for use in python to do targeted peptide quantification (DI2A_targeted_protein_quantification.ipynb). 
5. **make_targetedMS2_lists_for_instrument_method.R** - code to read the python-processed MSPLIT-DIA results (processed with "DI2A_make_quant_targets_file_for_R.ipynb"), do protein FDR filtration, and generate the targeted MS2 lists per CV value needed to generate instrument data collection method on the orbitrap fusion lumos. 

## Supplemental Tables
1. **Table S1** gives the list of 395 proteins identified from the best scouting experiment conditions from MCF7 cells
2. **Table S2** gives the MSPLIT-DIA search output from the qualitative analysis of A549 cell samples corresponding to 2,248 unique peptides at <1% peptide-level FDR
3. **Table S3** gives the list of 552 peptides targeted to quantify 552 proteins
