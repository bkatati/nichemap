# Title:
# Niche Partitioning Association of Fungal Genera Correlated with Lower _Fusarium_ and Fumonisin-B1 levels in Maize
(Nichemapping of Fusarium antagonists for Peer Review)

Support files to manuscript:

# i)  *phyloseq.rds* 

is a phyloseq object, which contains the internal and external mycobiome Amplicon Sequence Variants (ASVs) for 2018/2019 maize season. The 'rds' file is the output file generated from the raw read output data of Illumina amplicon sequencing of samples (materials and methods). The processing of the raw reads was done through the pipeline Divisive Amplicon Denoising Algorithym version 2 (DADA2) which generated the ASVs. Run code for the DADA2 may be found at the external link https://benjjneb.github.io/dada2/tutorial_1_8.html. Assignment of the taxa to the generated ASVs was also done in DADA2 using the Unite Fungal Database (UNITE_Community, 2019; DOI: 10.15156/BIO/786353). The raw sequencing data of the samples for the external mycobiome is based on the previous study (Katati et. al., 2023, doi: 10.1128/aem.00078-23). The internal mycobiome data is from this study.

"phyloseq.rds" is run in section-A of the R code "Nichemap" to generate the relative abundances of the internal mycobiome. This generates a file of the internal + external mycobiome for season-1 (2018/2019 cropping season). External mycobiome is later excluded from the csv file, then the file transposed to create the internal mycobiome csv file "InternoBiome.csv."

# ii) "*InternoBiome.csv*"

is full internal fungal microbiome ASV relative abundances (%), which has been generated from step (i). Variable 'wet' implies samples collected from a wetter northerly agroecological zone (AEZ) of Zambia, while 'dry' implies samples from a southerly drier AEZ.

# iii) "*Spearman.csv*"

contains relative abundances (%) of genera identified to have strong negative correlations with each other using Spearman Rank Correlation rho. The file is run in the R code "Nichemap." Abbreviations: FB1 = Fumonisin-B1 (mg/kg); Fus = *Fusarium*; FusGib = *Fusarium* + *Gibberella*; Sten = *Stenocarpella*; Sar = *Sarocladium.*
Note: *Fusarium* (*equiseti*) is detected as a distinct anamorph of the former *Gibberella* (*intricans*) based on the UNITE fungal database (UNITE Community, 2019). *Talaromyces* is assessed as a distinct clade from its asexual anamorph *Penicillium*.

# iv) "*Questionnaire_dataframe.xlsx*" 

contains the questionnaire matrix and dataframe for the agronomic factors (Pest_incidence; No_crop_rotation; No_field_burning; Seed maturing) with potential influence on *Fusarium* relative abundance and FB1 (See manuscript's Materials and Methods on generation of the matrix).

# v) "*dataframe_FusB1.csv*"

is dataframe for agronomic factors (Pest_incidence; No_crop_rotation; No_field_burning; Seed maturing) with potential influence on *Fusarium* relative abundance and FB1 and has been derived from step ‘iv’ above. (See manuscript's Materials and Methods on ordination of *Fusarium* / FB1 levels).

# vi) "*S3_ExternoBiome.csv*"

contains part of external mycobiome data based on 2018/2019 cropping season from the previous study (Katati et. al., 2023, Doi: 10.1128/aem.00078-23). This is for simulation of the correlation of *Sarocladium* (sar) with *Fusarium* (fus). S1 refers to southerly (dry) agro-ecological Zone of Zambia under the 2018/2019 cropping season, whereas N1 refers to northerly (wet) agro-ecological zone under the same cropping season. The simulation is done in R code "Nichemap" under section "Supplemental Data S3."
