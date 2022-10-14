library('ChemmineR')
library('ChemmineDrugs')

setwd('/home/ah310/L1000/R_scripts/')

QUIZC.CIDs <- read.csv('QUIZC_cids_inchi_smiles_allcids.csv')

# compounds.list <- vector(mode = "list", length = length(QUIZC.CIDs$pubchem_cid))
# 
# start <- Sys.time()
# for (i in seq(length(QUIZC.CIDs$pubchem_cid))){
#   print(i)
#   compounds.list[i] <- getIds(as.numeric(QUIZC.CIDs$pubchem_cid[i]))
# }
# Sys.time() - start
# 
# save(compounds.list, file = "QUIZC_cids_inchi_smiles_allcompounds_SDFs.RData")

# load("QUIZC_cids_inchi_smiles_allcompounds_SDFs.RData")
# 
# ap.list <- vector(mode = "list", length = length(QUIZC.CIDs$pubchem_cid))
# apfp.list <- vector(mode = "list", length = length(QUIZC.CIDs$pubchem_cid))
# # fp.list <- vector(mode = "list", length = length(QUIZC.CIDs$pubchem_cid))
# 
# for (i in seq(length(QUIZC.CIDs$pubchem_cid))){
# 
#   # compute atom pair descriptors for the given compound SDF
#   apset <- sdf2ap(compounds.list[[i]])
#     # generate fingerprints from atom pairs of the given compound SDF
#   apfpset <- desc2fp(apset)
# 
#   # generate PubChem fingerprints from PubChems SDFs
#   # fpset <- fp2bit(compounds.list[[i]], type=3)
# 
#   ap.list[i] <- apset
#   apfp.list[i] <- apfpset
#   # fp.list[i] <- fpset
# 
# }
# 
# save(ap.list, file = "QUIZC_cids_inchi_smiles_allcompounds_APs.RData")
# save(apfp.list, file = "QUIZC_cids_inchi_smiles_allcompounds_APFPs.RData")
# # save(fp.list, file = "QUIZC_cids_inchi_smiles_allcompounds_FPs.RData")

# load("QUIZC_cids_inchi_smiles_allcompounds_APs.RData")
load("QUIZC_cids_inchi_smiles_allcompounds_APFPs.RData")

PubChem.CID.combs <- combn(seq(length(QUIZC.CIDs$pubchem_cid)), 2)

start <- Sys.time()
# APsim.df <- data.frame(matrix(ncol = 3, nrow = length(PubChem.CID.combs[1, ])))
APFPsim.df <- data.frame(matrix(ncol = 3, nrow = length(PubChem.CID.combs[1, ])))
# FPsim.df <- data.frame(matrix(ncol = 3, nrow = length(PubChem.CID.combs[1, ])))

# colnames(APsim.df) <- c('Compound A', 'Compound B', 'Tanimoto coefficient')
colnames(APFPsim.df) <- c('Compound A', 'Compound B', 'Tanimoto coefficient')
# colnames(FPsim.df) <- c('Compound A', 'Compound B', 'Tanimoto coefficient')

for (i in seq(length(APFPsim.df$`Compound A`))) {
  # APsim.df[i, ] <- list(QUIZC.CIDs[PubChem.CID.combs[1, i], 'pubchem_cid'],
  #                       QUIZC.CIDs[PubChem.CID.combs[2, i], 'pubchem_cid'],
  #                       cmp.similarity(ap.list[[PubChem.CID.combs[1, i]]], ap.list[[PubChem.CID.combs[2, i]]]))
  APFPsim.df[i, ] <- list(QUIZC.CIDs[PubChem.CID.combs[1, i], 'pubchem_cid'],
                          QUIZC.CIDs[PubChem.CID.combs[2, i], 'pubchem_cid'],
                          fpSim(apfp.list[[PubChem.CID.combs[1, i]]], apfp.list[[PubChem.CID.combs[2, i]]], method="Tanimoto"))
  # FPsim.df[i, ] <- list(QUIZC.CIDs[PubChem.CID.combs[1, i], 'pubchem_cid'],
  #                         QUIZC.CIDs[PubChem.CID.combs[2, i], 'pubchem_cid'],
  #                         fpSim(fp.list[[PubChem.CID.combs[1, i]]], fp.list[[PubChem.CID.combs[2, i]]], method="Tanimoto"))
  
  if (i %% 1000000 == 0) {
    print(i)
    write.csv(APFPsim.df, sprintf('APFPsim_Tanimoto_scores_QUIZC_cids_inchi_smiles_%s.csv', i/1000000), row.names = FALSE)
  }
}
Sys.time() - start
# write.csv(APsim.df, 'APsim_Tanimoto_scores_QUIZC_cids_inchi_smiles.csv', row.names = FALSE)
write.csv(APFPsim.df, 'APFPsim_Tanimoto_scores_QUIZC_cids_inchi_smiles.csv', row.names = FALSE)
# write.csv(FPsim.df, 'FPsim_Tanimoto_scores_QUIZC_cids_inchi_smiles.csv', row.names = FALSE)

