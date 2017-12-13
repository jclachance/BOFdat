Metabolites
===========

Similar to lipids, metabolites are intricate as they are specie-specific. BOFdat allows to identify the relevant metabolites to be added to the biomass objective function. 2 filter functions are implemented in the metabolite sub-package. 

The *filter_for_model_metab* function compares the metabolomic data with the provided model.
The *filter_for_universal_biomass_metab* function compares the list of metabolites to a list of metabolites previously added in the BOF of GEMs. This table was extracted from the supplementary materials of: Xavier JC, Patil KR, Rocha I. Integration of Biomass Formulations of Genome-Scale Metabolic Models with Experimental Data Reveals Universally Essential Cofactors in Prokaryotes. Metab Eng. 2017;39: 200–208.

Before generating coefficients for metabolites and adding those to the BOF, it is strongly advised to use the filter functions.

Once these filters have been applied, the 2 files that the *generate_coefficients* function take as input are the "abundance" file, which is your raw metabolomic file in 2 columns and the "conversion" file which is the conversion of the name of each of the metabolite present in the metabolomic to your identifiers in BiGG format.

