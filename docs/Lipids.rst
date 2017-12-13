Lipids
======

Lipids are intricate metabolites to model in GEMs. Identifying all lipid species composing the lipidome of an organism can be a daunting task. Modern lipidomics method have allowed to make this task easier by identifying all lipid present in the cell in a single experiment. BOFdat supports the use of lipidomic to generate BOFsc. The filter function also allows to compare lipidomics to the existing model. 

The "abundance" file is your raw lipidomic file in 2 columns. The first column should include the name of each compound as they appear in the lipidomic results, the second column is the abundance of each of these molecules:

The "conversion" file is the conversion of the name of each of the lipid species present in the lipidomic (column one of the "abundance" file) to your identifiers in BiGG format. These identifiers are used in the model.


