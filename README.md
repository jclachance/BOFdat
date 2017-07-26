# biomass
## Package to generate biomass coefficients from OMICs data.

Author: Jean-Christophe Lachance
Date: 05-16-2017
Version: 0.1

Most bacterial species do not have all of their components measured experimentally. Genome-scale metabolic models rely both on a defined media and on a precise biomass objective function to generate reliable predictions on flux-states and gene essentiality. Generally, the biomass objective function is fetched from another organism, namely E.coli. This package aims to produce an easy way to generate biomass coefficients that reflect experimental reality by incorporating basic dry weight measurments and relative abundances of macromolecules obtained from multiple OMICs datasets. 

### Molecule classes supported by this package:
- DNA
- RNA
- Proteins
- Lipids
- Metabolites --> to be developped

### Experimental dry weight measurement:
- DNA: user_defined = X, default = 0.03
- RNA: user_defined = X, default = 0.2 
- Proteins: user_defined = X, default = 0.55
- Lipids: user_defined = X, default = 0.1
- Metabolites --> to be developped

Default values are obtained from:
http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/

### General input files

- Organism GenBank annotation, format = genbank, .gb,.gbff (any format supported by BioPython)
- Experimental dry weight measurement, optional: see default values

### Input files per molecule class:
- DNA: whole-genome sequence, fasta format
- RNA:

  1- mRNA: messenger abundance, default = 0.05 and transcriptomic data, format = .csv 2 columns: gene_id, abundance
  
  2- rRNA: ribosomal abundance, default = 0.85
  
  3- tRNA: transfer RNA abundance, default = 0.05
  
- Proteins: proteomics data, format = .csv 2 columns: gene_id, abundance
- Lipids: lipidomics data, format = .csv 2 columns: gene_id, abundance
- Metabolites: to be devleopped, from metabolomic data

### Note

DNA, RNA and protein components are ubiquitous in life and their components should be present in all biomass function in BiGG identifiers.
Lipids and metabolites types differ from a specie to another and should first be extracted qualitatively from metabolomics and lipidomics data. 



