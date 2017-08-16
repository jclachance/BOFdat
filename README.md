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
- Metabolites

### Experimental dry weight measurement:
In patch 1 the default values are those for E.coli in M9 minimal media in exponential phase in optimal growth conditions.

CELL_WEIGHT = 280
DNA_RATIO = 0.05
RNA_RATIO = 0.20
rRNA_RATIO = 0.8
tRNA_RATIO = 0.1
mRNA_RATIO = 0.1
PROTEIN_RATIO = 0.6
LIPID_RATIO = 0.1
R_WEIGHT = 284.486
METAB_RATIO = 0.05

Default values are obtained from:
http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/

### General input files

- Organism GenBank annotation, format = genbank, .gb,.gbff (any format supported by BioPython)
- Organism fasta DNA, format may be .fa .fna, .faa (any format supported by BioPython)
- Experimental dry weight measurement, optional: see default values

### Input files per molecule class:
The package takes a path to the files as input.

- DNA: a path to the DNA fasta

- RNA: path to the genbank annotation, path to the model and a path to the transcriptomic data. Transcriptomic should be a comma-separated file with 2 columns. The first column is the gene identifier and the second is the abundance.

- Proteins: proteomics data, a comma-separated file with 2 columns where the first column is your gene_id and the second column is the abundance obtained from experimental data.

- Lipids: lipidomic data. 2 comma-separated files and a path to the model. The first file contains the lipidomic results, 2 column, one is the raw lipid name from the lipidomic results the second is the abundance. The second file is the conversion of the lipid names to a BiGG identifiers. Note that lipid names should be unique. 

- Metabolites: metabolomic data. 2 comma-separated files and a path to the model. The first file contains the metabolomic results, 2 column, one is the raw lipid name from the metabolomic results the second is the abundance. The second file is the conversion of the metabolite names to a BiGG identifiers. Note that metabolite names should be unique. 

### Note

DNA, RNA and protein components are ubiquitous in life and their components should be present in all biomass function in BiGG identifiers.
Lipids and metabolites types differ from a specie to another and should first be extracted qualitatively from metabolomics and lipidomics data. The identified components should be added to the model through manual curation. It is the responsibility of the modeller to ensure the quality of his/her model. 



