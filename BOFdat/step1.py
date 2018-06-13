"""
Step1
=====

This module adds the macromolecules and calculates their stoichiometric coefficients.

"""

from BOFdat.core import dna
from BOFdat.core import rna
from BOFdat.core import protein
from BOFdat.core import lipid
from BOFdat.core import maintenance
from BOFdat.util import update


def generate_dna_coefficients(path_to_fasta, path_to_model , DNA_WEIGHT_FRACTION=0.031):
    """
    Generates a dictionary of metabolite:coefficients for the 4 DNA bases from the organism's
    DNA fasta file and the weight percentage of DNA in the cell.

    :param path_to_fasta: a path to the DNA fasta file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param DNA_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """
    dna_coefficients = dna.generate_coefficients(path_to_fasta, path_to_model , DNA_WEIGHT_FRACTION=0.031)
    return dna_coefficients

def generate_rna_coefficients(path_to_genbank, path_to_model, path_to_transcriptomic,
                         RNA_WEIGHT_FRACTION=0.205,
                         rRNA_WEIGHT_FRACTION=0.9,
                         tRNA_WEIGHT_FRACTION=0.05,
                         mRNA_WEIGHT_FRACTION=0.05,
                         identifier='locus_tag'):
    """
    Generates a dictionary of metabolite:coefficients for the 4 RNA bases from the organism's
    GenBank annotated file, total RNA weight percentage, transcriptomic. Alternately, ribosomal,
    transfer and messenger RNA relative abundances can be incorporated otherwise the default 80% rRNA, 10% tRNA and
    10% mRNA are used.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param path_to_transcriptomic: a two column pandas dataframe (gene_id, abundance)

    :param RNA_WEIGHT_FRACTION: the weight fraction of RNA in the entire cell

    :param rRNA_WEIGHT_FRACTION: the fraction of rRNA to total

    :param tRNA_WEIGHT_FRACTION: the fraction of tRNA to total

    :param mRNA_WEIGHT_FRACTION: the fraction of mRNA to total

    :param identifier: the type of identifier in the input file, 'locus_tag' or 'geneID'

    :return: a dictionary of metabolites and coefficients
    """
    rna_coefficients = rna.generate_coefficients(path_to_genbank, path_to_model, path_to_transcriptomic,
                         RNA_WEIGHT_FRACTION,
                         rRNA_WEIGHT_FRACTION,
                         tRNA_WEIGHT_FRACTION,
                         mRNA_WEIGHT_FRACTION,
                         identifier)

    return rna_coefficients

def generate_protein_coefficients(path_to_genbank, path_to_model, path_to_proteomic, PROTEIN_WEIGHT_FRACTION=0.55):
    """
    Generates a dictionary of metabolite:coefficients for the 20 amino acids contained in proteins from the organism's
    GenBank annotated file, total Protein weight percentage and proteomic data.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param path_to_proteomic: a two column pandas dataframe (protein_id, abundance)

    :param PROTEIN_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """
    protein_coefficients = protein.generate_coefficients(path_to_genbank, path_to_model, path_to_proteomic, PROTEIN_WEIGHT_FRACTION)
    return protein_coefficients

def generate_lipid_coefficients(path_to_lipidomic,path_to_conversion_file,
                     path_to_model,
                     LIPID_WEIGHT_FRACTION=0.091,
                     R_WEIGHT = 284.486):
    """
    Generates a dictionary of metabolite:coefficients for the lipid content of the cell. Lipids vary from a specie to another.
    The lipidomic data provides the relative abundance of each lipid specie while the to_bigg_dict allows to convert identifiers given in the lipidomic data to BiGG identifiers for which
    the metabolite weight is known and can be added easily to the biomass.

    :param path_to_lipidomic: a dataframe of metabolites identified in the lipidomic experiment

    :param path_to_conversion_file: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param LIPID_RATIO: measured lipid ratio of the cell, otherwise default

    :param R_WEIGHT: weight of a carbon chain, otherwise default. If the weight of the lipid is not known it will be inferred based on the number of R chains and this given weight.

    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.

    """
    lipid_coefficients = lipid.generate_coefficients(path_to_lipidomic,path_to_conversion_file,
                     path_to_model,
                     LIPID_WEIGHT_FRACTION,
                     R_WEIGHT)
    return lipid_coefficients

def generate_maintenance_costs(path_to_data, path_to_model, show_GAM = False):
    """

    Growth-associated maintenance (GAM) is the ATP cost of assembling macromolecules in the organism.
    This function calculates GAM from provided path to experimental data. This data includes growth rates on
    different carbon sources, the associated uptake rate for each carbon source and the secretion rates of metabolic
    wastes. More information on the format in which to provide the experimental data is available on GitHub.

    :param path_to_data: The data file is the outcome of the HPLC growth, uptake and secretion rate experiment.

    :param path_to_model: The path to the model, json or sbml formats supported

    :param show_GAM: bool, will associate colors with carbon sources for easier display later

    :return: a dictionary {GAM:value, NGAM:value}
    """
    maintenance_costs = maintenance.experimental_maintenance(path_to_data, path_to_model,show_GAM)
    return maintenance_costs
