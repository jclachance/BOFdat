from biomass import Update

"""
Usage: create a dna_update object and apply the function get_coefficients to generate a dictionary of
metabolites and coefficients based on experimental measurements. Use update_biomass function to update
the biomass objective function with the generated coefficients.

Inherits: Update
"""
# Mandatory imports

# from biomass import Update

# Objects
# none so far

# Methods

def get_coefficients(path_to_fasta, path_to_model, CELL_WEIGHT=280, DNA_RATIO=0.031):
    """
    Generates a dictionary of coefficients from the fasta file, experimental data
    and user inputted cell weight and DNA ratio (otherwise default value)
    ========================
    Parameters
    :param path_to_fasta: a path to the DNA fasta file of the organism,
    format should be compatible with BioPython SeqIO
    :param path_to_model: a path to the model, format supported are json and xml
    :param CELL_WEIGHT: experimentally measured cell weight in femtograms, float
    :param DNA_RATIO: the ratio of DNA in the entire cell
    ========================
    Return
    :return: a dictionary of metabolites and coefficients
    """
    BASES = ['A', 'T', 'C', 'G']

    # Imports the genome as a fasta file
    # PATH provided by user
    def import_genome(fasta):
        from Bio import SeqIO
        genome = SeqIO.read(fasta, 'fasta')
        return genome

    # Imports the model
    # PATH provided by user
    # Checks for format
    def import_model(path_to_model):
        # Determine format to import
        model_format = path_to_model.split('.')[-1]
        import cobra
        # Import model 2 formats possible so far json and xml
        if model_format == 'json':
            model = cobra.io.load_json_model(path_to_model)
        elif model_format == 'xml':
            model = cobra.io.read_sbml_model(path_to_model)

        return model

    def get_number_of_bases(genome):
        # Counts the number of each letter in the genome
        base_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for element in genome:
            value = base_genome.get(element)
            new_value = value + 1
            base_genome[element] = new_value

        return base_genome

    def get_ratio(base_genome, genome):
        # Get the ratios for each letter in the genome
        ratio_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

        for letter in BASES:
            number_of_base = float(base_genome.get(letter))
            total = float(len(genome))
            ratio = number_of_base / total
            ratio_genome[letter] = ratio

        return ratio_genome

    def convert_to_coefficient(model, ratio_genome, CELL_WEIGHT, DNA_RATIO):
        # Transform the ratios into mmol/gDW
        DNA_WEIGHT = CELL_WEIGHT * DNA_RATIO

        base_to_bigg = {'A': model.metabolites.datp_c, 'T': model.metabolites.dttp_c,
                        'C': model.metabolites.dctp_c, 'G': model.metabolites.dgtp_c}
        coefficients,metabolites = [],[]

        # Calculate the biomass coefficient for each metabolite
        for letter in BASES:
            ratio = ratio_genome.get(letter)
            total_weight = ratio * DNA_WEIGHT
            metab = base_to_bigg.get(letter)
            mol_weight = metab.formula_weight
            mmols_per_cell = (total_weight / mol_weight) * 1000
            mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
            coefficients.append(mmols_per_gDW)
            metabolites.append(base_to_bigg.get(letter))

        DNA_coefficients = dict(zip(metabolites,coefficients))
        return DNA_coefficients

    # Get the biomass coefficients
    genome = import_genome(path_to_fasta)
    base_in_genome = get_number_of_bases(genome)
    ratio_in_genome = get_ratio(base_in_genome, genome)
    biomass_coefficients = convert_to_coefficient(import_model(path_to_model), ratio_in_genome, CELL_WEIGHT,
                                                  DNA_RATIO)
    return biomass_coefficients

'''
The option to update the coefficients of the metabolites in the biomass objective function is left to the user
'''
def update_biomass_coefficients(dict_of_coefficients,model):
    """
    Updates the biomass coefficients given the input dictionary.
    ========================
    Parameters
    :param dict_of_coefficients: dictionary of metabolites and coefficients
    :param model: model to update
    ========================
    Return
    :return: none
    """
    Update.update_biomass(dict_of_coefficients,model)