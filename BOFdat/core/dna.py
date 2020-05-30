"""
DNA
===

This module generates BOFsc for the 4 bases of DNA (dATP, dTTP, dCTP and dGTP)

"""
import warnings

BASES = ['A', 'T', 'C', 'G']

# Methods
def _import_genome(fasta):
    from Bio import SeqIO
    try:
        #Import as a single handle genome    
        genome = list(SeqIO.parse(fasta,'fasta'))
        if len(genome) > 1:
            warnings.warn('%s handles in the genome file.This may indicate that your genome is not completely assembled. \nBOFdat will parse the contigs but the stoichiometric coefficients may not be accurate.'%(len(genome),))
    except:
        raise ImportError('The file provided cannot be imported.')
        
    return genome

def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _get_number_of_bases(genome):
    # Counts the number of each letter in the genome
    base_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for record in genome:
        for element in record:
            value = base_genome.get(element.upper())
            if value == None:
                continue
            else:
                new_value = value + 1
                base_genome[element] = new_value

    return base_genome

def _get_ratio(base_genome):
    # Get the ratios for each letter in the genome
    ratio_genome = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    #DNA is double strand so the number of A = number of T and C = G
    at_number = base_genome.get('A') + base_genome.get('T')
    base_genome['A'] = at_number
    base_genome['T'] = at_number
    gc_number = base_genome.get('C') + base_genome.get('G')
    base_genome['G'] = gc_number
    base_genome['C'] = gc_number

    for letter in BASES:
        number_of_base = float(base_genome.get(letter))
        total = 2*(sum(base_genome.values()))
        ratio = number_of_base / total
        ratio_genome[letter] = ratio

    return ratio_genome

def _convert_to_coefficient(model, ratio_genome, CELL_WEIGHT, DNA_RATIO):
    # Transform the ratios into mmol/gDW
    DNA_WEIGHT = CELL_WEIGHT * DNA_RATIO

    DIPHOSPHATE_WEIGHT = 174.951262

    base_to_bigg = {'A': model.metabolites.datp_c, 'T': model.metabolites.dttp_c,
                    'C': model.metabolites.dctp_c, 'G': model.metabolites.dgtp_c}
    coefficients,metabolites = [],[]

    # Calculate the biomass coefficient for each metabolite
    for letter in BASES:
        ratio = ratio_genome.get(letter)
        total_weight = ratio * DNA_WEIGHT
        metab = base_to_bigg.get(letter)
        mol_weight = metab.formula_weight - DIPHOSPHATE_WEIGHT
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        coefficients.append(mmols_per_gDW)
        metabolites.append(base_to_bigg.get(letter))

    DNA_coefficients = dict(zip(metabolites,[-i for i in coefficients]))
    return DNA_coefficients

def generate_coefficients(path_to_fasta, path_to_model , DNA_WEIGHT_FRACTION=0.031):
    """
    Generates a dictionary of metabolite:coefficients for the 4 DNA bases from the organism's
    DNA fasta file and the weight percentage of DNA in the cell.

    :param path_to_fasta: a path to the DNA fasta file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param DNA_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """
    CELL_WEIGHT = 280
    if DNA_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    #Operations
    genome = _import_genome(path_to_fasta)
    base_in_genome = _get_number_of_bases(genome)
    ratio_in_genome = _get_ratio(base_in_genome)
    model = _import_model(path_to_model)
    biomass_coefficients = _convert_to_coefficient(model, ratio_in_genome, CELL_WEIGHT,
                                                  DNA_WEIGHT_FRACTION)
    #Add Pyrophosphate synthesis as the sum of the coefficients
    ppi_coeff = sum(biomass_coefficients.values())
    ppi_dict = {model.metabolites.get_by_id('ppi_c'):-ppi_coeff}
    biomass_coefficients.update(ppi_dict)

    return biomass_coefficients

'''
The option to update the coefficients of the metabolites in the biomass objective function is left to the user
'''
def update_biomass_coefficients(dict_of_coefficients,model):
    """
    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.
    """

    from BOFdat import update
    update.update_biomass(dict_of_coefficients,model)
