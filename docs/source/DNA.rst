DNA
===

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
