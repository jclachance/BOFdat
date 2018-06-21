"""
Step2
=====

This module finds the co-enzymes and inorganic ions and calculates their stoichiometric coefficients.

"""
from BOFdat.core import coenzymes_and_ions
from BOFdat.util import update

def find_coenzymes_and_ions(path_to_model):
    """
    This function finds both coenzymes and inorganic ions in the model.
    The coenzymes are found based on the level of connectivity of the metabolites.
    The inorganic ions are found based on prior knowledge of cell ionic composition.

    :param path_to_model: The path to the model, json or sbml formats supported
    :param WEIGHT_FRACTION: The expected weight fraction of the soluble pool

    :return: Dictionary of metabolites and stoichiometric coefficients
    """
    return coenzymes_and_ions.find_coenzymes_and_ions(path_to_model)
