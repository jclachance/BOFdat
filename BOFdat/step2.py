"""
Step2
=====

This module finds the co-enzymes and inorganic ions and calculates their stoichiometric coefficients.

"""
from BOFdat.core import coenzymes_and_ions
from BOFdat.core import update

def find_coenzymes_and_ions(model):
    soluble_pool = coenzymes_and_ions.connectivity_analysis(model)
    return soluble_pool