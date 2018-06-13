"""
Step3
=====

This module finds the cluster of specie-specific metabolic end goals and calculates their stoichiometric coefficients.

"""
from BOFdat.core import initial_population
from BOFdat.core import chose_metab
#from BOFdat.core import metab_end_goals --> write this code
from BOFdat.core import update

def generate_initial_population(population_name, model, base_biomass, exp_essentiality,number_of_populations=3):
    """
        This function generates the initial population to run the genetic algorithm on.

        :param population_name: The name and path to write the populations to
        :param pop_size: The number of populations to be generated, default=3
        :param model: Model object
        :param base_biomass: The output of step 1 and 2 of BOFdat
        :param exp_essentiality: Experimental essentiality as a 2 columns csv file the output of the
        :return: write population files to the name provided
    """
    initial_population.make_initial_population(population_name, model, base_biomass, exp_essentiality, number_of_populations)

def find_metabolites(model_path, init_pop_path, exp_essentiality_path, base_biomass=True,
                    logbook=True, hall_of_fame=True, history=False, processes=None, **kwargs):
    """
        This function treats the inputs for the genetic algorithm.

        :param model_path: Path to the model for which the biomass objective function is defined
        :param init_pop_path: Path to the initial population on which to run the algorithm. The population should be generated with the initial_population module
        :param exp_essentiality_path: Path to the experimental essentiality data. Two columns csv file, for each gene a 0 in the essentiality column indicates a non-essential gene, a 1 an essential one.
        :param base_biomass: default=True,
        :param logbook: default=True, generates a logbook of fitness data over generations to the path in kwargs
        :param hall_of_fame: default=True, generates a Hall Of Fame of the best individuals of all time to the path in kwargs
        :param history: default=False, NOT FUNCTIONNAL AS OF 0.3
        :param processes: defaul=None, the number of cores to use
        :param kwargs:
        :return:
    """
    chose_metab.qual_definition(model_path, init_pop_path, exp_essentiality_path, base_biomass,
                    logbook, hall_of_fame, history, processes, **kwargs)

#def cluster_metabolites():