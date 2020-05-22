"""
Step3
=====

This module finds the cluster of specie-specific metabolic end goals and calculates their stoichiometric coefficients.

"""
from BOFdat.core import initial_population
from BOFdat.core import metab_end_goals
from BOFdat.core import group_end_goals
from BOFdat.core import update

def generate_initial_population(population_write_path, model_path, base_biomass_path, exp_essentiality,
                                number_of_populations=3,WEIGHT_FRACTION=0.05,**kwargs):
    """
        This function generates a given number of initial populations. An initial population is matrix of metabolites
        and individuals where each individual is a list of 0 and 1 corresponding to the presence or absence of a given
        metabolite.

        :param population_write_path: The path to write the population to. BOFdat will add "pop_N.csv" to the given path,
        where N is the Nth population generated.
        :param pop_size: The number of populations to generate, default is 3.
        :param model_path: The path to the model.
        :param base_biomass_path: The path to the output of step 1 and 2 of BOFdat in 2 column ".csv" file.
        :param exp_essentiality: Experimental essentiality as a 2 columns ".csv" file.
        First column is the gene identifiers, second column is the binary experimental essentiality ("E":0,"NE":1).
        :param WEIGHT_FRACTION: weight fraction of the category represented, between 0 and 1
        :param kwargs: {'metab_index': the list of metabolites used to generate the population}
        :return: write population files to the name provided
    """
    """
    Add robustness element here
    """
    initial_population.make_initial_population(population_write_path, model_path, base_biomass_path, exp_essentiality,
                                               number_of_populations,WEIGHT_FRACTION,kwargs)

def find_metabolites(model_path, init_pop_path, exp_essentiality_path, base_biomass=True,
                    logbook=True, hall_of_fame=True, history=False, processes=None, **kwargs):
    """
        This function is the core of BOFdat Step 3. It runs the genetic algorithm on a given population.
        The algorithm will optimize the biomass composition so that its metabolite composition provides a gene
        essentiality prediction by the model that best matches experimental essentiality. An initial population of
        individuals generated with the "generate_initial_population" function is evolved for a given number of
        generations by systematically applying the genetic operators (mutation, crossovers and selection).
        The individuals for which the matthews correlation coefficient (MCC) best matches the predicted essentiality
        are selected.

        The output of this function are the logbook and the hall of fame. The logbook shows the progression of the
        evolution, showing the maximum, mean and minimum MCC for each generation. The Hall of Fame contains the
        metabolite composition of the 1000 best individuals generated through an entire evolution.

        :param model_path: Path to the model for which the biomass objective function is defined
        :param init_pop_path: Path to the initial population on which to run the algorithm. The population should be generated with the initial_population module
        :param exp_essentiality_path: Path to the experimental essentiality data. Two columns csv file, for each gene a 0 in the essentiality column indicates a non-essential gene, a 1 an essential one.
        :param base_biomass: default=True, if True a list of metabolites and their coefficients will always be added
        to the individuals generated throughout the evolution. A pre-determined dictionary of metabolites and coefficients
        can be added with kwargs.
        :param logbook: default=True, generates a logbook of fitness data over generations to the path in kwargs
        :param hall_of_fame: default=True, generates a Hall Of Fame of the best individuals of all time to the path in kwargs
        :param history: default=False, NOT FUNCTIONNAL AS OF 0.1.7
        :param processes: defaul=None, the number of cores to use
        :param kwargs:

    """
    metab_end_goals.qual_definition(model_path, init_pop_path, exp_essentiality_path, base_biomass,
                    logbook, hall_of_fame, history, processes, **kwargs)

def cluster_metabolites(outpath,model_path,CONNECTIVITY_THRESHOLD=15,HOF_PERCENT_THRESHOLD=0.2,eps=6,show_frequency=True,show_matrix=True,**kwargs):
    """
    This function uses a clustering algorithm based on metabolite network distance to find the metabolic objectives
    of the cell from the output of the genetic algorithm.

    :param outpath: Path to the outputs (Hall of fame) of the genetic algorithm
    :param model_path: Path to the model for which the biomass objective function is defined
    :param CONNECTIVITY_THRESHOLD: The threshold above which to remove metabolites for network distance calculation
    :param show_frequency: boolean, if *True* will display the frequency of each metabolite in the hall of fames
    :param show_matrix: boolean, if *True* will display the seaborn cluster map for the reduced distance matrix
    :return:
    """
    if show_frequency:
        if kwargs.get('frequency_fig_name'):
            frequency_fig_name = kwargs.get('frequency_fig_name')
        else:
            frequency_fig_name = 'frequency_fig.svg'
    else:
        frequency_fig_name = 'frequency_fig.svg'

    if show_matrix:
        if kwargs.get('matrix_fig_name'):
            matrix_fig_name = kwargs.get('matrix_fig_name')
        else:
            matrix_fig_name = 'matrix_fig.svg'
    else:
        matrix_fig_name = 'matrix_fig.svg'
    keywords = {'frequency_fig_name':frequency_fig_name,
                'matrix_fig_name':matrix_fig_name}
    print(keywords)
    return group_end_goals.cluster_metabolites(outpath,
                                        model_path,
                                        CONNECTIVITY_THRESHOLD,
                                        HOF_PERCENT_THRESHOLD,
	                                    eps,
                                        show_frequency,
                                        show_matrix,
                                        frequency_fig_name,matrix_fig_name)


