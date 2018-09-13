"""
Step3
=====

This module finds the cluster of specie-specific metabolic end goals and calculates their stoichiometric coefficients.

"""
from BOFdat.core import initial_population
from BOFdat.core import metab_end_goals
from BOFdat.core import group_end_goals
from BOFdat.core import update

def generate_initial_population(population_name, model, base_biomass, exp_essentiality,
                                number_of_populations=3,WEIGHT_FRACTION=0.05,**kwargs):
    """
        This function generates the initial population to run the genetic algorithm on.

        :param population_name: The name and path to write the populations to
        :param pop_size: The number of populations to be generated, default=3
        :param model: Model object
        :param base_biomass: The output of step 1 and 2 of BOFdat
        :param exp_essentiality: Experimental essentiality as a 2 columns csv file the output of the
        :param WEIGHT_FRACTION: weight fraction of the category represented, between 0 and 1
        :param kwargs: {'metab_index': the list of metabolites used to generate the population}
        :return: write population files to the name provided
    """
    """
    Add robustness element here
    """
    initial_population.make_initial_population(population_name, model, base_biomass, exp_essentiality,
                                               number_of_populations,WEIGHT_FRACTION,kwargs)

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
    metab_end_goals.qual_definition(model_path, init_pop_path, exp_essentiality_path, base_biomass,
                    logbook, hall_of_fame, history, processes, **kwargs)

def cluster_metabolites(outpath,model_path,CONNECTIVITY_THRESHOLD=15,HOF_PERCENT_THRESHOLD=0.2,eps=6,show_frequency=True,show_matrix=True,**kwargs):
    """
    This function analyzes the outputs of the genetic algorithm

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


