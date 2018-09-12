"""
Initial population
======

This module generates initial population for the genetic algorithm.

"""
from BOFdat.util.update import _import_csv_file,_import_base_biomass,_import_model,_import_essentiality
from BOFdat.util.update import _get_biomass_objective_function, determine_coefficients
import warnings
import random
from random import shuffle, randint
import cobra
import pandas as pd
import numpy as np
from itertools import repeat
from deap import creator, base, tools
from deap.tools import History, HallOfFame
from cobra import Reaction
from cobra.flux_analysis import single_gene_deletion
from cobra.util.solver import linear_reaction_coefficients
from sklearn.metrics import matthews_corrcoef

# Parallel
import multiprocessing
# Timeout imports and definitions
from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


class Individual:
    biomass_name = ''
    biomass = {}
    solvability = True

"""
DEPRECATED --> FUNCTIONS MOVED TO BOFdat.util.update

def _get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return list(linear_reaction_coefficients(model).keys())[0]

def _import_model(path_to_model):
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _import_csv_file(path):
    csv_file = pd.read_csv(path)
    # 1- Verify number of columns
    if len(csv_file.columns) > 2:
        raise Exception("Your file format is not appropriate, more than 2 columns")
    # 2- Verify presence of header
    if type(csv_file.iloc[0:0, 0]) == str and type(csv_file.iloc[0:0, 1]) == str:
        csv_file = csv_file.iloc[1:]
    # 3- Remove null data
    if csv_file.isnull().values.any():
        csv_file = csv_file.dropna()

    return csv_file
def _import_base_biomass(path):
    two_col_df = _import_csv_file(path)
    metabolites = [str(i) for i in two_col_df.iloc[0:, 0]]
    coefficients = [float(i) for i in two_col_df.iloc[0:, 1]]
    base_biomass_df = pd.DataFrame({'Metabolites':metabolites,'Coefficients':coefficients},
                                   columns=['Metabolites','Coefficients'])
    return base_biomass_df

DEPRECATED --> function not used anymore

def _make_metab_ind(m,metab_index):
    # Generates an individual with  metabolites
    ind_dict = {}
    for i in metab_index:
        if i.id == m.id:
            ind_dict[i.id] = 1
        else:
            ind_dict[i.id] = 0
    return ind_dict

"""

def _branching_analysis(model):
    metab, number_of_rxn = [], []
    for m in model.metabolites:
        metab.append(m.id)
        number_of_rxn.append(len(m.reactions))
        branching_df = pd.DataFrame({'Metab': metab, 'Number of metab': number_of_rxn})
        # Define threshold using inner stats about the data
        THRESHOLD = branching_df['Number of metab'].median() + branching_df['Number of metab'].std()
        branching_df = branching_df[branching_df['Number of metab'] > THRESHOLD]
        branching_df.sort_values('Number of metab', inplace=True, ascending=False)

        return [m for m in branching_df['Metab']]

def _eval_metab(metab, model, exp_ess):
    """
    This function is used to evaluate the fitness of each metabolite individually
    :param metab:
    :param model:
    :param exp_ess:
    :return:
    """
    # Set this as warning
    model.solver = 'gurobi'
    old_biomass = list(linear_reaction_coefficients(model).keys())[0]  # index removed
    old_biomass.remove_from_model()
    # Make a biomass reaction and optimize for it
    biomass = Reaction('BIOMASS')
    model.add_reaction(biomass)
    biomass.add_metabolites({model.metabolites.get_by_id(metab): -0.1})
    biomass.objective_coefficient = 1.

    # Generate deletion results --> BOTTLENECK FOR SURE
    deletion_results = single_gene_deletion(model, model.genes, processes=1)

    # Filter the results to get a boolean result
    a = [(str(next(iter(i))), 1) for i in deletion_results[deletion_results['growth'] > 1e-3].index]
    b = [(str(next(iter(i))), 0) for i in deletion_results[deletion_results['growth'] <= 1e-3].index]
    c = a + b
    pred_ess = pd.DataFrame(c, columns=['Genes', 'Predicted_growth'])

    compare_df = pd.merge(right=exp_ess, left=pred_ess, on='Genes', how='inner')

    # Apply mcc
    u = np.array([f for f in compare_df.Measured_growth])
    v = np.array([x for x in compare_df.Predicted_growth])

    return matthews_corrcoef(u, v)

def _assess_solvability(m, model):
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    model.solver = 'gurobi'
    biomass = _get_biomass_objective_function(model)
    biomass.remove_from_model()
    BIOMASS = Reaction('BIOMASS')
    model.add_reactions([BIOMASS])
    model.reactions.BIOMASS.add_metabolites({m: -1.})
    model.reactions.BIOMASS.objective_coefficient = 1.
    solution = model.slim_optimize()
    # If the model can produce that metabolite
    if solution > 1e-9:
        return (m,True)
    else:
        return (m,False)

def _make_pop_ind(ind_size, index):
    # Generates an individual with 100 metabolites
    ind_dict = {}
    for i in index:
        if sum(ind_dict.values()) < ind_size:
            ind_dict[i] = randint(0, 1)
        else:
            ind_dict[i] = 0
    return ind_dict


def _assess_ind_solvability(ind, model):
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    model.solver = 'gurobi'
    biomass = _get_biomass_objective_function(model)
    biomass.remove_from_model()
    BIOMASS = Reaction('BIOMASS')
    model.add_reactions([BIOMASS])
    model.reactions.BIOMASS.add_metabolites(ind.biomass)
    model.reactions.BIOMASS.objective_coefficient = 1.
    solution = model.slim_optimize()
    # If the model can produce that metabolite
    if solution > 1e-9:
        ind.solvability = True
    else:
        ind.solvability = False
    return ind

def _parallel_init(eval_func, iterable, metab_index,base_biomass,model,weight_fraction):
    """
    This function runs the evaluation function in parallel with 3 arguments.
    It is used twice: first to get the metabolite that the model can produce,
    second to verify the solvability of the generated individuals (multiple metabolites)

    """
    processes = 4
    metab_index_iter = repeat(metab_index)
    base_biomass_iter = repeat(base_biomass)
    model_iter = repeat(model)
    weight_fraction_iter = repeat(weight_fraction)
    with ProcessPool(max_workers=processes,max_tasks=4) as pool:
        future = pool.map(eval_func, iterable, metab_index_iter,
                          base_biomass_iter, model_iter, weight_fraction_iter, timeout=400)
        iterator = future.result()
        all_results = []
        while True:
            try:
                result = next(iterator)
                all_results.append(result)
            except StopIteration:
                break
            except TimeoutError as error:
                print("function took longer than %d seconds" % error.args[1])
                result = 0, 100
                all_results.append(result)
            except ProcessExpired as error:
                print("%s. Exit code: %d" % (error, error.exitcode))
            except Exception as error:
                print("function raised %s" % error)
                print(error.traceback)  # Python's traceback of remote process

    return all_results

def _parallel_assess(eval_func, iterable, model):
    """
    This function runs the evaluation function in parallel with 3 arguments.
    It is used twice: first to get the metabolite that the model can produce,
    second to verify the solvability of the generated individuals (multiple metabolites)

    """
    processes = 10
    model_iter = repeat(model)
    with ProcessPool(max_workers=processes,max_tasks=4) as pool:
        future = pool.map(eval_func, iterable, model_iter, timeout=40)
        iterator = future.result()
        all_results = []
        while True:
            try:
                result = next(iterator)
                print(result)
                all_results.append(result)
            except StopIteration:
                break
            except TimeoutError as error:
                print("function took longer than %d seconds" % error.args[1])
                result = 0, 100
                all_results.append(result)
            except ProcessExpired as error:
                print("%s. Exit code: %d" % (error, error.exitcode))
            except Exception as error:
                print("function raised %s" % error)
                print(error.traceback)  # Python's traceback of remote process

    return all_results

def _pebble_eval(eval_func,iterable,model,exp_ess):
    """
    This function runs the evaluation function in parallel with 4 parameters.
    It is used once by the population generator function to score the metabolites in order of mcc.

    """
    processes = 4
    exp_ess_iter = repeat(exp_ess)
    model_iter = repeat(model)
    with ProcessPool(max_workers=processes,max_tasks=4) as pool:
        future = pool.map(eval_func, iterable, model_iter, exp_ess_iter,timeout=40)
        iterator = future.result()
        all_results = []
        while True:
            try:
                result = next(iterator)
                all_results.append(result)
            except StopIteration:
                break
            except TimeoutError as error:
                print("function took longer than %d seconds" % error.args[1])
                result = 0, 100
                all_results.append(result)
            except ProcessExpired as error:
                print("%s. Exit code: %d" % (error, error.exitcode))
            except Exception as error:
                print("function raised %s" % error)
                print(error.traceback)  # Python's traceback of remote process

    return all_results

def _generate_metab_index(model, base_biomass,exp_essentiality):
    #exp_ess = pd.read_csv(exp_essentiality, index_col=0)
    metab_index = [m for m in model.metabolites]
    # 1- Remove metabolites present in the base biomass
    base_biomass_metab = [k.id for k in base_biomass.keys()]
    metab_index = [m for m in metab_index if m.id not in base_biomass_metab]

    # 2- Remove highly branched metabolites
    highly_branched_metab = _branching_analysis(model)
    metab_index = [m for m in metab_index if m.id not in highly_branched_metab]

    #3- Remove metabolites from atp hydrolysis reaction
    atp_hydrolysis = ['atp', 'h2o', 'adp', 'pi', 'h', 'ppi']
    metab_index = [m for m in metab_index if m.id not in atp_hydrolysis]

    #4- Remove unsolvable metabolites
    print('Assessing individual metabolite solvability')
    solvability = _parallel_assess(_assess_solvability,metab_index, model)
    metab_index = [t[0] for t in solvability if t[1] == True]

    #5- Find the most relevant metabolites for a maximum gene essentiality prediction
    # Generate a population to test mcc of each metabolite one by one
    # This allows to remove irrelevant metabolites from the selection
    metab_id = [m.id for m in metab_index]
    result = _pebble_eval(_eval_metab, metab_id, model, exp_essentiality)
    result_df = pd.DataFrame({'metab': metab_id, 'mcc': result})
    result_df.sort_values('mcc', ascending=False, inplace=True)
    THRESHOLD = result_df['mcc'].std() + result_df['mcc'].median()

    return [m for m in result_df['metab'][result_df['mcc'] > THRESHOLD]]


def _generate_initial_populations(population_name, metab_index, base_biomass, model, WEIGHT_FRACTION):
    """
    This function generates one initial population
    """
    # Population size = number of individuals in the population
    # Individual size = number of metabolites in an individual
    # Define the population size as a function of the coverage and individual size
    # This is an option where the individual size is fixed and the population size varies
    IND_SIZE = 20  # --> chosen as such so that the final individuals are easy to analyze for modellers
    COVERAGE = 10  # --> empirical but generally suggested in literature
    POP_SIZE = COVERAGE * len(metab_index) / IND_SIZE

    # This is an option where the population size is fixed and the individual size varies
    # POP_SIZE = 100
    # 5- Generate appropriate number of initial population to obtain a significant result after running the genetic algorithm
    # Save the index before modification madness
    df_index = [m for m in metab_index]

    biomass_list = []
    it = 1
    print("Generating %s populations" % (POP_SIZE,))
    while len(biomass_list) < POP_SIZE:
        print('Loop %s, %s valid individuals' % (it, len(biomass_list)))
        for n in range(POP_SIZE):
            # Generate an ordered index with any metabolites but those present in the base biomass
            # index = [m for m in metab_list if m.id not in [v.id for v in base_biomass.keys()]]
            # Generate an initial population from any metabolite
            index = [m for m in df_index]

            shuffle(index)
            # Make the individual
            # ind_size = randint(1,float(len(metab_index)/3))
            ind_dict = _make_pop_ind(IND_SIZE, index)
            biomass = {}
            biomass.update(base_biomass)
            # Find the metabolites
            list_of_metab = []
            for i in index:
                # Add metabolites to the temporary biomass
                if ind_dict.get(i) == 1:
                    list_of_metab.append(i)
            # Determine coefficients based on weight fraction and molecular weight
            new_elements = determine_coefficients(list_of_metab, model, WEIGHT_FRACTION)

            # Add to the biomass equation
            biomass.update({model.metabolites.get_by_id(k): v for k, v in new_elements.iteritems()})

            individual = Individual()
            individual.biomass_name = 'biomass' + str(it) + '_' + str(n)
            individual.biomass = biomass
            individual.solvability = False
            ind = _assess_ind_solvability(individual, model)
            if ind.solvability != False:
                biomass_list.append(individual)

        # Verify if the generated individual produces a solution that is in
        shuffle(df_index)
        it += 1

    # Convert biomass list into binary individual
    df = pd.DataFrame(index=df_index)
    for ind in biomass_list:
        metab_list = []
        for m in df_index:
            if m in [d.id for d in list(ind.biomass.keys())]:
                metab_list.append(1)
            else:
                metab_list.append(0)
        df[ind.biomass_name] = metab_list

    # Write the population
    df.to_csv(population_name)


def make_initial_population(population_name, path_to_model, base_biomass_path,
                            exp_essentiality_path,number_of_populations,WEIGHT_FRACTION,kwargs):
    """
    This function generates the initial population to run the genetic algorithm on.

    :param population_name: The name and path to write the populations to
    :param pop_size: The number of populations to be generated, default=3
    :param model: Model object
    :param base_biomass: The path to a 2 column csv file of metabolite ID and stoichiometric coefficients that represents the output of step1 and 2
    :param exp_essentiality_path: The path to a 2 column csv file of gene identifiers (same as the model identifiers) and the binary essentiality
    :return:
    """
    #Import all necessary elements
    model = _import_model(path_to_model)
    exp_essentiality = _import_essentiality(exp_essentiality_path)
    base_biomass_df = _import_base_biomass(base_biomass_path)
    # Convert base_biomass dataframe to dictionary
    base_biomass = dict(zip([model.metabolites.get_by_id(k) for k in base_biomass_df['Metabolites']],
                            [v for v in base_biomass_df['Coefficients']]))
    #1- Make the metabolite index
    if kwargs.get('metab_index'):
        metab_index = kwargs.get('metab_index')
    else:
        metab_index = _generate_metab_index(model, base_biomass,exp_essentiality)

    # 2- Make the initial populations in parallel
    pop_names = [population_name + '_' + str(n) + '.csv' for n in range(number_of_populations)]
    _parallel_init(_generate_initial_populations, pop_names, metab_index, base_biomass, model,WEIGHT_FRACTION)



