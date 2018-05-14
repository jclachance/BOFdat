"""
Initial population
======

This module generates initial population for the genetic algorithm.

"""
import random
from random import shuffle, randint
# import cobra
import pandas as pd
import numpy as np
from itertools import repeat
from deap import creator, base, tools
from deap.tools import History, HallOfFame
from cobra import Reaction
from cobra.flux_analysis import single_gene_deletion
from cobra.util.solver import linear_reaction_coefficients
from sklearn.metrics import matthews_corrcoef

# Recording data and stats
import matplotlib.pyplot as plt
import networkx
# Parallel
import multiprocessing
# Timeout imports and definitions
from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired

class Individual:
    biomass_name = ''
    biomass = {}
    solvability = True

def _get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return list(linear_reaction_coefficients(model).keys())[0]

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


def _make_metab_ind(m,metab_index):
    # Generates an individual with  metabolites
    ind_dict = {}
    for i in metab_index:
        if i.id == m.id:
            ind_dict[i.id] = 1
        else:
            ind_dict[i.id] = 0
    return ind_dict


def _eval_metab(metab, model, exp_ess):
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

    # Apply hamming distance
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

def _pebble_map(eval_func, iterable, model):
    processes = multiprocessing.cpu_count()
    print(processes)
    model_iter = repeat(model)
    with ProcessPool(processes) as pool:
        future = pool.map(eval_func, iterable, model_iter, timeout=40)
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

def _pebble_eval(eval_func,iterable,model,exp_ess):
    processes = 4
    exp_ess_iter = repeat(exp_ess)
    model_iter = repeat(model)
    with ProcessPool(processes) as pool:
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
    exp_ess = pd.read_csv(exp_essentiality, index_col=0)
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
    print(metab_index)
    # 3- Remove unsolvable metabolites
    print('Going to assess solvability')
    solvability = _pebble_map(_assess_solvability,metab_index, model)
    metab_index = [t[0] for t in solvability if t[1] == True]
    # 4- Find the most relevant metabolites for a maximum gene essentiality prediction
    # Generate a population to test mcc of each metabolite one by one
    metab_id = [m.id for m in metab_index]

    result = _pebble_eval(_eval_metab, metab_id, model, exp_ess)
    result_df = pd.DataFrame({'metab': metab_id, 'mcc': result})
    result_df.sort_values('mcc', ascending=False, inplace=True)
    THRESHOLD = result_df['mcc'].std() + result_df['mcc'].median()

    return [m for m in result_df['metab'][result_df['mcc'] > THRESHOLD]]


def _generate_initial_populations(population_name, metab_index, base_biomass, model):
    # Define the population size as a function of the coverage and individual size
    IND_SIZE = 20  # --> chosen as such so that the final individuals are easy to analyze for modellers
    COVERAGE = 10  # --> empirical but generally suggested in literature
    pop_size = COVERAGE * len(metab_index) / IND_SIZE
    # 5- Generate appropriate number of initial population to obtain a significant result after running the genetic algorithm
    # Save the index before modification madness
    df_index = [m for m in metab_index]

    biomass_list = []
    it = 1
    while len(biomass_list) < pop_size:
        print('Im in loop %s and I have %s valid individuals' % (it, len(biomass_list)))
        ind_list = []
        for n in range(pop_size):
            # Generate an ordered index with any metabolites but those present in the base biomass
            # index = [m for m in metab_list if m.id not in [v.id for v in base_biomass.keys()]]
            # Generate an initial population from any metabolite
            index = [m for m in df_index]

            shuffle(index)
            # Make the individual
            ind_dict = _make_pop_ind(IND_SIZE, index)
            biomass = {}
            biomass.update(base_biomass)

            for i in index:
                # Add metabolites to the temporary biomass
                if ind_dict.get(i) == 1:
                    biomass.update({model.metabolites.get_by_id(i): -0.1})

            individual = Individual()
            individual.biomass_name = 'biomass' + str(it) + '_' + str(n)
            individual.biomass = biomass
            individual.solvability = False
            ind_list.append(individual)

        solutions = _pebble_map(_assess_ind_solvability, ind_list, model)
        for ind in solutions:
            if ind.solvability != False:
                biomass_list.append(ind)

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
    print(df)
    # Write the initial population to file
    df.to_csv(population_name)


def make_initial_population(population_name, model, base_biomass, exp_essentiality,number_of_populations=3):
    """
    This function generates the initial population to run the genetic algorithm on.

    :param population_name: The name and path to write the populations to
    :param pop_size: The number of populations to be generated, default=3
    :param model: Model object
    :param base_biomass: The output of step 1 and 2 of BOFdat
    :param exp_essentiality: Experimental essentiality as a 2 columns csv file the output of the
    :return:
    """
    # Convert base_biomass dataframe to dictionary
    base_biomass = dict(zip([model.metabolites.get_by_id(k) for k in base_biomass['Metabolites']],
                            [v for v in base_biomass['Coefficients']]))
    #1- Make the metabolite index
    metab_index = _generate_metab_index(model, base_biomass,exp_essentiality)
    #2- Make the initial populations
    for n in range(number_of_populations):
        pop_name = population_name + '_' + str(n) + '.csv'
        _generate_initial_populations(pop_name, metab_index, base_biomass)
