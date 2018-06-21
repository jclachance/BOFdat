"""
Metabolic end goals
===================

This module takes the outputs of the genetic algorithm and groups them into clusters of metabolic end goals.

"""

import random
from random import shuffle, randint
from BOFdat.util.update import _import_csv_file,_import_model,_import_base_biomass,_import_essentiality
from BOFdat.util.update import _get_biomass_objective_function
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
#from scipy.spatial.distance import hamming
# Recording data and stats
import matplotlib.pyplot as plt
#import networkx
# Parallel
import multiprocessing
# Timeout imports and definitions
from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


"-----------------------------"


def initIndividual(icls, content):
    return icls(content)


"-----------------------------"


def initPopulation(pcls, ind_init, filename):
    contents = pd.read_csv(filename, index_col=0)
    return pcls(ind_init(contents[c]) for c in contents.columns)


"-----------------------------"


def make_ind(index):
    # Generates an individual with 100 metabolites
    ind_dict = {}
    ind_size = randint(0,float(len(index))/3)
    for i in index:
        if sum(ind_dict.values()) < ind_size:
            ind_dict[i] = randint(0, 1)
        else:
            ind_dict[i] = 0
    return ind_dict


"-----------------------------"


def feasible(growth_rate):
    """Feasibility function for the individual. Returns True if feasible False
	otherwise."""
    if growth_rate > 1e-6:
        return True
    return False


"-----------------------------"


def generate_new_individual(universal_index, model, base_biomass):
    solvable = False
    while solvable == False:
        shuffle_index = [m for m in universal_index if m not in [v.id for v in base_biomass.keys()]]
        shuffle(shuffle_index)
        ind = make_ind(shuffle_index)
        # Remove old biomass and add new one
        old_biomass = _get_biomass_objective_function(model)
        old_biomass.remove_from_model()
        biomass = Reaction('BIOMASS')
        model.add_reaction(biomass)

        # Add a constraint to the biomass
        # Meaning that the final model has to produce DNA and RNA
        biomass.add_metabolites(base_biomass)

        for i in shuffle_index:
            # Add metabolites to the temporary biomass
            if ind.get(i) == 1:
                biomass.add_metabolites({i: -0.1})

        # Set new biomass and test individual
        biomass.objective_coefficient = 1.
        solvable = feasible(model.slim_optimize())

        if solvable == True:
            # Match the individual's metabolites with the original index
            # and make an ordered list accordingly
            ind_list = []
            for i in universal_index:
                # Add metabolites to the temporary biomass
                if ind.get(i) == 1:
                    ind_list.append(1)
                else:
                    ind_list.append(0)
            return ind_list
        else:
            pass


"-----------------------------"


def eval_ind(individual, initial_pop, model, base_biomass, exp_ess, distance):
    # Set this as warning
    model.solver = 'gurobi'
    old_biomass = list(linear_reaction_coefficients(model).keys())[0]  # index removed
    old_biomass.remove_from_model()
    # Make a biomass reaction and optimize for it
    biomass = Reaction('BIOMASS')
    model.add_reaction(biomass)
    index = initial_pop.index
    for i in range(len(index)):
        if individual[i] == 1:
            biomass.add_metabolites({initial_pop.index[i]: -0.1})
    biomass.add_metabolites(base_biomass)
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
    '''
    if distance == 'hd':
        dist = hamming(u, v)
    '''
    if distance == 'mcc':
        dist = matthews_corrcoef(u, v)
    else:
        print('Error: Invalid distance metric')

    return dist, sum(individual)


"-----------------------------"


def pebble_map(toolbox_evaluate, pop, initial_pop, model, base_biomass, exp_ess, distance, processes):
    print(processes)
    initial_pop_iter = repeat(initial_pop)
    model_iter = repeat(model)
    base_biomass_iter = repeat(base_biomass)
    exp_ess_iter = repeat(exp_ess)
    distance_iter = repeat(distance)
    with ProcessPool(processes) as pool:
        future = pool.map(toolbox_evaluate, pop, initial_pop_iter, model_iter, base_biomass_iter, exp_ess_iter, distance_iter, timeout=40)
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


"-----------------------------"


def record_stats(pop, gen, multi_level):
    # Fitness data for generation
    mean_fit = np.mean([ind.fitness.values[0] for ind in pop])
    std_fit = np.std([ind.fitness.values[0] for ind in pop])
    max_fit = np.max([ind.fitness.values[0] for ind in pop])
    min_fit = np.min([ind.fitness.values[0] for ind in pop])

    # Number of metabolites data
    sizes = [sum(ind) for ind in pop]

    current_stats = pd.DataFrame([[mean_fit, std_fit, max_fit, min_fit, sizes]],
                                 index=[gen], columns=multi_level)

    return current_stats


"-----------------------------"


# Use history from the toolbox instead
def get_ind_composition(hof, index):
    fitness_values, biomass_composition = [], []
    for ind in hof:
        fitness_values.append(ind.fitness.values[0])
        small_list = []
        # Find all the metabolites in the individual
        for v, m in zip(ind, index):
            if v == 1:
                small_list.append(m)
        biomass_composition.append(small_list)

    return pd.DataFrame({'Fitness': fitness_values, 'Biomass composition': biomass_composition})


"-----------------------------"


def _genetic_algo(model, initial_pop, exp_ess, base_biomass, toolbox, GA_param, distance, processes,
                  outputs):
    index = [i for i in initial_pop.index]
    print('made it to the GA')
    print(outputs)
    # Initial population evaluation
    pop = toolbox.population_guess()
    if 'history_name' in outputs.keys():
        history.update(pop)
    print('I will evaluate the initial population')
    fits = toolbox.map(toolbox.evaluate, pop, initial_pop, model, base_biomass,exp_ess, distance, processes)
    print('Fitnesses that will be attributed to individuals %s' % (fits,))
    for ind, fit in zip(pop, fits):
        ind.fitness.values = fit
    print('Finished evaluating initial population')

    # Record outputs
    if 'hall_of_fame_name' in outputs.keys():
        # Hall of fame
        hof = HallOfFame(maxsize=outputs.get('hall_of_fame_size'))
        hof.update(pop)

    if 'logbook_name' in outputs.keys():
        multi_level = pd.MultiIndex(levels=[['Fitness', 'Size'], ['Mean', 'Std', 'Max', 'Min', 'Sizes']],
                                    labels=[[0, 0, 0, 0, 1, ], [0, 1, 2, 3, 4]])
        data_logbook = []
        #valid_logbook = []

    # The GA itself
    for g in range(GA_param.get('NGEN')):
        print('Starting generation %s' % (g,))
        offspring = toolbox.select(pop, k=len(pop))
        offspring = list(map(toolbox.clone, offspring))
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            # Keep the parents for comparison
            parent1 = list(map(toolbox.clone, [child1]))[0]
            parent2 = list(map(toolbox.clone, [child2]))[0]
            rand_nb = random.random()
            if rand_nb < GA_param.get('CXPB'):
                print('doing a crossover')
                if child1 == child2:
                    print('generating a new individual')
                    l = toolbox.generate_new_individual(initial_pop.index, model, base_biomass)
                    child2 = creator.Individual(l)
                    if child2 == parent2:
                        print('generating new individual didnt work')

                toolbox.mate(child1, child2)
                # Assess the efficiency of crossover
                if child1 == parent1 or child1 == parent2:
                    print('crossover yielded identical individuals')
                    l = toolbox.generate_new_individual(initial_pop.index, model, base_biomass)
                    child2 = creator.Individual(l)
                    toolbox.mate(child1, child2)
                else:
                    print('crossover actually worked')

                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            rand_nb = random.random()
            if rand_nb < GA_param.get('MUTPB'):
                toolbox.mutate(mutant)
                del mutant.fitness.values

        # The individuals that have been crossed or mutated will have an invalid fitness
        # and will need to be r-evaluated (saves on computing time)
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        print('I will evaluate %s invalid individuals, wish me luck!' % (len(invalid_ind),))
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind, initial_pop, model, base_biomass, exp_ess, distance, processes)
        print(fitnesses)
        for ind, fit in zip(invalid_ind, fitnesses):
            print(sum(ind))
            ind.fitness.values = fit

        print("  Evaluated %s individuals" % (len(invalid_ind),))
        pop[:] = offspring
        # Save stats of the evolution
        print('Trying to compile stats')
        if 'logbook_name' in outputs.keys():
            current_stats = record_stats(pop, g, multi_level)
            # Manually update the logbook
            data_logbook.append(current_stats)
        if 'hall_of_fame_name' in outputs.keys():
            hof.update(pop)
        if 'history_name' in outputs.keys():
            history.update(pop)

    # Outputs
    # Saving the logbook
    if 'logbook_name' in outputs.keys():
        final_logbook = pd.concat(data_logbook)
        final_logbook.to_csv(outputs.get('logbook_name'))

    # Get the best individuals from the hall of fame (metabolites)
    if 'hall_of_fame_name' in outputs.keys():
        hof_df = get_ind_composition(hof, index)
        hof_df.to_csv(outputs.get('hall_of_fame_name'))




"-----------------------------"
"""
This function should be use to increase stability in import of initial population..
Initial population generation tool should be provided with package
def _import_initial_pop(init_pop_path):
	import pandas as pd
	import warnings
	initial_pop =pd.read_csv(init_pop_path,index_col=1)
	#1- Verify number of columns
	if len(proteomics.columns) > 2:
		raise Exception("Your file format is not appropriate, more than 2 columns")
	#2- Verify presence of header
	if type(proteomics.loc[0, 0]) == str and type(proteomics.loc[0, 1]) == str:
		proteomics = proteomics.iloc[1:]
	#3- Remove null data
	if proteomics.isnull().values.any():
		proteomics = proteomics.dropna()
	#4- Verify column order (identifiers first, abundance second)
	try:
		try:
			abundances = [float(i) for i in proteomics.iloc[0:, 1]]
			identifiers = [str(i) for i in proteomics.iloc[0:, 0]]
		except:
			abundances = [float(i) for i in proteomics.iloc[0:, 0]]
			identifiers = [str(i) for i in proteomics.iloc[0:, 1]]
	except Exception:
		raise Exception('The abundances cannot be converted to float.')

	conform_df = pd.DataFrame({'identifiers':identifiers,'abundances':abundances},columns=['identifiers','abundances'])
	#5- Verify redundancy in protein identifiers
	if len(set(conform_df['identifiers'])) == len(conform_df['identifiers']):
		pass
	else:
		warnings.warn('Redundancy in dataset identifiers')
	#6- Make sure that protein id are used
	if len(list(set(conform_df['identifiers']).intersection(set(seq_dict.keys())))) > 0:
		pass
	else:
		raise Exception('Identifiers not protein_id')

	keys = [k for k in conform_df.identifiers]
	values = [v for v in conform_df.abundances]
	return dict(zip(keys, values))
"""
"-----------------------------"



def qual_definition(model_path, init_pop_path, exp_essentiality_path, base_biomass=True,
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
    # Required arguments
    model = _import_model(model_path)
    exp_ess = _import_essentiality(exp_essentiality_path)
    """Will have to make the import of the initial population and experimental essentiality data safer"""
    initial_pop = pd.read_csv(init_pop_path, index_col=0)

    #valid_ess = pd.read_csv(valid_essentiality_path, index_col=0)
    # GA parameters
    default_param = {'CXPB': 0.5,
                     'MUTPB': 0.1,
                     'NGEN': 20,
                     'indpb': 0.005,
                     'distance_measure': 'mcc',
                     'fitness_distance': 1.0,
                     'fitness_size': -0.25}

    if 'GA_param' in kwargs:
        default_param.update(kwargs.get('GA_param'))
    else:
        print('Using default GA parameters')
    # Non required inputs
    if processes is None:
        try:
            processes = multiprocessing.cpu_count()
        except NotImplementedError:
            warn("Number of cores could not be detected - assuming 1.")
            processes = 1

    if base_biomass:
        if 'base_biomass_path' in kwargs:
            base_biomass = _import_base_biomass(base_biomass_path)
            base_biomass = dict(zip([model.metabolites.get_by_id(k) for k in base_biomass['Metabolites']],
                                    [v for v in base_biomass['Coefficients']]))
        else:
            # Give a standard biomass reaction
            bb_id = {'ala__L_c': -0.8582035429224959,
                     'arg__L_c': -0.1968950229490204,
                     'asn__L_c': -0.2046276924644766,
                     'asp__L_c': -0.3120212125365338,
                     'atp_c': -75.55223000000001,
                     'ctp_c': -0.16308309312316666,
                     'cys__L_c': -0.04349107533373362,
                     'datp_c': -0.02365734018153735,
                     'dctp_c': -0.0264672933260632,
                     'dgtp_c': -0.023187313048617483,
                     'dttp_c': -0.024331654502082398,
                     'gln__L_c': -0.17647268259672946,
                     'glu__L_c': -0.33695369012459897,
                     'gly_c': -0.880181177186217,
                     'gtp_c': -0.15908917528888614,
                     'his__L_c': -0.08628587541929372,
                     'ile__L_c': -0.3172325107365821,
                     'leu__L_c': -0.450016536630228,
                     'lys__L_c': -0.29929818581613993,
                     'met__L_c': -0.11698875330288233,
                     'phe__L_c': -0.13737116178622735,
                     'pro__L_c': -0.2469302131726492,
                     'ser__L_c': -0.33144864820741754,
                     'thr__L_c': -0.3337027605043413,
                     'trp__L_c': -0.027235586441322963,
                     'tyr__L_c': -0.0965789766332646,
                     'utp_c': -0.15004465345468998,
                     'val__L_c': -0.4875735097180578}

            base_biomass = {}
            for k, v in bb_id.items():
                base_biomass[model.metabolites.get_by_id(k)] = v

    else:
        base_biomass = {}

    # Output files
    outputs = {}
    if logbook:
        if 'logbook_name' in kwargs:
            outputs['logbook_name'] = kwargs.get('logbook_name')
        else:
            outputs['logbook_name'] = 'logbook.csv'
    if hall_of_fame:
        if 'hall_of_fame_name' in kwargs:
            outputs['hall_of_fame_name'] = kwargs.get('hall_of_fame_name')
        else:
            outputs['hall_of_fame_name'] = 'hall_of_fame.csv'
        if 'hall_of_fame_size' in kwargs:
            outputs['hall_of_fame_size'] = kwargs.get('hall_of_fame_size')
        else:
            outputs['hall_of_fame_size'] = 1000
    if history:
        if 'history_name' in kwargs:
            outputs['history_name'] = kwargs.get('history_name')
        else:
            outputs['history_name'] = 'history.csv'
    print('initializing toolbox')
    # Initialize DEAP toolbox
    toolbox = base.Toolbox()
    creator.create("FitnessMulti", base.Fitness, weights=(1.0, -0.25))
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    toolbox.register("individual_guess", initIndividual, creator.Individual)
    toolbox.register("population_guess",
                     initPopulation,
                     list,
                     toolbox.individual_guess,
                     init_pop_path)
    toolbox.register('generate_new_individual', generate_new_individual)

    # Get the constants for the evaluation function
    toolbox.register("evaluate", eval_ind)
    toolbox.register('map', pebble_map)

    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=default_param.get('indpb'))
    toolbox.register("select", tools.selTournament, tournsize=3)
    # Register history and Decorate the variation operators
    history = History()
    toolbox.decorate("mate", history.decorator)
    toolbox.decorate("mutate", history.decorator)

    # Run genetic algorithm
    print(default_param)
    distance = default_param.get('distance_measure')
    _genetic_algo(model, initial_pop, exp_ess, base_biomass, toolbox, default_param, distance, processes,
                  outputs)

