"""
Coenzymes
=========

This module performs metabolite branching analysis to find potential coenzymes worthy of adding to the BOF.

"""
import pandas as pd

def _get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return list(linear_reaction_coefficients(model).keys())[0]


def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')


def _assess_solvability(metabolite_list, model):
    from cobra import Reaction
    print('Generating list of solvable metabolites')
    solvable_metab = []
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    atp_hydrolysis = ['atp', 'h2o', 'adp', 'pi', 'h', 'ppi']
    gases = ['o2', 'co2']
    for m in metabolite_list:
        biomass = _get_biomass_objective_function(model)
        biomass.remove_from_model()
        BIOMASS = Reaction('BIOMASS')
        model.add_reactions([BIOMASS])
        # Exclude mass balanced metabolite atp_c
        if m.id[:-2] not in atp_hydrolysis and m.id[:-2] not in gases:
            model.reactions.BIOMASS.add_metabolites({m: -1.})
            model.reactions.BIOMASS.objective_coefficient = 1.
            solution = model.optimize()
            # If the model can produce that metabolite
            if solution.f > 1e-9:
                solvable_metab.append(m.id)
        else:
            model.reactions.BIOMASS.objective_coefficient = 1.

    return solvable_metab


def _find_inorganic_ions(model):
    # List of inorganic ions from literature
    #IONS = ['ca2','cl','cobalt2','cu2','fe2','fe3','h','k','mg2','mn2','mobd','na1','nh4','ni2','so4','zn2']
    import os.path
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.abspath(os.path.join(current_dir,os.pardir))
    file_path = os.path.join(parent_dir,'data/BIOMASS_universal_components.csv')
    universal_components = pd.read_csv(file_path)
    IONS = [m[1:] for m in universal_components.iloc[1:, 2].dropna()]
    inorganic_ions = [m.id for m in model.metabolites if m.id[:-2] in IONS and m.id.endswith('_c')]
    return inorganic_ions


def _find_coenzymes(model):
    # 1- Analyze metabolite connectivity
    metab, number_of_rxn = [], []
    for m in model.metabolites:
        metab.append(m.id)
        number_of_rxn.append(len(m.reactions))
    branching_df = pd.DataFrame({'Metab': metab, 'Number of metab': number_of_rxn})

    # 2- Define threshold using inner stats about the data
    THRESHOLD = branching_df['Number of metab'].mean() + branching_df['Number of metab'].std()
    branching_df = branching_df[branching_df['Number of metab'] > THRESHOLD]
    branching_df.sort_values('Number of metab', inplace=True, ascending=False)
    metabolite_list = [model.metabolites.get_by_id(m) for m in branching_df['Metab']]
    solvable_coenzymes = _assess_solvability(metabolite_list, model)

    return solvable_coenzymes


def find_coenzymes_and_ions(path_to_model):
    # 1- Import model
    model = _import_model(path_to_model)
    # 2- Find coenzymes
    coenzymes = _find_coenzymes(model)
    # 3- Find inorganic ions
    inorganic_ions = _find_inorganic_ions(model)
    # 4- Merge
    coenzymes_and_ions = coenzymes + inorganic_ions

    return coenzymes_and_ions