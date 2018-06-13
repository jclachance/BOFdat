"""
Coenzymes
=========

This module performs metabolite branching analysis to find potential coenzymes worthy of adding to the BOF.

"""

def _get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return list(linear_reaction_coefficients(model).keys())[0]


def _assess_solvability(metabolite_list, model):
    from cobra import Reaction
    print('Generating list of solvable metabolites')
    solvable_metab = []
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    atp_hydrolysis = ['atp', 'h2o', 'adp', 'pi', 'h', 'ppi']
    gases = ['o2','co2']
    for m in metabolite_list:
        biomass = get_biomass_objective_function(model)
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
                solvable_metab.append(m)
        else:
            model.reactions.BIOMASS.objective_coefficient = 1.

    return solvable_metab


def _determine_coefficients(list_of_metab, model, WEIGHT_FRACTION=0.05):
    RATIO = float(1) / len(list_of_metab)
    dict_of_coefficients = {}

    for m in list_of_metab:
        total_weight = RATIO * WEIGHT_FRACTION
        mol_weight = model.metabolites.get_by_id(m).formula_weight
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell
        dict_of_coefficients[m] = mmols_per_gDW

    return dict_of_coefficients


def connectivity_analysis(model,**kwargs):
    """

    :param model:
    :param kwargs:
    :return:
    """
    metab, number_of_rxn = [], []
    for m in model.metabolites:
        metab.append(m.id)
        number_of_rxn.append(len(m.reactions))

    branching_df = pd.DataFrame({'Metab': metab, 'Number of metab': number_of_rxn})
    # Define threshold using inner stats about the data
    THRESHOLD = branching_df['Number of metab'].mean() + branching_df['Number of metab'].std()
    branching_df = branching_df[branching_df['Number of metab'] > THRESHOLD]
    branching_df.sort_values('Number of metab', inplace=True, ascending=False)
    metabolite_list = [model.metabolites.get_by_id(m) for m in branching_df['Metab']]
    solvable_metab = assess_solvability(metabolite_list, model)
    #Determine the stoichiometric coefficient
    dict_of_coeff = _determine_coefficients(solvable_metab,model)

    return dict_of_coeff

