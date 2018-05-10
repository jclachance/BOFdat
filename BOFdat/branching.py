"""
Branching
======

This module performs metabolite branching analysis to find potential co-enzymes worthy of adding to the BOF.

"""

def get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return list(linear_reaction_coefficients(model).keys())[0]


def assess_solvability(metabolite_list, model):
    from cobra import Reaction
    print('Generating list of solvable metabolites')
    solvable_metab = []
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    atp_hydrolysis = ['atp', 'h2o', 'adp', 'pi', 'h', 'ppi']
    for m in metabolite_list:
        biomass = get_biomass_objective_function(model)
        biomass.remove_from_model()
        BIOMASS = Reaction('BIOMASS')
        model.add_reactions([BIOMASS])
        # Exclude mass balanced metabolite atp_c
        if m.id[:-2] not in atp_hydrolysis:
            model.reactions.BIOMASS.add_metabolites({m: -1.})
            model.reactions.BIOMASS.objective_coefficient = 1.
            solution = model.optimize()
            # If the model can produce that metabolite
            if solution.f > 1e-9:
                solvable_metab.append(m)
        else:
            model.reactions.BIOMASS.objective_coefficient = 1.

    return solvable_metab


def connectivity_analysis(model):
    metab, number_of_rxn = [], []
    for m in model.metabolites:
        metab.append(m.id)
        number_of_rxn.append(len(m.reactions))

    branching_df = pd.DataFrame({'Metab': metab, 'Number of metab': number_of_rxn})
    # Define threshold using inner stats about the data
    THRESHOLD = branching_df['Number of metab'].median() + branching_df['Number of metab'].std()
    branching_df = branching_df[branching_df['Number of metab'] > THRESHOLD]
    branching_df.sort_values('Number of metab', inplace=True, ascending=False)
    metabolite_list = [model.metabolites.get_by_id(m) for m in branching_df['Metab']]
    solvable_metab = assess_solvability(metabolite_list, model)

    dummy_coeff = -0.00223

    return dict(zip(solvable_metab, [dummy_coeff for m in solvable_metab]))

