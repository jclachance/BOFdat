"""
The standard input for this class is a dictionnart of metabolites and coefficients where:
keys = metabolites
values = coefficients
CobraPy supports easily this kind of format to update reactions.
Since no real standard exists in BiGG to exists to annotate the Biomass Objective Function:
- Either the name of the biomass objective function should be provided by the user or
- Default = the model takes the current objective reaction and updates it with the current metabolites

"""

def update_biomass(dict_of_coefficients, model):
    """
    Updates the biomass coefficients given the input dictionary.
    ========================
    Parameters
    :param dict_of_coefficients: a dictionary of metabolites and coefficients
    :param model: the model to be modified

    ========================
    Return
    :return: none
    """

    def get_biomass_objective_function(model):
        from cobra.util.solver import linear_reaction_coefficients
        return linear_reaction_coefficients(model).keys()[0]

    def remove_biomass_metab(dict_of_metab, biomass):
        biomass.add_metabolites(dict_of_metab)
        #print('The actual model %s solves in %s and its biomass contains %s metabolites'
        #      % (model.name, model.optimize(), len(biomass.metabolites)))

    def add_biomass_metab(dict_of_metab, biomass):
        biomass.add_metabolites(dict_of_metab)
        #print('The actual model %s solves in %s and its biomass contains %s metabolites'
        #      % (model.name, model.optimize(), len(biomass.metabolites)))
        #print(biomass.reaction)

    def find_metabolites_in_biomass(biomass, dict_of_coefficients):
        # 1- Find metabolites in BIOMASS
        reactants_biomass, biomass_coefficients, reactant_removal = [], [], []
        for k, v in dict_of_coefficients.iteritems():
            metab = k
            coefficient = v
            dict_addition = {k: -v}
            for r in biomass.reactants:
                if metab.id == r.id:
                    #print('Found %s in biomass reaction' % (metab.id,))
                    # Get the current coefficient for metab in biomass
                    current_coefficient = biomass.get_coefficient(r.id)
                    dict_removal = {r: -current_coefficient}
                    remove_biomass_metab(dict_removal, biomass)
                    add_biomass_metab(dict_addition, biomass)

    biomass = get_biomass_objective_function(model)
    find_metabolites_in_biomass(biomass, dict_of_coefficients)