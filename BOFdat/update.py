"""
Update
======

This module updates BOFsc in the model.

"""
#def change_gam():
def get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return linear_reaction_coefficients(model).keys()[0]

def remove_biomass_metab(dict_of_metab, biomass):
    biomass.add_metabolites(dict_of_metab)
    # print('The actual model %s solves in %s and its biomass contains %s metabolites'
    #      % (model.name, model.optimize(), len(biomass.metabolites)))

def add_biomass_metab(dict_of_metab, biomass):
    biomass.add_metabolites(dict_of_metab)
    # print('The actual model %s solves in %s and its biomass contains %s metabolites'
    #      % (model.name, model.optimize(), len(biomass.metabolites)))
    # print(biomass.reaction)

def find_in_biomass(biomass, dict_of_coefficients):
    # 1- Find metabolites in BIOMASS
    for k, v in dict_of_coefficients.iteritems():
        metab = k
        dict_addition = {k: -v}
        for r in biomass.reactants:
            if metab.id == r.id:
                #print('Found %s in biomass reaction' % (metab.id,))
                # Get the current coefficient for metab in biomass
                current_coefficient = biomass.get_coefficient(r.id)
                dict_removal = {r: -current_coefficient}
                remove_biomass_metab(dict_removal, biomass)
                add_biomass_metab(dict_addition, biomass)

def find_metabolites_in_biomass(biomass, dict_of_coefficients):
    # 1- Find metabolites in BIOMASS
    for k, v in dict_of_coefficients.iteritems():
        metab = k
        dict_addition = {k: -v}
        for r in biomass.reactants:
            if metab.id[:-2] == r.id[:-2]:
                print('Found %s in biomass reaction' % (metab.id,))
                # Get the current coefficient for metab in biomass
                current_coefficient = biomass.get_coefficient(r.id)
                dict_removal = {r: -current_coefficient}
                remove_biomass_metab(dict_removal, biomass)
                add_biomass_metab(dict_addition, biomass)

def update_biomass(dict_of_coefficients, model):
    """
    Updates the biomass coefficients given the input dictionary.

    :param dict_of_coefficients: a dictionary of metabolites and coefficients

    :param model: the model to be modified

    :return: none
    """
    biomass = get_biomass_objective_function(model)
    find_in_biomass(biomass, dict_of_coefficients)

def update_biomass_metabolites(dict_of_coefficients, model):
    """
    Updates the biomass coefficients given the input dictionary.

    :param dict_of_coefficients: a dictionary of metabolites and coefficients

    :param model: the model to be modified

    :return: none
    """
    biomass = get_biomass_objective_function(model)
    find_metabolites_in_biomass(biomass, dict_of_coefficients)

def update_maintenance(gams,model,RNA_atp):
    #Get the biomass objective function
    biomass = get_biomass_objective_function(model)
    # GAM is updated in the BOF as the new ATP coefficient in the BOF
    dict_of_coefficients = {biomass.metabolites.get_by_id('atp_c'):gams['GAM'] + RNA_atp}
    find_in_biomass(biomass,dict_of_coefficients)
    #NGAM is the new objective value of the ATP maintenance reaction
    model.reactions.ATPM.lower_bound = abs(gams['NGAM'])


def make_new_BOF(model, update_model=False, update_NGAM=True, *args, **kwargs):
    """
    Generates a complete new biomass objective function.

    :param model: Model object
    :param update_model: Boolean. If True the provided model object will be updated directly, removing the previous objective.
    :param update_NGAM: Boolean. If True the ATPM function will have its lower bound set to NGAM value.
    :param args: Dictionaries of metabolite object and stoichiometric coefficients
    :param kwargs: Maintenance generated with maintenance.py
    :return: If update_model is False returns the new biomass objective function as a dictionary of metabolite objects and stoichiometric coefficients
    """
    from cobra import Reaction
    base_biomass = {}
   
    for d in args:
        base_biomass.update(d)
    ppi_coeff, h2o_coeff, atp_coeff = [], [], []
    remove_keys = []
    for k, v in base_biomass.iteritems():
        if k.id == 'ppi_c':
            ppi_coeff.append(v) 
            remove_keys.append(k)
        elif k.id == 'h2o_c':
            h2o_coeff.append(v) 
            remove_keys.append(k)
        elif k.id == 'atp_c':
            atp_coeff.append(v)
            remove_keys.append(k)
        else:
            pass

    # Calculate the total for the maintenance related metabolites
    ppi_coeff = sum(ppi_coeff)
    # H2O is produced in the hydrolysis of ATP for GAM
    # ATP is produced in RNA polymerisation
    try:
        GAM = kwargs.get('maintenance').get('GAM')
        NGAM = kwargs.get('maintenance').get('NGAM')
        h2o_coeff = -GAM + sum(h2o_coeff)
        atp_coeff = -GAM + sum(atp_coeff)
        if update_NGAM:
            #Force flux through ATP maintenance reaction
            model.reactions.ATPM.lower_bound = abs(NGAM)
        else:
            print('Not changing NGAM')
    except:
        raise (KeyError('Maintenance costs not provided'))
    # Make ATP hydrolysis reaction
    atp_hydrolysis = {model.metabolites.atp_c: atp_coeff,
                      model.metabolites.h2o_c: h2o_coeff,
                      model.metabolites.adp_c: GAM,
                      model.metabolites.pi_c: GAM,
                      model.metabolites.h_c: GAM,
                      model.metabolites.ppi_c: ppi_coeff
                      }

    # Remove the keys that were found before
    for k in remove_keys:
        base_biomass.pop(k, None)

    base_biomass.update(atp_hydrolysis)
    if update_model:
        old_biomass = get_biomass_objective_function(model)
        old_biomass.remove_from_model()
        biomass = Reaction('BIOMASS')
        model.add_reaction(biomass)
        biomass.add_metabolites(base_biomass)
        biomass.objective_coefficient = 1.
    else:
        return base_biomass
