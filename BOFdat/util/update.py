"""
Update
======

This module updates BOFsc in the model.

"""
import pandas as pd
import warnings

#------------------------------#
#Functions with multiple use elsewhere in the package

def _import_csv_file(path):
    # 1- Verify if the number of columns is 3
    # in which case the index column is probably the first one
    # 1- Verify number of columns
    csv_file = pd.read_csv(path)
    if len(csv_file.columns) == 2:
        pass
    elif len(csv_file.columns) == 3:
        warnings.warn('File has 3 columns, assuming index on first column')
        csv_file = pd.read_csv(path, index_col=0)
    else:
        raise ImportError('The file format is inappropriate')
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

def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')


def _import_essentiality(path):
    two_col_df = _import_csv_file(path)
    gene_id = [str(i) for i in two_col_df.iloc[0:, 0]]
    essentiality = [str(i) for i in two_col_df.iloc[0:, 1]]
    conform_df = pd.DataFrame({'Genes': gene_id, 'Measured_growth': essentiality},
                              columns=['Genes', 'Measured_growth'])
    return conform_df

def _get_biomass_objective_function(model):
    from cobra.util.solver import linear_reaction_coefficients
    return linear_reaction_coefficients(model).keys()[0]

#------------------------------#

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
    biomass = _get_biomass_objective_function(model)
    find_in_biomass(biomass, dict_of_coefficients)

def update_biomass_metabolites(dict_of_coefficients, model):
    """
    Updates the biomass coefficients given the input dictionary.

    :param dict_of_coefficients: a dictionary of metabolites and coefficients

    :param model: the model to be modified

    :return: none
    """
    biomass = _get_biomass_objective_function(model)
    find_metabolites_in_biomass(biomass, dict_of_coefficients)

def update_maintenance(gams,model,RNA_atp):
    #Get the biomass objective function
    biomass = _get_biomass_objective_function(model)
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
        old_biomass = _get_biomass_objective_function(model)
        old_biomass.remove_from_model()
        biomass = Reaction('BIOMASS')
        model.add_reaction(biomass)
        biomass.add_metabolites(base_biomass)
        biomass.objective_coefficient = 1.
    else:
        return pd.DataFrame.from_records([(k.id,v) for k,v in base_biomass.iteritems()],columns=['Metabolite','Coefficient'])

def determine_coefficients(list_of_metab, model, weight_fraction):
    RATIO = float(1) / len(list_of_metab)
    dict_of_coefficients = {}

    for m in list_of_metab:
        total_weight = RATIO * weight_fraction
        mol_weight = model.metabolites.get_by_id(m).formula_weight
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell
        dict_of_coefficients[m] = -mmols_per_gDW

    return dict_of_coefficients

def convert_to_dictionary(path_to_biomass):
    """
    Function that allows to convert the output of any step of BOFdat into a dictionary

    :param path_to_biomass: path to a 2 column csv file of metabolite identifiers and coefficients
    :param path_to_model: path to the model in supported format
    :return:
    """
    biomass_df = _import_csv_file(path_to_biomass)

    return {biomass_df.iloc[i, 0]: biomass_df.iloc[i, 1] for i in range(len(biomass_df))}

def save_biomass(biomass,file_path):
    """
    Save the result of a step of BOFdat as a standard 2 column csv file

    :param biomass: a dictionary of metabolite identifiers and coefficients
    :param file_path: file name and path under which to save the biomass
    :return:
    """
    df = pd.DataFrame.from_records([(k,v) for k,v in biomass.iteritems()],
                                     columns=['Metabolites','Coefficients'])
    df.to_csv(file_path)


