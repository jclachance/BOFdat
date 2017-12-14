"""
Metabolite
==========

This module generates BOFsc for the metabolite content of the cell.

"""
# Operations to screen data for Universal biomass table and model should be included before the get_coefficient method
#Public
def filter_for_model_metab(path_to_conversion_file, path_to_model):
    """
    Finds metabolites identified through manual curation in model.

    :param path_to_conversion_file: a dictionary converting from the name present in the metabolomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param path_to_model: a path to the model, format supported are json and xml

    :return: updated dictionary with metabolites found in the model
    """
    # Get the model
    model = _import_model(path_to_model)
    # Get the metabolites in model
    # Remove the compartment
    model_metab_id = [m.id[:-2] for m in model.metabolites]
    # Get to_bigg_dict
    import pandas as pd
    to_bigg_df = pd.read_csv(path_to_conversion_file)
    if len(to_bigg_df.columns) > 2:
        del to_bigg_df[to_bigg_df.columns[0]]
    to_bigg_dict = dict(zip([i for i in to_bigg_df[to_bigg_df.columns[0]]],
                            [i for i in to_bigg_df[to_bigg_df.columns[1]]]))

    # Get the metabolites that are in the model
    model_metab = {k: v for k, v in to_bigg_dict.iteritems() if v in model_metab_id}
    print(model_metab)
    # Get the metabolites that are not in the model but present in OMICs data
    non_model_metab = [k for k,v in to_bigg_dict.iteritems() if v not in model_metab_id]
    if len(non_model_metab) != 0:
        print("These metabolites were not found in the model but were present in your metabolomic data, "
          "consider adding them to your model: %s " % (', '.join([metab for metab in non_model_metab]),))

    return pd.DataFrame({'metab_name': model_metab.keys(), 'metab_id': model_metab.values()}, columns=['metab_name', 'metab_id'])

#Private
def _import_universal_table():
    import pandas as pd
    Universal_biomass = pd.read_csv('BIOMASS_universal_components.csv', skiprows=1)
    all_compounds = []
    for compound in Universal_biomass['Biomolecular components - 518']:
        compound = compound[1:]
        all_compounds.append(compound)

    for compound in Universal_biomass['tRNA charged - 20']:
        if type(compound) == str:
            compound = compound[1:]
            all_compounds.append(compound)

    for compound in Universal_biomass['Inorganic ions - 16']:
        if type(compound) == str:
            compound = compound[1:]
            all_compounds.append(compound)
    return all_compounds

#Public
def filter_for_universal_biomass_metab(path_to_conversion_file):
    """
    Filters the list of metabolites by comparing it to all metabolites that have been added to metabolic models.
    Source: Joana C. Xavier, Kiran Raosaheb Patil and Isabel Rocha, Integration of Biomass Formulations of
    Genome-Scale Metabolic Models with Experimental Data Reveals Universally Essential Cofactors in Prokaryotes,
    Metabolic Engineering, http://dx.doi.org/10.1016/j.ymben.2016.12.002

    :param path_to_bigg_dict:

    :return:
    """
    # Get the metabolites
    all_compounds = _import_universal_table()
    # Get to_bigg_dict
    import pandas as pd
    to_bigg_df = pd.read_csv(path_to_conversion_file)
    to_bigg_dict = dict(zip([i for i in to_bigg_df[to_bigg_df.columns[0]]],
                            [i for i in to_bigg_df[to_bigg_df.columns[1]]]))
    # Get the metabolites that are in the universal biomass table
    universal_metab = []
    for m in to_bigg_df[to_bigg_df.columns[1]]:
        if type(m) == str:
            if m in all_compounds:
                universal_metab.append(m)
    print("These metabolites were found in the metabolomic data and universal table of biomass components, "
          "consider adding them to the biomass objective function: %s " % (
              [metab for metab in universal_metab]))
    # Generate a dataframe
    metab_name, metab_id = [], []
    for k, v in to_bigg_dict.iteritems():
        if v in universal_metab:
            metab_name.append(k)
            metab_id.append(v)

    return pd.DataFrame({'metab_name': metab_name, 'metab_id': metab_id}, columns=['metab_name', 'metab_id'])

def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        model = cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        model = cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')
    return model


def _import_metabolomic(path_to_metabolomic):
    import pandas as pd
    import warnings
    metabolomic = pd.read_csv(path_to_metabolomic, header=None)
    # 1- Verify number of columns
    if len(metabolomic.columns) > 2:
        raise Exception("Your file format is not appropriate, more than 2 columns")
    # 2- Verify presence of header
    if type(metabolomic.loc[0, 0]) == str and type(metabolomic.loc[0, 1]) == str:
        metabolomic = metabolomic.iloc[1:]
    # 3- Remove null data
    if metabolomic.isnull().values.any():
        metabolomic = metabolomic.dropna()
    # 4- Verify column order (identifiers first, abundance second)
    try:
        try:
            abundances = [float(i) for i in metabolomic.iloc[0:, 1]]
            identifiers = [str(i) for i in metabolomic.iloc[0:, 0]]
        except:
            abundances = [float(i) for i in metabolomic.iloc[0:, 0]]
            identifiers = [str(i) for i in metabolomic.iloc[0:, 1]]
    except Exception:
        raise Exception('The abundances cannot be converted to float.')

    conform_df = pd.DataFrame({'identifiers': identifiers, 'abundances': abundances},
                              columns=['identifiers', 'abundances'])
    # 5- Verify redundancy in protein identifiers
    if len(set(conform_df['identifiers'])) == len(conform_df['identifiers']):
        pass
    else:
        warnings.warn('Redundancy in dataset identifiers')

    return conform_df

def _import_conversion(path_to_conversion_file):
    import pandas as pd
    import warnings
    conversion_file = pd.read_csv(path_to_conversion_file, header=None)
    # 1- Verify number of columns
    if len(conversion_file.columns) > 2:
        raise Exception("Your file format is not appropriate, more than 2 columns")
    # 2- Verify presence of header
    if type(conversion_file.loc[0, 0]) == str and type(conversion_file.loc[0, 1]) == str:
        conversion_file = conversion_file.iloc[1:]
    # 3- Remove null data
    if conversion_file.isnull().values.any():
        conversion_file = conversion_file.dropna()
    # 4- Assume column order correct and generate dataframe accordingly
    metab_name = [str(i) for i in conversion_file.iloc[0:, 1]]
    model_id = [str(i) for i in conversion_file.iloc[0:, 0]]
    conform_df = pd.DataFrame({'metab_name': metab_name, 'model_id': model_id},
                              columns=['metab_name', 'model_id'])
    # 5- Verify redundancy in identifiers
    if len(set(conform_df['metab_name'])) == len(conform_df['metab_name']):
        pass
    else:
        warnings.warn('Redundancy in dataset identifiers')

    return conform_df

def _remove_molecules(bigg_abundance):
    amino_acids = ['ala__L', 'cys__L', 'asp__L', 'glu__L', 'phe__L', 'gly', 'his__L', 'ile__L', 'lys__L',
                   'leu__L', 'met__L', 'asn__L', 'pro__L', 'gln__L', 'arg__L', 'ser__L', 'thr__L',
                   'val__L', 'trp__L',
                   'tyr__L']
    nucleotides = ['datp', 'dttp', 'dctp', 'dgtp', 'atp', 'utp', 'ctp', 'gtp']

    new_dict = {}
    for k,v in bigg_abundance.iteritems():
        if k in amino_acids or k in nucleotides:
            pass
        else:
            new_dict[k] = v

    return new_dict

def _convert_metabolomic_to_bigg(metabolomic, conversion,remove_DNA_RNA_prot):
    """

    This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
    with their relative abundances.

    :param metabolomic: a two column csv file that contains original lipid names and relative abundances

    :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers

    :return: a dictionary containing BiGG identifiers and their relative abundances
    """
    import pandas as pd

    # Generate the dictionary of lipid_id and relative abundances
    df = pd.merge(left=metabolomic, right=conversion, left_on='identifiers', right_on='metab_name',how='outer')
    # 1- Verify if strictly no match (specifically doubles the length
    if len(df.index) == (len(metabolomic['identifiers']) + len(conversion['metab_name'])):
        # Change column order in conversion file
        conversion.columns = ['model_id', 'metab_name']
        df = pd.merge(left=metabolomic, right=conversion, left_on='identifiers', right_on='metab_name',how='outer')
        if len(df.index) == (len(metabolomic['identifiers'])):
            pass
        else:
            raise Exception('No match was found between metabolomic data and conversion file, review your conversion file.')

    # 2- Verify if mismatch
    elif len(df.index) != (len(metabolomic['identifiers'])):
        raise Exception('Incomplete match, assign model identifier to each metab identified in metabolomic dataset')

    # Combine abundances for metabolites
    df1 = pd.concat([df['model_id'], df['abundances']], axis=1)
    grouped = df1.groupby('model_id').agg(lambda x: sum(x))
    unfiltered_dict = dict(zip([i for i in grouped.index], [i for i in grouped.abundances]))
    # Remove DNA, RNA and proteins associated molecules if requested
    if remove_DNA_RNA_prot == False:
        return unfiltered_dict
    elif remove_DNA_RNA_prot == True:
        return _remove_molecules(unfiltered_dict)

def _get_relative_abundance(bigg_abundance,model):
    # Calculate relative abundances
    total_peak = sum(bigg_abundance.values())
    keys,values = [],[]
    #Convert to model identifiers
    for k, v in bigg_abundance.iteritems():
        for m in model.metabolites:
            if m.id.startswith(k):
                metab_id = m.id
        keys.append(metab_id)
        values.append(v/total_peak)

    return dict(zip(keys,values))

def _get_metabolite_weight(bigg_abundance, model):
    import warnings
    # Import Universal BiGG model
    '''
    Deprecated
    def import_universal_model():
        import cobra
        import os
        print(os.getcwd())
        return cobra.io.load_json_model('universal_BiGG_zac.json')

    # universal_model = import_universal_model()
    '''
    keys, values = [], []
    for k, v in bigg_abundance.iteritems():
        #Find k in the model
        for m in model.metabolites:
            if m.id.startswith(k):
                try:
                    # Get metabolite in BiGG ID
                    model_metab = model.metabolites.get_by_id(m.id)
                    # Get molecular weight of the compound from chemical formula
                    formula_weight = model_metab.formula_weight
                    metab_id = m.id
                except:
                    warnings.warn('Metabolite %s not found in model, consider using the filter function' % (k,))
        values.append(formula_weight)
        keys.append(metab_id)

    return dict(zip(keys, values))

def _calculate_coefficient(weight_dict, relative_abundance, METAB_WEIGHT, CELL_WEIGHT, model):
    keys, values = [], []
    for k, v in weight_dict.iteritems():
        # Generate the total weight of the compound in the cell
        total_weight = METAB_WEIGHT * relative_abundance.get(k)
        # Get molarity of the compound in the cell
        mmols_per_cell = (total_weight / v) * 1000
        # Convert to usable units in BIOMASS equation
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        # Save value
        values.append(mmols_per_gDW)
        keys.append(model.metabolites.get_by_id(k))

    return dict(zip(keys, values))

#Public
def generate_coefficients(path_to_metabolomic, path_to_conversion_file, path_to_model, METAB_WEIGHT_FRACTION=0.029,
                                            remove_DNA_RNA_prot=True):
    """
    This functions generates a dictionary of stoichiometric coefficients for the uodate of the biomass
    objective function from experimental data.
    The metabolomic data taken as an input is supposed to be filtered using other functions.
    For filtering see functions filter_for_model_metab and/or filter_for_biomass_metab

    :param path_to_metabolomic: the path to a csv file of BiGG identifiers for metabolites and relative abundances

    :param path_to_model: the path to the metabolic model of interest

    :param TOTAL_METAB: the soluble fraction of the cell

    :param CELL_WEIGHT: total cell weight

    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.
    """
    CELL_WEIGHT = 280
    if METAB_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    # Operation 0.1
    # Get the total lipid weight in the cell
    METAB_WEIGHT = METAB_WEIGHT_FRACTION * CELL_WEIGHT
    # Operations
    # 1- Get model
    model = _import_model(path_to_model)
    # 2- Import data
    metabolomic = _import_metabolomic(path_to_metabolomic)
    conversion = _import_conversion(path_to_conversion_file)
    # 3- Convert names to BiGG
    bigg_abundance = _convert_metabolomic_to_bigg(metabolomic, conversion,remove_DNA_RNA_prot)
    # 4- Get the relative abundance of each metabolite
    rel_abundance = _get_relative_abundance(bigg_abundance,model)
    # 3- Get the weight of each lipid specie
    weight_dict = _get_metabolite_weight(bigg_abundance, model)
    # 4- Calculate biomass coefficients
    biomass_coefficients = _calculate_coefficient(weight_dict, rel_abundance, METAB_WEIGHT, CELL_WEIGHT, model)

    return biomass_coefficients

def update_biomass_coefficients(dict_of_coefficients, model):
    """

    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.

    """
    from BOFdat import update
    update.update_biomass_metabolites(dict_of_coefficients, model)







































#######################
# Section under progress
#######################
'''
def get_estimated_coefficients(UNIVERSAL_TABLE,model, TOTAL_METAB= 0.05):
    """
    This function generates the coefficients for metabolites that are present in the model and
    in the prokaryotic biomass universal table generated in Rocha et.al. 2017.
    ===================
    Sources:
    [1] Rocha paper 2017
    Parameters
    ===================
    :param universal_table:
    :param model:
    ===================
    Return
    :return:
    """
    def search_model(model,UNIVERSAL_TABLE):
        """
        Search the model metabolites in the universal table of metabolites.
        :param model: the model
        :param UNIVERSAL_TABLE: the universal table
        :return:
        """
        biomass_metabs = []
        for m in model.metabolites:
            if m.id in UNIVERSAL_TABLE.values():
                biomass_metabs.append(m)

        return biomass_metabs

    def attribute_coefficients():
    #Operations
    #1- Get the metabolites that are common between the universal table and the model
    biomass_metabs = search_model(model,UNIVERSAL_TABLE)
    #2- Split the soluble pool into the number of metabolites that were found
    # and attribute resulting coefficient for these metabolites
    dict_of_coefficients = attribute_coefficients(biomass_metabs, )

#Deprecated
def convert_to_bigg(metab_names):
    """
    This function queries the BiGG database for metabolite names and attempts
    to convert the names into BiGG identifiers. The output is a dictionary of
    metabolite_names: bigg_identifier and a message indicating the number of
    metabolites that were successfully converted to BiGG.

    Parameters
    ===================
    :param metab_names: a list of metabolite names to be converted to bigg identifiers
    ===================
    Return
    :return: a message indicating the number of metabolites succesfully converted,
    the non-converted metabolites and a dictionary of metabolite names and associated bigg identifiers.
    """
    #Query BiGG database
    name_to_bigg = {}

    return name_to_bigg

'''
