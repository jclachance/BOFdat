"""
Metabolite
==========

This module generates BOFsc for the metabolite content of the cell.

"""
def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        model = cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        model = cobra.io.read_sbml_model(path_to_model)
    else:
        print('Model format type not supported')
    return model

def _import_metabolomic(path_to_lipidomic):
    import pandas as pd
    df = pd.read_csv(path_to_lipidomic)
    df.columns = ['metab_name', 'abundance']
    return df

def _import_conversion(path_to_conversion_file):
    import pandas as pd
    df = pd.read_csv(path_to_conversion_file)
    df.columns = ['metab_name', 'metab_id']
    keys = [i for i in df.metab_name]
    values = [i for i in df.metab_id]
    return dict(zip(keys, values))


# Operations to screen data for Universal biomass table and model should be included before the get_coefficient method
def filter_for_model_metab(path_to_conversion_file, path_to_model):
    """

    :param path_to_conversion_file: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

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
        del to_bigg_df[to_bigg_dict.columns[0]]
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


def filter_for_universal_biomass_metab(path_to_conversion_file):
    """
    Filters the list of metabolites by comparing it to all metabolites that have been added to metabolic models.
    Source: Joana C. Xavier, Kiran Raosaheb Patil and Isabel Rocha, Integration of Biomass Formulations of
    Genome-Scale Metabolic Models with Experimental Data Reveals Universally Essential Cofactors in Prokaryotes,
    Metabolic Engineering, http://dx.doi.org/10.1016/j.ymben.2016.12.002

    :param path_to_bigg_dict:

    :return:
    """
    def import_universal_table():
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

    # Get the model
    all_compounds = import_universal_table()
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


def generate_coefficients(path_to_metabolomic, path_to_conversion_file, path_to_model, METAB_WEIGHT_FRACTION=0.029,
                                            CELL_WEIGHT=280, remove_DNA_RNA_prot=True):
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
    # Operation 0.1
    # Get the total lipid weight in the cell
    METAB_WEIGHT = METAB_WEIGHT_FRACTION * CELL_WEIGHT

    # Operation 0.2

    # Remove molecules previously calculated
    def remove_molecules(bigg_abundance):
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

    # Operation 1
    def convert_metabolomic_to_bigg(metabolomic, to_bigg_dict,remove_DNA_RNA_prot):
        """

        This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
        with their relative abundances.

        :param lipidomic: a two column csv file that contains original lipid names and relative abundances

        :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers

        :return: a dictionary containing BiGG identifiers and their relative abundances
        """
        # Generate the dictionary
        keys, values = [], []
        for i, row in metabolomic.iterrows():
            if type(to_bigg_dict.get(row.metab_name)) != str:
                pass
            else:
                keys.append(to_bigg_dict.get(row.metab_name))
                values.append(row.abundance)
        unfiltered_dict = dict(zip(keys, values))
        if remove_DNA_RNA_prot == False:
            return unfiltered_dict
        elif remove_DNA_RNA_prot == True:
            return remove_molecules(unfiltered_dict)

    # Operation 2
    def get_relative_abundance(bigg_abundance,model):
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

    # Operation 3
    def get_metabolite_weight(bigg_abundance, model):
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
                        print('Metabolite %s not found in model' % (k,))
            values.append(formula_weight)
            keys.append(metab_id)

        return dict(zip(keys, values))

    # Operation 4
    def calculate_coefficient(weight_dict, relative_abundance, METAB_WEIGHT, CELL_WEIGHT, model):
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

    # Operations
    # 0.1- Calculate the total metabolite weight
    # 0.2- Import data
    metabolomic_compliant = _import_metabolomic(path_to_metabolomic)
    bigg_compliant = _import_conversion(path_to_conversion_file)
    # 0.3- Get model
    model = _import_model(path_to_model)
    # 1- Convert names to BiGG
    bigg_abundance = convert_metabolomic_to_bigg(metabolomic_compliant, bigg_compliant,remove_DNA_RNA_prot)
    # 2- Get the relative abundance of each metabolite
    rel_abundance = get_relative_abundance(bigg_abundance,model)
    # 3- Get the weight of each lipid specie
    weight_dict = get_metabolite_weight(bigg_abundance, model)
    # 4- Calculate biomass coefficients
    biomass_coefficients = calculate_coefficient(weight_dict, rel_abundance, METAB_WEIGHT, CELL_WEIGHT, model)

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
