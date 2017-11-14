"""
Metabolite
==========

Usage: create a dna_update object and apply the function get_coefficients to generate a dictionary of
metabolites and coefficients based on experimental measurements. Use update_biomass function to update
the biomass objective function with the generated coefficients.

Inherits: Update

"""

# Operations to screen data for Universal biomass table and model should be included before the get_coefficient method
def filter_for_model_metab(path_to_bigg_dict, path_to_model):
    def import_model(path_to_model):
        import cobra
        extension = path_to_model.split('.')[-1]
        if extension == 'json':
            model = cobra.io.load_json_model(path_to_model)
        elif extension == 'xml':
            model = cobra.io.read_sbml_model(path_to_model)
        else:
            print('Model format type not supported')
        return model

    # Get the model
    model = import_model(path_to_model)
    # Get the metabolites in model
    model_metab_id = [m.id for m in model.metabolites]
    # Get to_bigg_dict
    import pandas as pd
    to_bigg_df = pd.read_csv(path_to_bigg_dict)
    to_bigg_dict = dict(zip([i for i in to_bigg_df[to_bigg_df.columns[0]]],
                            [i for i in to_bigg_df[to_bigg_df.columns[1]]]))
    # Get the metabolites that are in the model
    model_metab = {k: v for k, v in to_bigg_dict.iteritems() if v in model_metab_id}
    # Get the metabolites that are not in the model but present in OMICs data
    non_model_metab = [k for k in to_bigg_dict.keys() if k not in model_metab_id]
    print("These metabolites were not found in the model but were present in your metabolomic data, "
          "consider adding them to your model: %s " % ([metab for metab in non_model_metab]))

    return model_metab


def filter_for_biomass_metab(path_to_bigg_dict):
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
    to_bigg_df = pd.read_csv(path_to_bigg_dict)
    to_bigg_dict = dict(zip([i for i in to_bigg_df[to_bigg_df.columns[0]]],
                            [i for i in to_bigg_df[to_bigg_df.columns[1]]]))
    # Get the metabolites that are in the universal biomass table
    universal_metab = []
    for m in to_bigg_df[to_bigg_df.columns[1]]:
        if type(m) == str:
            if m[:-2] in all_compounds:
                universal_metab.append(m)
    print("These metabolites were found in the metabolomic data and universal table of biomass components, "
          "consider adding them to your the biomass objective function: %s " % (
              [metab for metab in universal_metab]))
    # Generate a dataframe
    metab_name, metab_id = [], []
    for k, v in to_bigg_dict.iteritems():
        if v in universal_metab:
            metab_name.append(k)
            metab_id.append(v)

    return pd.DataFrame({'metab_name': metab_name, 'metab_id': metab_id}, columns=['metab_name', 'metab_id'])


def generate_coefficients_from_experimental_data(path_to_metabolomic, path_to_bigg_dict, path_to_model, METAB_RATIO=0.029,
                                            CELL_WEIGHT=280):
    """
    This functions generates a dictionary of stoichiometric coefficients for the uodate of the biomass
    objective function from experimental data.
    The metabolomic data taken as an input is supposed to be filtered using other functions.
    THIS FUNCTION DOES FILTER FOR METABOLITES IN THE MODEL OR THE UNIVERSAL TABLE.

    Parameters
    ===================
    :param path_to_metabolomic: the path to a csv file of BiGG identifiers for metabolites and relative abundances
    :param path_to_model: the path to the metabolic model of interest
    :param TOTAL_METAB: the soluble fraction of the cell
    :param CELL_WEIGHT: total cell weight
    ===================
    Return
    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.
    """
    # Operation 0.1
    # Get the total lipid weight in the cell
    METAB_WEIGHT = METAB_RATIO * CELL_WEIGHT

    # Operation 0.2
    def make_compliant_metabolomic(path_to_lipidomic):
        import pandas as pd
        df = pd.read_csv(path_to_lipidomic)
        df.columns = ['lipid_name', 'abundance']
        return df

    def make_compliant_bigg(path_to_bigg_dict):
        import pandas as pd
        df = pd.read_csv(path_to_bigg_dict)
        df.columns = ['lipid_name', 'lipid_id']
        keys = [i for i in df.lipid_name]
        values = [i for i in df.lipid_id]
        return dict(zip(keys, values))

    # Operation 0.3
    def import_model(path_to_model):
        import cobra
        extension = path_to_model.split('.')[-1]
        if extension == 'json':
            model = cobra.io.load_json_model(path_to_model)
        elif extension == 'xml':
            model = cobra.io.read_sbml_model(path_to_model)
        return model

    # Operation 1
    def convert_metabolomic_to_bigg(metabolomic, to_bigg_dict):
        """
        This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
        with their relative abundances.
        :param lipidomic: a two column csv file that contains original lipid names and relative abundances
        :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers
        :return: a dictionary containing BiGG identifiers and their relative abundances
        """
        import pandas as pd
        # Generate the dictionary
        keys, values = [], []
        for i, row in metabolomic.iterrows():
            keys.append(to_bigg_dict.get(row.lipid_name))
            values.append(row.abundance)

        return dict(zip(keys, values))

    # Operation 2
    def get_relative_abundance(bigg_abundance):
        # Calculate relative abundances
        total_peak = sum(bigg_abundance.values())
        return {k: v / total_peak for k, v in bigg_abundance.iteritems()}

    # Operation 3
    def get_metabolite_weight(bigg_abundance, model):
        # Function to change the molecular weight of a lipid that entails R chains of undefined weight
        # Import Universal BiGG model
        def import_universal_model():
            import cobra
            import os
            print(os.getcwd())
            return cobra.io.load_json_model('universal_BiGG_zac.json')

        # universal_model = import_universal_model()
        keys, values = [], []
        for k, v in bigg_abundance.iteritems():
            try:
                # Get metabolite in BiGG ID
                model_metab = model.metabolites.get_by_id(k)
                # Get molecular weight of the compound from chemical formula
                values.append(model_metab.formula_weight)
                keys.append(k)
            except:
                print('Metabolite %s not found in model' % (k,))
        return dict(zip(keys, values))

    # Operation 4
    def calculate_coefficient(weight_dict, relative_abundance, METAB_WEIGHT, CELL_WEIGHT,model):
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
    # 0.2- Make data compliant for the rest of the functions
    metabolomic_compliant = make_compliant_metabolomic(path_to_metabolomic)
    bigg_compliant = make_compliant_bigg(path_to_bigg_dict)
    # 0.3- Get model
    model = import_model(path_to_model)
    # 1- Convert names to BiGG
    bigg_abundance = convert_metabolomic_to_bigg(metabolomic_compliant, bigg_compliant)
    # 2- Get the relative abundance of each metabolite
    rel_abundance = get_relative_abundance(bigg_abundance)
    # 3- Get the weight of each lipid specie
    weight_dict = get_metabolite_weight(bigg_abundance, model)
    # 4- Calculate biomass coefficients
    biomass_coefficients = calculate_coefficient(weight_dict, rel_abundance, METAB_WEIGHT, CELL_WEIGHT,model)

    return biomass_coefficients


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



def update_biomass_coefficient(dict_of_coefficients, model):
    """
    Updates the biomass coefficients given the input dictionary.
    ========================
    Parameters
    :param dict_of_coefficients: dictionary of metabolites and coefficients
    :param model: model to update
    ========================
    Return
    :return: none
    """

    Update.update_biomass(dict_of_coefficients, model)

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
