"""
Lipid
=====

This module generates BOFsc for the lipid content of the cell.

"""
#Private non-API functions
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

def _import_lipidomic(path_to_lipidomic):
        import pandas as pd
        import warnings
        lipidomic  = pd.read_csv(path_to_lipidomic,header=None)
        # 1- Verify number of columns
        if len(lipidomic.columns) > 2:
            raise Exception("Your file format is not appropriate, more than 2 columns")
        # 2- Verify presence of header
        if type(lipidomic.loc[0, 0]) == str and type(lipidomic.loc[0, 1]) == str:
            lipidomic = lipidomic.iloc[1:]
        # 3- Remove null data
        if lipidomic.isnull().values.any():
            lipidomic = lipidomic.dropna()
        # 4- Verify column order (identifiers first, abundance second)
        try:
            try:
                abundances = [float(i) for i in lipidomic.iloc[0:, 1]]
                identifiers = [str(i) for i in lipidomic.iloc[0:, 0]]
            except:
                abundances = [float(i) for i in lipidomic.iloc[0:, 0]]
                identifiers = [str(i) for i in lipidomic.iloc[0:, 1]]
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
        conversion_file = pd.read_csv(path_to_conversion_file,header=None)
        # 1- Verify number of columns
        if len(conversion_file.columns) > 2:
            raise Exception("Your file format is not appropriate, more than 2 columns")
        # 2- Verify presence of header
        if type(conversion_file.loc[0, 0]) == str and type(conversion_file.loc[0, 1]) == str:
            conversion_file = conversion_file.iloc[1:]
        # 3- Remove null data
        if conversion_file.isnull().values.any():
            conversion_file = conversion_file.dropna()
        #4- Assume column order correct and generate dataframe accordingly
        lipid_name = [str(i) for i in conversion_file.iloc[0:, 1]]
        model_id = [str(i) for i in conversion_file.iloc[0:, 0]]
        conform_df = pd.DataFrame({'lipid_name': lipid_name, 'model_id': model_id},
                                  columns=['lipid_name', 'model_id'])
        # 5- Verify redundancy in identifiers
        if len(set(conform_df['lipid_name'])) == len(conform_df['lipid_name']):
            pass
        else:
            warnings.warn('Redundancy in dataset identifiers')

        return conform_df

def _convert_lipidomics_to_bigg(lipidomic,conversion):
    """
    This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
    with their relative abundances.

    :param lipidomic: a dataframe that contains original lipid names and relative abundances

    :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers

    :return: a dictionary containing BiGG identifiers and their relative abundances
    """
    import pandas as pd
    #Generate the dictionary of lipid_id and relative abundances
    df = pd.merge(left=lipidomic, right=conversion, left_on='identifiers',right_on='lipid_name',how='outer')
    #1- Verify if strictly no match (specifically doubles the length
    if len(df.index) == (len(lipidomic['identifiers'])+len(conversion['lipid_name'])):
        #Change column order in conversion file
        conversion.columns = ['model_id','lipid_name']
        df = pd.merge(left=lipidomic, right=conversion, left_on='identifiers',right_on='lipid_name',how='outer')
        if len(df.index) == (len(lipidomic['identifiers'])):
            pass
        else:
            raise Exception('No match was found between lipidomic data and conversion file, review your conversion file.')
    else:
        pass

    #2- Verify if mismatch
    if len(df.index) != (len(lipidomic['identifiers'])):
        raise Exception('Incomplete match, assign model identifier to each lipid identified in lipidomic dataset')

    #Combine abundances for metabolites
    df1 = pd.concat([df['model_id'], df['abundances']], axis=1)
    grouped = df1.groupby('model_id').agg(lambda x: sum(x))

    return dict(zip([i for i in grouped.index], [i for i in grouped.abundances]))

def _get_relative_abundance(bigg_abundance):
    # Calculate relative abundances
    total_peak = sum(bigg_abundance.values())
    return {k: float(v) / total_peak for k, v in bigg_abundance.items()}

def _get_lipid_weight(model,compound_list,R_WEIGHT):
    #Function to change the molecular weight of a lipid that entails R chains of undefined weight
    def change_mol_weight(metab, R_WEIGHT):
        # R_WEIGHT is the estimated weight of a fatty acid chain R of length 18 carbons
        formula_list = []
        number_of_R = []
        # Remove R from formula
        for element in metab.formula:
            if element != 'R':
                formula_list.append(element)
            elif element == 'R':
                number_of_R.append('R')

        r_weight = R_WEIGHT * len(number_of_R)
        new_formula = ''.join(formula_list)
        metab.formula = new_formula
        total_weight = metab.formula_weight + r_weight
        return total_weight

    keys,values = [],[]
    for i in compound_list:
        keys.append(i)
        # Get metabolite in BiGG ID
        model_metab = model.metabolites.get_by_id(i)
        # Get molecular weight of the compound from chemical formula
        try:
            values.append(model_metab.formula_weight)

        except:
            print('Metabolite %s has no weight. Verifying if R in formula' % (model_metab.id,))
            if 'R' in model_metab.formula:
                print('%s formula is %s and contains R. Attributing R_WEIGHT' % (
                    model_metab.name, model_metab.formula))
                values.append(change_mol_weight(model_metab,R_WEIGHT))
            else:
                raise Exception('No R chain was found in that compound, please update your the molecular weight or formula of %s'
                      % (model_metab.id,))

    return dict(zip(keys,values))

def _calculate_coefficients(weight_dict,relative_abundance,LIPID_WEIGHT,CELL_WEIGHT,model):
    keys,values = [],[]
    for k,v in weight_dict.items():
        #Generate the total weight of the compound in the cell
        total_weight = LIPID_WEIGHT * relative_abundance.get(k)
        # Get molarity of the compound in the cell
        mmols_per_cell = (total_weight / v) * 1000
        # Convert to usable units in BIOMASS equation
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        # Save value
        values.append(-mmols_per_gDW)
        keys.append(model.metabolites.get_by_id(k))

    return dict(zip(keys,values))

#User accessible API functions
def generate_coefficients(path_to_lipidomic,path_to_conversion_file,
                     path_to_model,
                     LIPID_WEIGHT_FRACTION=0.091,
                     R_WEIGHT = 284.486):
    """

    Generates a dictionary of metabolite:coefficients for the lipid content of the cell. Lipids vary from a specie to another.
    The lipidomic data provides the relative abundance of each lipid specie while the to_bigg_dict allows to convert identifiers given in the lipidomic data to BiGG identifiers for which
    the metabolite weight is known and can be added easily to the biomass.

    :param path_to_lipidomic: a dataframe of metabolites identified in the lipidomic experiment

    :param path_to_conversion_file: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param LIPID_RATIO: measured lipid ratio of the cell, otherwise default

    :param R_WEIGHT: weight of a carbon chain, otherwise default. If the weight of the lipid is not known it will be inferred based on the number of R chains and this given weight.

    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.

    """
    CELL_WEIGHT = 280
    if LIPID_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    # Operation 0.1
    #Get the total lipid weight in the cell
    LIPID_WEIGHT = LIPID_WEIGHT_FRACTION * CELL_WEIGHT

    #1- Import model
    model = _import_model(path_to_model)
    #2- Import lipidomic and conversion
    lipidomic = _import_lipidomic(path_to_lipidomic)
    conversion = _import_conversion(path_to_conversion_file)
    #3- Generate a dictionary of BiGG IDs and relative abundances
    bigg_abundance = _convert_lipidomics_to_bigg(lipidomic, conversion)
    #4- Get the relative abundance of each lipid
    rel_abundance = _get_relative_abundance(bigg_abundance)
    #5- Get the weight of each lipid specie
    weight_dict = _get_lipid_weight(model,bigg_abundance.keys(),R_WEIGHT)
    #6- Calculate biomass stoichiometric coefficient
    biomass_coefficients = _calculate_coefficients(weight_dict,rel_abundance,LIPID_WEIGHT,CELL_WEIGHT,model)

    return biomass_coefficients

def filter_for_model_lipid(path_to_conversion_file, path_to_model):
    """
    Finds lipids identified through manual curation in model.

    :param path_to_conversion_file: a two column csv file converting from the name present in the lipidomic data to
    BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param path_to_model: a path to the model, format supported are json and xml

    :return: updated dataframe with metabolites found in the model
    """
    import pandas as pd
    # Get the model
    model = _import_model(path_to_model)
    # Get the metabolites in model
    model_metab_id = [m.id for m in model.metabolites]
    # Get to_bigg_dict
    to_bigg_df = pd.read_csv(path_to_conversion_file)
    to_bigg_dict = dict(zip([i for i in to_bigg_df[to_bigg_df.columns[0]]],
                            [i for i in to_bigg_df[to_bigg_df.columns[1]]]))

    # Get the metabolites that are in the model
    model_metab = {k: v for k, v in to_bigg_dict.items() if v in model_metab_id}

    # Get the metabolites that are not in the model but present in OMICs data
    non_model_metab = [k for k,v in to_bigg_dict.items() if v not in model_metab_id]
    if len(non_model_metab) != 0:
        print("These lipids were not found in the model but were present in your lipidomic data, "
                     "consider adding them to your model: %s " % ([metab for metab in non_model_metab]))

    model_metab_df = pd.DataFrame({'lipid_name':model_metab.keys(),'lipid_id':model_metab.values()},columns=['lipid_name','lipid_id'])

    return model_metab_df

def update_biomass_coefficients(dict_of_coefficients,model):
    """

    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.

    """
    from BOFdat import update
    update.update_biomass(dict_of_coefficients, model)
