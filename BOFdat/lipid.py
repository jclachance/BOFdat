"""
Lipid
=====

This module generates BOFsc for the lipid content of the cell.

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

def filter_for_model_lipid(path_to_conversion_file, path_to_model):
    """

    :param path_to_conversion_file: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

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
    model_metab = {k: v for k, v in to_bigg_dict.iteritems() if v in model_metab_id}

    # Get the metabolites that are not in the model but present in OMICs data
    non_model_metab = [k for k,v in to_bigg_dict.iteritems() if v not in model_metab_id]
    if len(non_model_metab) != 0:
        print("These lipids were not found in the model but were present in your lipidomic data, "
                     "consider adding them to your model: %s " % ([metab for metab in non_model_metab]))

    model_metab_df = pd.DataFrame({'lipid_name':model_metab.keys(),'lipid_id':model_metab.values()},columns=['lipid_name','lipid_id'])

    return model_metab_df

def generate_coefficients(path_to_lipidomic,path_to_conversion_file,
                     path_to_model,
                     CELL_WEIGHT=280,
                     LIPID_WEIGHT_FRACTION=0.091,
                     R_WEIGHT = 284.486):
    """

    Generates a dictionary of metabolite:coefficients for the lipid content of the cell. Lipids vary from a specie to another.
    The lipidomic data provides the relative abundance of each lipid specie while the to_bigg_dict allows to convert identifiers given in the lipidomic data to BiGG identifiers for which
    the metabolite weight is known and can be added easily to the biomass.

    :param lipidomic: a dictionary or dataframe of metabolites identified in the lipidomic experiment

    :param path_to_conversion_file: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param CELL_WEIGHT: measured cell weight in femtograms, otherwise default

    :param LIPID_RATIO: measured lipid ratio of the cell, otherwise default

    :param R_WEIGHT: weight of a carbon chain, otherwise default. If the weight of the lipid is not known it will be inferred based on the number of R chains and this given weight.

    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.

    """
    if LIPID_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    # Operation 0.1
    #Get the total lipid weight in the cell
    LIPID_WEIGHT = LIPID_WEIGHT_FRACTION * CELL_WEIGHT

    # Operation 0.2
    def make_compliant_lipidomic(path_to_lipidomic):
        import pandas as pd
        return pd.read_csv(path_to_lipidomic,names=['lipid_name', 'abundance'],skiprows=1)

    def make_compliant_bigg(path_to_bigg_dict):
        import pandas as pd
        return pd.read_csv(path_to_bigg_dict, names=['lipid_name','lipid_id'],skiprows=1)

    # Operation 1
    def convert_lipidomics_to_bigg(lipid_abun,lipid_conv):
        """
        This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
        with their relative abundances.

        :param lipidomic: a two column csv file that contains original lipid names and relative abundances

        :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers

        :return: a dictionary containing BiGG identifiers and their relative abundances
        """
        import pandas as pd
        #Generate the dictionary of lipid_id and relative abundances
        df = pd.merge(left=lipid_conv, right=lipid_abun, on='lipid_name')
        df1 = pd.concat([df.lipid_id, df.abundance], axis=1)
        grouped = df1.groupby('lipid_id').agg(lambda x: sum(x))

        return dict(zip([i for i in grouped.index], [i for i in grouped.abundance]))

    # Operation 2
    def get_relative_abundance(bigg_abundance):
        # Calculate relative abundances
        total_peak = sum(bigg_abundance.values())
        return {k: float(v) / total_peak for k, v in bigg_abundance.iteritems()}

    # Operation 3
    def get_lipid_weight(model,compound_list,R_WEIGHT):
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

            try:
                # Get molecular weight of the compound from chemical formula
                values.append(model_metab.formula_weight)

            except:
                print('metabolite has no weight %s. Verifying if R in formula' % (model_metab.id,))
                if 'R' in model_metab.formula:
                    print('%s formula is %s and contains R. Attributing R chain weight to 284.486' % (
                        model_metab.name, model_metab.formula))
                    values.append(change_mol_weight(model_metab,R_WEIGHT))

                else:
                    print('no R chain was found in that compound, please update your the molecular weight of %s'
                          % (model_metab.id,))

        return dict(zip(keys,values))

    # Operation 4
    def calculate_coefficients(weight_dict,relative_abundance,LIPID_WEIGHT,CELL_WEIGHT,model):
        keys,values = [],[]
        for k,v in weight_dict.iteritems():
            #Generate the total weight of the compound in the cell
            total_weight = LIPID_WEIGHT * relative_abundance.get(k)
            # Get molarity of the compound in the cell
            mmols_per_cell = (total_weight / v) * 1000
            # Convert to usable units in BIOMASS equation
            mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
            # Save value
            values.append(mmols_per_gDW)
            keys.append(model.metabolites.get_by_id(k))

        return dict(zip(keys,values))

    #Operations
    #0.1- Calculate the total lipid weight
    #0.2- Make data compliant for the rest of the functions
    lipidomic_compliant = make_compliant_lipidomic(path_to_lipidomic)
    bigg_compliant = make_compliant_bigg(path_to_conversion_file)

    #0.3- Get the model
    model = _import_model(path_to_model)
    #1- Generate a dictionary of BiGG IDs and relative abundances
    bigg_abundance = convert_lipidomics_to_bigg(lipidomic_compliant, bigg_compliant)
    #2- Get the relative abundance of each lipid
    rel_abundance = get_relative_abundance(bigg_abundance)
    #3- Get the weight of each lipid specie
    weight_dict = get_lipid_weight(model,bigg_abundance.keys(),R_WEIGHT)
    #4- Calculate biomass stoichiometric coefficient
    biomass_coefficients = calculate_coefficients(weight_dict,rel_abundance,LIPID_WEIGHT,CELL_WEIGHT,model)

    return biomass_coefficients

def update_biomass_coefficients(dict_of_coefficients,model):
    """

    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.

    """
    from BOFdat import update
    update.update_biomass(dict_of_coefficients, model)