"""
Lipid
=====

This module generates BOFsc for the lipid content of the cell.

"""

def get_coefficients(path_to_lipidomic,path_to_bigg_dict,
                     path_to_model,
                     CELL_WEIGHT=280,
                     LIPID_RATIO=0.091,
                     R_WEIGHT = 284.486):
    """

    Generates a dictionary of metabolite:coefficients for the lipid content of the cell. Lipids vary from a specie to another.
    The lipidomic data provides the relative abundance of each lipid specie while the to_bigg_dict allows to convert identifiers given in the lipidomic data to BiGG identifiers for which
    the metabolite weight is known and can be added easily to the biomass.

    :param lipidomic: a dictionary or dataframe of metabolites identified in the lipidomic experiment

    :param to_bigg_dict: a dictionary converting from the name present in the lipidomic data to BiGG identifiers. This dictionary is generated through manual curation from the modeller.

    :param CELL_WEIGHT: measured cell weight in femtograms, otherwise default

    :param LIPID_RATIO: measured lipid ratio of the cell, otherwise default

    :param R_WEIGHT: weight of a carbon chain, otherwise default. If the weight of the lipid is not known it will be inferred based on the number of R chains and this given weight.

    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.

    """

    # Operation 0.1
    #Get the total lipid weight in the cell
    LIPID_WEIGHT = LIPID_RATIO * CELL_WEIGHT

    # Operation 0.2
    def make_compliant_lipidomic(path_to_lipidomic):
        import pandas as pd
        return pd.read_csv(path_to_lipidomic,names=['lipid_name', 'abundance'],skiprows=1)

    def make_compliant_bigg(path_to_bigg_dict):
        import pandas as pd
        df = pd.read_csv(path_to_bigg_dict, names=['lipid_name','lipid_id'],skiprows=1)
        keys = [i for i in df.lipid_name]
        values = [i for i in df.lipid_id]
        return dict(zip(keys,values))

    #Operation 0.3
    def import_model(path_to_model):
        import cobra
        extension = path_to_model.split('.')[-1]
        if extension == 'json':
            model = cobra.io.load_json_model(path_to_model)
        elif extension == 'xml':
            model = cobra.io.read_sbml_model(path_to_model)
        return model

    # Operation 1
    def convert_lipidomics_to_bigg(lipidomic,to_bigg_dict):
        """
        This function generates a dictionary of BiGG identifiers that were generated through manual curation of the user
        with their relative abundances.
        :param lipidomic: a two column csv file that contains original lipid names and relative abundances
        :param to_bigg_dict: a dictionary converting names in the experimental data to BiGG identifiers
        :return: a dictionary containing BiGG identifiers and their relative abundances
        """
        import pandas as pd
        #Generate the dictionary
        keys,values = [],[]
        for i,row in lipidomic.iterrows():
            keys.append(to_bigg_dict.get(row.lipid_name))
            values.append(row.abundance)

        return dict(zip(keys,values))

    # Operation 2
    def get_relative_abundance(bigg_abundance):
        # Calculate relative abundances
        total_peak = sum(bigg_abundance.values())
        return {k: v / total_peak for k, v in bigg_abundance.iteritems()}

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
    bigg_compliant = make_compliant_bigg(path_to_bigg_dict)
    #0.3- Get the model
    model = import_model(path_to_model)
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
    from biomass import Update
    Update.update_biomass(dict_of_coefficients, model)