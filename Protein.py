"""
Protein
=======

This module generates BOFsc for the 20 amino acids contained in proteins.

"""

# Methods
def get_coefficients(path_to_genbank, path_to_model, path_to_proteomic, CELL_WEIGHT=280, PROTEIN_RATIO=0.55):
    """

    Generates a dictionary of metabolite:coefficients for the 20 amino acids contained in proteins from the organism's
    GenBank annotated file, total Protein weight percentage and proteomic data.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param path_to_proteomic: a two column pandas dataframe (gene_id, abundance)

    :param CELL_WEIGHT: experimentally measured cell weight in femtograms, float

    :param PROTEIN_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """

    import pandas as pd
    AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                   'Y']

    def import_genome(genbank):
        from Bio import SeqIO
        return SeqIO.read(genbank, 'genbank')

    def import_model(path_to_model):
        import cobra
        extension = path_to_model.split('.')[-1]
        if extension == 'json':
            return cobra.io.load_json_model(path_to_model)
        elif extension == 'xml':
            return cobra.io.read_sbml_model(path_to_model)

    def import_proteomic(path_to_proteomic):
        import pandas as pd
        return pd.read_csv(path_to_proteomic, names=['gene_ID','Mean'],skiprows=1)

    def get_protein_sequence(genome_record):

        # Get the sequence of each protein, locus tag and prot id from genbank file
        seq_list, protID_list, locus_list = [], [], []
        for element in genome_record.features:
            if element.type == 'CDS' and 'protein_id' in element.qualifiers:
                locus_list.append(element.qualifiers['locus_tag'][0])
                protID_list.append(element.qualifiers['protein_id'][0])
                seq_list.append(element.qualifiers['translation'][0])

        return pd.DataFrame({'Locus_tag': locus_list, 'Protein_ID': protID_list, 'Sequence': seq_list})

    def get_aa_composition(df1):

        # For each protein find the amino acid composition
        # Outputs a dictionnary of dictionnaries where:
        # Keys = locus_tag
        # Values = A dictionary for each amino acid
        # This dictionary contains:
        # Keys = amino acid by letter code
        # Values = the occurence of that amino acid
        list_of_dict = []
        for i, row in df1.iterrows():
            list_of_occurences = []
            # Get the occurence for each letter
            for letter in AMINO_ACIDS:
                protein_sequence = row.Sequence
                occurence_of_letter = protein_sequence.count(letter)
                list_of_occurences.append(occurence_of_letter)
            # Generate dictionary of occurences for a given gene
            dict_of_occurences = dict(zip(AMINO_ACIDS, list_of_occurences))
            # Generate dict for each gene
            dict_per_locus = {row.Locus_tag: dict_of_occurences}
            # Store the amount of each amino acid per gene in a list
            list_of_dict.append(dict_per_locus)

        return list_of_dict

    def make_coeff_dict(proteomics):
        keys = [k for k in proteomics.gene_ID]
        values = [v for v in proteomics.Mean]
        return dict(zip(keys,values))

    def normalize_aa_composition(list_of_dict,path_to_proteomic):
        # Normalize the value of each amino acid per protein following proteomic data
        normalized_dict = {'A': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'G': 0., 'H': 0., 'I': 0.,
                           'K': 0., 'L': 0., 'M': 0., 'N': 0., 'P': 0., 'Q': 0., 'R': 0., 'S': 0., 'T': 0., 'V': 0.,
                           'W': 0., 'Y': 0.}
        #Import proteomic data into dictionnary
        proteomics = import_proteomic(path_to_proteomic)
        coeff_dict = make_coeff_dict(proteomics)
        for d in list_of_dict:
            # Get the coefficient from proteomics
            text = str(list(d.keys())[0])
            coeff = coeff_dict.get(text)
            # If no protein abundance coefficient is 0.
            try:
                coeff_number = float(coeff)
            except:
                coeff_number = 0.

            # Multiply each amino acid by the coefficient
            amino_acids = list(d.values())
            for letter in AMINO_ACIDS:
                value = float(amino_acids[0].get(letter))
                # Update the normalized value
                normalized_value = value * coeff_number
                new_value = normalized_dict.get(letter) + normalized_value
                normalized_dict[letter] = new_value

        return normalized_dict

    # Operations
    # 1- Parse the genome, extract protein sequence, count and store amino acid composition of each protein
    genome_record = import_genome(path_to_genbank)
    df1 = get_protein_sequence(genome_record)
    list_of_dict = get_aa_composition(df1)
    normalized_dict = normalize_aa_composition(list_of_dict,path_to_proteomic)

    def get_norm_sum(normalized_dict):
        # 1- Sum normalized ratios
        norm_sum = 0.
        for letter in AMINO_ACIDS:
            value = normalized_dict.get(letter)
            norm_sum = value + norm_sum

        return norm_sum

    def get_ratio(normalized_dict, norm_sum, PROTEIN_RATIO):
        # 2- Divide letter to norm_sum to get ratio of each amino acid in the cell
        # based on proteomic data
        ratio_dict = {'A': 0., 'C': 0., 'D': 0., 'E': 0., 'F': 0., 'G': 0., 'H': 0., 'I': 0.,
                      'K': 0., 'L': 0., 'M': 0., 'N': 0., 'P': 0., 'Q': 0., 'R': 0., 'S': 0., 'T': 0., 'V': 0.,
                      'W': 0., 'Y': 0.}

        # Constant for the amount of protein in the cell
        PROTEIN_WEIGHT = CELL_WEIGHT * PROTEIN_RATIO
        for letter in AMINO_ACIDS:
            value = normalized_dict.get(letter)
            ratio = value / norm_sum
            # Convert ratios to grams
            converted_ratio = ratio * PROTEIN_WEIGHT
            ratio_dict[letter] = converted_ratio

        return ratio_dict

    def convert_to_coefficient(path_to_model, CELL_WEIGHT):
        model = import_model(path_to_model)
        # 3- Convert gram ratios to mmol/g Dry weight
        '''
        To verify that the normalized to grams to get to the total amount of protein
        (here everything is converted to grams instead of femto grams)
        '''
        letter_to_bigg = {'A': model.metabolites.ala__L_c, 'C': model.metabolites.cys__L_c,
                          'D': model.metabolites.asp__L_c, 'E': model.metabolites.glu__L_c,
                          'F': model.metabolites.phe__L_c,
                          'G': model.metabolites.gly_c, 'H': model.metabolites.his__L_c,
                          'I': model.metabolites.ile__L_c, 'K': model.metabolites.lys__L_c,
                          'L': model.metabolites.leu__L_c,
                          'M': model.metabolites.met__L_c, 'N': model.metabolites.asn__L_c,
                          'P': model.metabolites.pro__L_c, 'Q': model.metabolites.gln__L_c,
                          'R': model.metabolites.arg__L_c,
                          'S': model.metabolites.ser__L_c, 'T': model.metabolites.thr__L_c,
                          'V': model.metabolites.val__L_c, 'W': model.metabolites.trp__L_c,
                          'Y': model.metabolites.tyr__L_c}

        metabolites, coefficients = [],[]
        # Get number of moles from number of grams
        for letter in AMINO_ACIDS:
            metab = letter_to_bigg.get(letter)
            mol_weight = metab.formula_weight
            grams = ratio_dict.get(letter)
            mmols_per_cell = (grams / mol_weight) * 1000
            mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
            coefficients.append(mmols_per_gDW)
            metabolites.append(letter_to_bigg.get(letter))

        Protein_biomass_coefficients = dict(zip(metabolites,coefficients))
        return Protein_biomass_coefficients

    # 2- Get coefficients from experimental proteomics data
    # Proteomics data should come in a 2 columns standard format gene_id:abundance
    norm_sum = get_norm_sum(normalized_dict)
    ratio_dict = get_ratio(normalized_dict, norm_sum, PROTEIN_RATIO)
    biomass_coefficients = convert_to_coefficient(path_to_model, CELL_WEIGHT)

    return biomass_coefficients

'''
The option to update the coefficients of the metabolites in the biomass objective function is left to the user
'''
def update_biomass_coefficients(dict_of_coefficients,model):
    """

    Updates the biomass coefficients given the input metabolite:coefficient dictionary.

    :param dict_of_coefficients: dictionary of metabolites and coefficients

    :param model: model to update

    :return: The biomass objective function is updated.
    """
    from biomass import Update
    Update.update_biomass(dict_of_coefficients,model)