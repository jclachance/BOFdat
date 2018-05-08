"""
RNA
===

This module generates BOFsc for the 4 bases of RNA (ATP, UTP, CTP and GTP)

"""
BASES = ['A', 'U', 'C', 'G']
import pandas as pd
from IPython import embed
# Methods
def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _import_transcriptomic(path_to_transcriptomic,all_locus):
    import pandas as pd
    import warnings
    transcriptomic =pd.read_csv(path_to_transcriptomic,header=None)
    #1- Verify number of columns
    if len(transcriptomic.columns) > 2:
        raise Exception("Your file format is not appropriate, more than 2 columns")
    #2- Verify presence of header
    if type(transcriptomic.loc[0, 0]) == str and type(transcriptomic.loc[0, 1]) == str:
        transcriptomic = transcriptomic.iloc[1:]
    #3- Remove null data
    if transcriptomic.isnull().values.any():
        transcriptomic = transcriptomic.dropna()
    #4- Verify column order (identifiers first, abundance second)
    try:
        try:
            abundances = [float(i) for i in transcriptomic.iloc[0:, 1]]
            identifiers = [str(i) for i in transcriptomic.iloc[0:, 0]]
        except:
            abundances = [float(i) for i in transcriptomic.iloc[0:, 0]]
            identifiers = [str(i) for i in transcriptomic.iloc[0:, 1]]
    except Exception:
        raise Exception('The abundances cannot be converted to float.')

    conform_df = pd.DataFrame({'identifiers':identifiers,'abundances':abundances},columns=['identifiers','abundances'])
    #5- Verify redundancy in protein identifiers
    if len(set(conform_df['identifiers'])) == len(conform_df['identifiers']):
        pass
    else:
        warnings.warn('Redundancy in dataset identifiers')
    #6- Make sure that locus_tag or gene ID are used
    if list(set(conform_df['identifiers']).intersection(set(all_locus))) > 0:
        pass
    else:
        raise Exception("Identifiers not 'locus_tag' or 'GeneID'")
    #7- Verify if given identifiers are provided in GenBank file
    if list(set(conform_df['identifiers']).intersection(set(all_locus))) == len(conform_df):
        pass
    else:
        warnings.warn('Some identifiers not found in provided annotation')
    return conform_df

def _get_number(seq):
    # This function is used if and only if the level of expression of each gene is assumed equal
    ratio_gene = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    for base in BASES:
        number = float(seq.count(base))
        ratio_gene[base] = number
    return ratio_gene

def _get_fraction(seq):
    # This function is to be used in case the level of expression of each gene cannot be assumed equal
    ratio_gene = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    for base in BASES:
        fraction = float(seq.count(base)) / len(seq)
        ratio_gene[base] = fraction
    return ratio_gene

def _get_RNA_sequence(location, strand,genome_record):
    # This function spits out the RNA sequence of a given gene
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    # Get the gene sequence
    sequence = []
    for number in range(location.start, location.end):
        sequence.append(genome_record[number])
    seq_str = ''.join(sequence)
    my_seq = Seq(seq_str, IUPAC.unambiguous_dna)
    try:
        if strand == 1:
            rna_seq = my_seq.transcribe()
        elif strand == -1:
            rna_seq = my_seq.reverse_complement().transcribe()
    except:
        raise Exception('No strand provided')

    return str(rna_seq)

def _make_number_df(number_list, locus_list, seq_list):
    # This function generates a dataframe with locus info, sequence and amount of each base
    A_number, U_number, C_number, G_number = [], [], [], []
    for d in number_list:
        for base in BASES:
            value = d.get(base)
            if base == 'A':
                A_number.append(value)
            elif base == 'U':
                U_number.append(value)
            elif base == 'C':
                C_number.append(value)
            elif base == 'G':
                G_number.append(value)

    generic_df = pd.DataFrame({'locus': locus_list, 'sequence': seq_list,
                               'A': A_number, 'U': U_number,
                               'C': C_number, 'G': G_number},
                              columns=['locus', 'sequence', 'A', 'U', 'C', 'G'])
    return generic_df

def _get_total_fractions(df):
    total_per_base = df.sum(axis=0, numeric_only=True)
    grand_total = float(total_per_base.sum(axis=0))

    fraction_dict = {'A': 0, 'U': 0, 'G': 0, 'C': 0}
    fraction_dict['A'] = total_per_base.A / grand_total
    fraction_dict['U'] = total_per_base.U / grand_total
    fraction_dict['G'] = total_per_base.G / grand_total
    fraction_dict['C'] = total_per_base.C / grand_total
    return fraction_dict

def _get_mRNA_fractions(df, path_to_transcriptomic,all_locus):
    transcriptomic = _import_transcriptomic(path_to_transcriptomic,all_locus)
    # Merge dataframes
    mean_abundance = pd.merge(left=df, right=transcriptomic, left_on='locus', right_on='identifiers')
    # Generate list of normalized values per gene per base by RPKM from transcriptomic data
    A_norm, U_norm, G_norm, C_norm = [], [], [], []
    for i, row in mean_abundance.iterrows():
        # 1- Multiply fraction by abundance for each gene
        A_mean = row['A'] * row['abundances']
        U_mean = row['U'] * row['abundances']
        G_mean = row['G'] * row['abundances']
        C_mean = row['C'] * row['abundances']
        A_norm.append(A_mean)
        U_norm.append(U_mean)
        G_norm.append(G_mean)
        C_norm.append(C_mean)

    # Get the mean of each list stored in a dictionary
    mRNA_fractions = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    mRNA_fractions['A'] = sum(A_norm) / sum(mean_abundance['abundances'])
    mRNA_fractions['U'] = sum(U_norm) / sum(mean_abundance['abundances'])
    mRNA_fractions['G'] = sum(G_norm) / sum(mean_abundance['abundances'])
    mRNA_fractions['C'] = sum(C_norm) / sum(mean_abundance['abundances'])

    return mRNA_fractions

def _process_record(path_to_genbank,path_to_transcriptomic,identifier):
    from Bio import SeqIO
    rRNA_seq, rRNA_locus, rRNA_number = [], [], []
    tRNA_seq, tRNA_locus, tRNA_number = [], [], []
    mRNA_seq, mRNA_locus, mRNA_fraction = [], [], []

    if identifier == 'locus_tag':
        genome_record = SeqIO.parse(path_to_genbank, 'genbank')
        for record in genome_record:
            for i,element in enumerate(record.features):
                if element.type == 'CDS':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    mRNA_seq.append(rna_seq)
                    mRNA_locus.append(element.qualifiers['locus_tag'][0])
                    mRNA_fraction.append(_get_fraction(rna_seq))

                if element.type == 'tRNA':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    tRNA_seq.append(rna_seq)
                    tRNA_number.append(_get_number(rna_seq))
                    try: tRNA_locus.append(element.qualifiers['locus_tag'][0])
                    except:
                        # It is possible that tRNA features lack the 'locus_tag' field!
                        # lines below are taken from Prokaryotic Annotation Guide!
                        # (https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/#RNA)
                        #
                        # RNA features (rRNA, tRNA, ncRNA) must include a corresponding gene feature with a locus_tag qualifier. 
                        # Please be sure to specify which amino acid the tRNA gene corresponds to. 
                        # If the amino acid of a tRNA is unknown, use tRNA-Xxx as the product, as in the example.
                        # Many submitters like to label the tRNAs such as tRNA-Gly1, etc. 
                        # If you wish to do this please include "tRNA-Gly1" as a note and not in /gene. 
                        # The use of /gene is reserved for the actual biological gene symbol such as "trnG".
                        # If a tRNA is a pseudogene, please use the /pseudo qualifier.
                        before = record.features[max(0,i-1)]
                        if before.type == 'gene':
                            try: tRNA_locus.append(before.qualifiers['locus_tag'][0])
                            except: raise Exception("Can't fetch tRNA locus_tag!")
                        else:
                            # raise Exception("Can't locate gene corresponding to a tRNA!") #TODO add a warning?
                            tRNA_locus.append("tRNA_%s" %(len(tRNA_locus)))

                if element.type == 'rRNA':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    rRNA_seq.append(rna_seq)
                    rRNA_number.append(_get_number(rna_seq))
                    try: rRNA_locus.append(element.qualifiers['locus_tag'][0])
                    except:
                        before = record.features[max(0,i-1)]
                        if before.type == 'gene':
                            try: rRNA_locus.append(before.qualifiers['locus_tag'][0])
                            except: raise Exception("Can't fetch rRNA locus_tag!")
                        else:
                            # raise Exception("Can't locate gene corresponding to a rRNA!") #TODO add a warning?
                            rRNA_locus.append("rRNA_%s" %(len(rRNA_locus)))



    elif identifier == 'geneID':
        genome_record = SeqIO.parse(path_to_genbank, 'genbank')
        for record in genome_record:
            for element in record.features:
                if element.type == 'mRNA':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    mRNA_seq.append(rna_seq)
                    mRNA_locus.append(j[7:] for j in element.qualifiers['db_xref'] if j.startswith('GeneID:'))
                    mRNA_fraction.append(_get_fraction(rna_seq))

                if element.type == 'tRNA':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    tRNA_seq.append(rna_seq)
                    tRNA_locus.append(j[7:] for j in element.qualifiers['db_xref'] if j.startswith('GeneID:'))
                    tRNA_number.append(_get_number(rna_seq))

                if element.type == 'rRNA':
                    rna_seq = _get_RNA_sequence(element.location, element.strand, record)
                    rRNA_seq.append(rna_seq)
                    rRNA_locus.append(j[7:] for j in element.qualifiers['db_xref'] if j.startswith('GeneID:'))
                    rRNA_number.append(_get_number(rna_seq))
    else:
        raise Exception("Unsupported identifier, use 'locus_tag' or 'GeneID'")

    # Make dataframe for each
    tRNA_df = _make_number_df(tRNA_number, tRNA_locus, tRNA_seq)
    rRNA_df = _make_number_df(rRNA_number, rRNA_locus, rRNA_seq)
    mRNA_df = _make_number_df(mRNA_fraction, mRNA_locus, mRNA_seq)
    rRNA_dict = _get_total_fractions(rRNA_df)
    tRNA_dict = _get_total_fractions(tRNA_df)
    all_locus = [i for i in rRNA_df['locus']] + [i for i in tRNA_df['locus']] + [i for i in mRNA_df['locus']]
    mRNA_dict = _get_mRNA_fractions(mRNA_df, path_to_transcriptomic,all_locus)

    return rRNA_dict, tRNA_dict, mRNA_dict

def _total_coefficients(mRNA_fractions, tRNA_fractions, rRNA_fractions, mRNA_RATIO, tRNA_RATIO, rRNA_RATIO):
    # Multiply mRNA,rRNA and tRNA
    # dict values by the coefficients
    # to get the average of each base from RNA in the cell
    RNA_total = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    for letter in BASES:
        mrna = mRNA_fractions[letter] * mRNA_RATIO
        trna = tRNA_fractions[letter] * tRNA_RATIO
        rrna = rRNA_fractions[letter] * rRNA_RATIO
        RNA_total[letter] = mrna + trna + rrna

    return RNA_total

def _convert_to_mmolgDW(RNA_coefficients, model, RNA_RATIO, CELL_WEIGHT):
    DIPHOSPHATE_WEIGHT = 174.951262
    # Get coefficients for BIOMASS
    # Transform the ratios into mmol/gDW
    RNA_WEIGHT = CELL_WEIGHT * RNA_RATIO

    rna_base_to_bigg = {'A': model.metabolites.atp_c, 'U': model.metabolites.utp_c,
                        'C': model.metabolites.ctp_c, 'G': model.metabolites.gtp_c}
    metabolites, coefficients = [], []
    # Get the total weight of each letter
    for letter in BASES:
        ratio = RNA_coefficients.get(letter)
        total_weight = ratio * RNA_WEIGHT
        metab = rna_base_to_bigg.get(letter)
        mol_weight = metab.formula_weight - DIPHOSPHATE_WEIGHT
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        coefficients.append(mmols_per_gDW)
        metabolites.append(rna_base_to_bigg.get(letter))

    RNA_biomass_ratios = dict(zip(metabolites, [-i for i in coefficients]))
    return RNA_biomass_ratios

def generate_coefficients(path_to_genbank, path_to_model, path_to_transcriptomic,
                         RNA_WEIGHT_FRACTION=0.205,
                         rRNA_WEIGHT_FRACTION=0.9,
                         tRNA_WEIGHT_FRACTION=0.05,
                         mRNA_WEIGHT_FRACTION=0.05,
                         identifier='locus_tag'):
    """
    Generates a dictionary of metabolite:coefficients for the 4 RNA bases from the organism's
    GenBank annotated file, total RNA weight percentage, transcriptomic. Alternately, ribosomal,
    transfer and messenger RNA relative abundances can be incorporated otherwise the default 80% rRNA, 10% tRNA and
    10% mRNA are used.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param path_to_transcriptomic: a two column pandas dataframe (gene_id, abundance)

    :param RNA_WEIGHT_FRACTION: the weight fraction of RNA in the entire cell

    :param rRNA_WEIGHT_FRACTION: the fraction of rRNA to total

    :param tRNA_WEIGHT_FRACTION: the fraction of tRNA to total

    :param mRNA_WEIGHT_FRACTION: the fraction of mRNA to total

    :param identifier: the type of identifier in the input file, 'locus_tag' or 'geneID'

    :return: a dictionary of metabolites and coefficients
    """
    CELL_WEIGHT = 280
    if RNA_WEIGHT_FRACTION > 1. or rRNA_WEIGHT_FRACTION > 1. or tRNA_WEIGHT_FRACTION > 1. or mRNA_WEIGHT_FRACTION > 1.:
        raise Exception('WEIGHT FRACTION should be a number between 0 and 1')
    # Operations
    model = _import_model(path_to_model)
    try: rRNA_dict, tRNA_dict, mRNA_dict = _process_record(path_to_genbank,path_to_transcriptomic,identifier)
    except: embed()
    RNA_coefficients = _total_coefficients(rRNA_dict, tRNA_dict, mRNA_dict,
                                           mRNA_WEIGHT_FRACTION, tRNA_WEIGHT_FRACTION, rRNA_WEIGHT_FRACTION)
    RNA_biomass_ratios = _convert_to_mmolgDW(RNA_coefficients,
                                            model, RNA_WEIGHT_FRACTION, CELL_WEIGHT)
    ppi_coeff = sum(RNA_biomass_ratios.values())
    ppi_dict = {model.metabolites.get_by_id('ppi_c'): ppi_coeff}
    RNA_biomass_ratios.update(ppi_dict)

    return RNA_biomass_ratios

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
    from BOFdat import update
    update.update_biomass(dict_of_coefficients,model)
