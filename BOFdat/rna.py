"""
RNA
===

This module generates BOFsc for the 4 bases of RNA (ATP, UTP, CTP and GTP)

"""
BASES = ['A', 'U', 'C', 'G']
import pandas as pd
# Methods
def _import_genome(path_to_genbank):
    # This function import the genome record from genbank
    from Bio import SeqIO
    genome_record = SeqIO.read(path_to_genbank, 'genbank')
    return genome_record


def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        model = cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        model = cobra.io.read_sbml_model(path_to_model)
    return model


def _import_transcriptomic(path_to_transcriptomic):
    import pandas as pd
    return pd.read_csv(path_to_transcriptomic, names=['gene_ID', 'Mean'], skiprows=1)


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
        raise ('error no strand provided')

    return str(rna_seq)


def _make_number_df(number_list, locus_list, seq_list):

    # This function gnerates a dataframe with locus info, sequence and amount of each base
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


def _get_mRNA_fractions(df, path_to_transcriptomic):
    transcriptomic = _import_transcriptomic(path_to_transcriptomic)
    # Merge dataframes
    mean_abundance = pd.merge(left=df, right=transcriptomic, left_on='locus', right_on='gene_ID')
    # Generate list of normalized values per gene per base by RPKM from transcriptomic data
    A_norm, U_norm, G_norm, C_norm = [], [], [], []
    for i, row in mean_abundance.iterrows():
        # 1- Multiply fraction by abundance for each gene
        A_mean = row.A * row.Mean
        U_mean = row.U * row.Mean
        G_mean = row.G * row.Mean
        C_mean = row.C * row.Mean
        A_norm.append(A_mean)
        U_norm.append(U_mean)
        G_norm.append(G_mean)
        C_norm.append(C_mean)

    # Get the mean of each list stored in a dictionary
    mRNA_fractions = {'A': 0, 'U': 0, 'C': 0, 'G': 0}
    mRNA_fractions['A'] = sum(A_norm) / sum(mean_abundance.Mean)
    mRNA_fractions['U'] = sum(U_norm) / sum(mean_abundance.Mean)
    mRNA_fractions['G'] = sum(G_norm) / sum(mean_abundance.Mean)
    mRNA_fractions['C'] = sum(C_norm) / sum(mean_abundance.Mean)

    return mRNA_fractions


def _process_record(genome_record,path_to_transcriptomic):
    rRNA_seq, rRNA_locus, rRNA_number = [], [], []
    tRNA_seq, tRNA_locus, tRNA_number = [], [], []
    mRNA_seq, mRNA_locus, mRNA_fraction = [], [], []

    for element in genome_record.features:
        if element.type == 'CDS':
            rna_seq = _get_RNA_sequence(element.location, element.strand,genome_record)
            mRNA_seq.append(rna_seq)
            mRNA_locus.append(element.qualifiers['locus_tag'][0])
            mRNA_fraction.append(_get_fraction(rna_seq))

        if element.type == 'tRNA':
            rna_seq = _get_RNA_sequence(element.location, element.strand,genome_record)
            tRNA_seq.append(rna_seq)
            tRNA_locus.append(element.qualifiers['locus_tag'][0])
            tRNA_number.append(_get_number(rna_seq))

        if element.type == 'rRNA':
            rna_seq = _get_RNA_sequence(element.location, element.strand,genome_record)
            rRNA_seq.append(rna_seq)
            rRNA_locus.append(element.qualifiers['locus_tag'][0])
            rRNA_number.append(_get_number(rna_seq))
    # Make dataframe for each
    rRNA_dict = _get_total_fractions(_make_number_df(rRNA_number, rRNA_locus, rRNA_seq))
    tRNA_dict = _get_total_fractions(_make_number_df(tRNA_number, tRNA_locus, tRNA_seq))
    mRNA_dict = _get_mRNA_fractions(_make_number_df(mRNA_fraction, mRNA_locus, mRNA_seq), path_to_transcriptomic)

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
        mol_weight = metab.formula_weight
        mmols_per_cell = (total_weight / mol_weight) * 1000
        mmols_per_gDW = mmols_per_cell / CELL_WEIGHT
        coefficients.append(mmols_per_gDW)
        metabolites.append(rna_base_to_bigg.get(letter))

    RNA_biomass_ratios = dict(zip(metabolites, coefficients))
    return RNA_biomass_ratios


def generate_coefficients(path_to_genbank, path_to_model, path_to_transcriptomic,
                     CELL_WEIGHT=280, TOTAL_RNA_RATIO=0.205,
                     rRNA_RATIO=0.8, tRNA_RATIO=0.1, mRNA_RATIO=0.1):
    """

    Generates a dictionary of metabolite:coefficients for the 4 RNA bases from the organism's
    GenBank annotated file, total RNA weight percentage, transcriptomic. Alternately, ribosomal,
    transfer and messenger RNA relative abundances can be incorporated otherwise the default 80% rRNA, 10% tRNA and
    10% mRNA are used.

    :param path_to_genbank: a path to the GenBank annotation file of the organism, format should be compatible with BioPython SeqIO

    :param path_to_model: a path to the model, format supported are json and xml

    :param proteomics: a two column pandas dataframe (gene_id, abundance)

    :param CELL_WEIGHT: experimentally measured cell weight in femtograms, float

    :param RNA_RATIO: the ratio of DNA in the entire cell

    :return: a dictionary of metabolites and coefficients
    """
    # Imports
    import pandas as pd

    # Operations
    genome_record = _import_genome(path_to_genbank)
    rRNA_dict, tRNA_dict, mRNA_dict = _process_record(genome_record,path_to_transcriptomic)
    RNA_coefficients = _total_coefficients(rRNA_dict, tRNA_dict, mRNA_dict, mRNA_RATIO, tRNA_RATIO, rRNA_RATIO)
    RNA_biomass_ratios = _convert_to_mmolgDW(RNA_coefficients,
                                            _import_model(path_to_model), TOTAL_RNA_RATIO, CELL_WEIGHT)

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
