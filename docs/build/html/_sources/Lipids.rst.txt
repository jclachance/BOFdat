Lipids
======

"""
    Lipids vary from a specie to another. The lipidomic data provides the relative abundance of each lipid specie
    while the to_bigg_dict allows to convert identifiers given in the lipidomic data to BiGG identifiers for which
    the metabolite weight is known and can be added easily to the biomass.
    =====================
    Parameters
    :param lipidomic: a dictionary or dataframe of metabolites identified in the lipidomic experiment
    :param to_bigg_dict: a dictionary converting from the name present in the lipidomic data to BiGG identifiers.
    This dictionary is generated through manual curation from the modeller.
    :param CELL_WEIGHT: measured cell weight in femtograms, otherwise default
    :param LIPID_RATIO: measured lipid ratio of the cell, otherwise default
    :param R_WEIGHT: weight of a carbon chain, otherwise default. If the weight of the lipid is not known it will be
    inferred based on the number of R chains and this given weight.
    =====================
    Return
    :return: a dictionary of metabolites and coefficients that can be used to update the biomass objective function.
"""
