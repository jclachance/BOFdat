def make_gene_flux_vector(fluxes, model):
    THRESHOLD = 1e-16
    fluxes[fluxes > THRESHOLD] = 1
    fluxes[fluxes <= THRESHOLD] = 0
    rxn_vector = pd.DataFrame({'Reactions': [i for i in fluxes.index],
                               'Flux': [i for i in fluxes]})
    # Make a dataframe in which each gene is a row with its flux value
    genes, flux_states = [], []
    for i, row in rxn_vector.iterrows():
        for g in model.reactions.get_by_id(row['Reactions']).genes:
            genes.append(g.id)
            flux_states.append(row['Flux'])

    gene_vector = pd.DataFrame({'Genes': genes, 'Flux': flux_states})
    # Group the genes together and get the total flux count
    grouped = gene_vector.groupby('Genes')
    temp_df = grouped.agg({'Flux': ['sum']})
    temp_df.columns = temp_df.columns.droplevel(level=1)
    # If the count is different than 0, leave the flux to on-state (1)
    temp_df[temp_df.Flux > 0] = 1
    # Make a dataframe to be compatible with compare_flux function
    gene_flux_vector = pd.DataFrame()
    gene_flux_vector['Genes'] = [i for i in temp_df.index]
    gene_flux_vector['Flux'] = [i for i in temp_df.Flux]

    return gene_flux_vector


def compare_flux_vectors(original_flux, new_flux):
    from scipy.spatial.distance import hamming
    compare_df = pd.merge(right=original_flux, left=new_flux, on='Genes', how='inner')
    # Apply hamming distance on proteomic data
    u = [f for f in compare_df.Flux_x]
    v = [x for x in compare_df.Flux_y]

    # Apply hamming distance on 2 vectors
    return hamming(u, v)