"""
Group end goals
===============

This module takes the outputs of the genetic algorithm and groups them into clusters of metabolic end goals.

"""
import warnings
import os
from os.path import join
from os import listdir
from BOFdat.util import dijkstra
from BOFdat.util.update import _import_model
from BOFdat.util.update import _get_biomass_objective_function
import cobra
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
#Visualization librairies
from pylab import rcParams
import matplotlib.pyplot as plt
import seaborn as sb


def _distance_metric(distance_df,fitness):
    distance_df = distance_df[distance_df[1]<1000]
    metrics = (fitness,distance_df[1].sum(),distance_df[1].mean(),distance_df[1].median())

    return metrics

def _get_hof_data(outpath):
    hofs=[]
    for f in listdir(outpath):
        if f.startswith('hof'):
            try:
                hofs.append(pd.read_csv(join(outpath, f),index_col=0))
            except:
                warnings.warn('Unable to load file %s'%(f,))

    return hofs

def _make_list(m):
    l = []
    for i in m.split(','):
        i = i.replace("'",'')
        i = i.replace("[",'')
        i = i.replace("]",'')
        i = i.replace(" ",'')
        l.append(i)

    return l

def _filter_hofs(hofs, THRESHOLD):
    """
    Selects the best individuals accross HOF generated from multiple evolutions.
    HOF are combined together and an arbitrary threshold is set to select the best individuals.
    The purpose of this operation is to ensure downstream analysis using metabolite frequency
    is relevant (the metabolites with higher frequency allow to obtain better MCC value).

    :param hofs: list of HOF dataframes with individuals and their associated fitness
    :param THRESHOLD: the percentage best indviduals to select (default is 20%)
    :return:
    """
    #Pool the hall of fame together
    pooled_hof = pd.concat(hofs)
    sorted_pooled_hof = pooled_hof.sort_values('Fitness',ascending=False)
    #Take the 20% best individuals
    final_index = int(len(pooled_hof)*THRESHOLD)
    metrics_df = sorted_pooled_hof.iloc[0:final_index,]
    #Get the list of metabolites contained in the best all of fames
    best_bof = []
    for i in range(len(metrics_df)):
        l = _make_list(metrics_df.iloc[i,0])
        for m in l:
            best_bof.append(m)

    return best_bof

def _make_freq_df(best_bof):
    def occurence(metab, best_bof):
        return float(len([m for m in best_bof if m == metab])) / len(best_bof)
    #Generate the frequency of each metabolite in hall of fames
    frequency_data = []
    for m in set(best_bof):
        f_m = occurence(m, best_bof)
        frequency_data.append((m, f_m))
    #Convert to DataFrame and sort by frequency
    freq_df = pd.DataFrame.from_records(frequency_data,
                                        columns=['Metab', 'Frequency'],
                                        index=[m for m in set(best_bof)])
    freq_df.sort_values('Frequency',inplace=True)

    return freq_df

def _generate_result_matrix(freq_df,dist_matrix):
    #Select the metabolites based on their apparition frequency in the HALL of FAME
    THRESHOLD = freq_df.mean()[0]
    selected_metab = [m for m in freq_df[freq_df.iloc[:, 1] > THRESHOLD].iloc[:, 0]]
    all_metab = set(selected_metab)
    #Generate the reduced distance matrix used for clustering of metabolic end goals
    best_bof_dist = dist_matrix[dist_matrix.index.isin(all_metab)]
    transposed_bbd = best_bof_dist.T
    r = transposed_bbd[transposed_bbd.index.isin(all_metab)]
    result_matrix = r.replace(np.nan,1000.)

    return result_matrix

def _dbscan_clustering(result_matrix,eps):
    from sklearn.cluster import dbscan
    # Cluster the distance with DBSCAN clustering algorithm
    dbscan_result = dbscan(result_matrix, eps=eps, min_samples=1)
    cluster_dict = {result_matrix.index[i]: dbscan_result[1][i] for i in range(len(dbscan_result[1]))}
    clusters = []
    for k,v in cluster_dict.iteritems():
        clusters.append((v,k))
        
    cluster_df = pd.DataFrame.from_records(clusters, columns=['cluster', 'metab'])
    
    #Group metabolites into their respective clusters
    grouped = cluster_df.groupby('cluster')
    grouped_clusters = grouped.agg(lambda x: ', '.join(x))
    
    return grouped_clusters,cluster_dict

def _display_occurence(freq_df,fig_name):
    """
    Plot the frequency of apparition of every metabolites in the Hall of Fames

    :param freq_df: dataframe of metabolites and their associated frequency
    :param fig_name: the name of the figure to generate

    """
    from pylab import rcParams
    #Set general ploting parameters
    default_param = {'figsize':(10,20),
                         'sb_style':'whitegrid',
                         'linestyle':':',
                         'marker':'',
                         'baseline_color':'r'
                         }

    rcParams['figure.figsize'] = default_param.get('figsize')

    freq_df.sort_values('Frequency',inplace=True,ascending=True)
    freq_df.plot(kind='barh')
    THRESHOLD = freq_df.mean()[0]
    plt.plot((THRESHOLD,THRESHOLD),(0,len(freq_df)),
                 default_param.get('baseline_color')
                 ,label=default_param.get('baseline_label'),lw=2.5,alpha=0.8,linestyle='-')
    number_of_evol = 24
    plt.title('Metabolite frequency in best individuals over %s evolutions'%(number_of_evol,),fontsize=24)
    plt.xlabel('Metabolite',fontsize=20)
    plt.xticks(fontsize=9)
    plt.ylabel('Number of apparitions',fontsize=20)
    plt.savefig(fig_name)

    plt.show()

def _find(x,my_palette,grouped_clusters):
    #1- Find the cluster to which x belongs
    def find_cluster(x,grouped_clusters):
        for i in range(len(grouped_clusters)):
            if x in grouped_clusters.loc[i][0].split(', '):
                return grouped_clusters.loc[i].name

    cluster_id = find_cluster(x,grouped_clusters)
    return my_palette.get(cluster_id)

def _display_result_matrix(cluster_matrix, cluster_dict, freq_df, fig_name):
    """
    Plot the result of spatial clustering with DBSCAN.

    :param grouped_clusters:
    :param cluster_matrix: the result of f
    :param figname: the name of the figure to generate
    :return:
    """
    from pylab import rcParams
    # General ploting parameters
    default_param = {'figsize': (15, 15),
                     'sb_style': 'whitegrid',
                     'fontsize': 12}
    rcParams['figure.figsize'] = default_param.get('figsize')
    rcParams['font.size'] = default_param.get('fontsize')

    # Arrange the cluster matrix so that the clusters from DBSCAN are together
    # Add a cluster column to the matrix
    cluster_matrix['cluster'] = [cluster_dict.get(cluster_matrix.index[i]) for i in range(len(cluster_matrix))]

    # Calculate each cluster's frequency
    # The cluster frequency is the sum of all individual metabolite frequency
    # Convert frequency dataframe to dictionary
    freq_dict = {row['Metab']: row['Frequency'] for i, row in freq_df.iterrows()}

    # Sort to ensure the rows are in order of cluster
    cluster_matrix.sort_values('cluster', inplace=True)
    metab_frequency = [freq_dict.get(i) for i in cluster_matrix.index]
    cluster_matrix['metabolite frequency'] = metab_frequency

    # Make a dictionary that maps cluster identifier number to their sum frequency
    cluster_frequency_dict = cluster_matrix.groupby('cluster').sum()['metabolite frequency'].to_dict()
    cluster_matrix['cluster frequency'] = [cluster_frequency_dict.get(i) for i in cluster_matrix.cluster]

    # Organize matrix for ploting
    cluster_matrix.sort_values('cluster frequency', ascending=False, inplace=True)
    cluster_matrix.replace(1000., 20., inplace=True)
    new_index = [i for i in cluster_matrix.index]

    # Transpose and associate metabolites to their respective clusters
    transposed_cluster_matrix = cluster_matrix.T
    transposed_cluster_matrix.drop(['cluster frequency',
                                    'metabolite frequency', 'cluster'],
                                   axis=0,
                                   inplace=True)
    ploting_matrix = transposed_cluster_matrix.reindex(new_index)

    # Initialize figure
    fig, ax = plt.subplots()

    # Plot the heatmap
    sb.heatmap(ploting_matrix.T, cbar=False, cmap="mako")

    ax.xaxis.set_ticks_position('top')

    i = 0
    for tick_label in ax.xaxis.get_ticklabels():
        # Change the xticks labels
        # tick_label.set_color(reversed_row_colors[i])
        tick_label.set_rotation(90)
        tick_label.set_fontsize(18)
        i += 1

    i = 0
    for tick_label in ax.yaxis.get_ticklabels():
        # Change the xticks labels
        # tick_label.set_color(reversed_row_colors[i])
        tick_label.set_fontsize(18)
        i += 1

    fig_path = os.path.join('/home/jean-christophe/Documents/Doctorat/BOFdat_paper/Figures_Drive/Figure 5/', fig_name)
    # plt.savefig(fig_path)
    plt.show()

def _select_metabolites(freq_df,grouped_clusters,best_bof):
    #Add selection of the clusters above average
    step3 = []
    for i in grouped_clusters.metab:
        cluster = i.split(', ')
        mini_step3 = []
        for m in cluster:
            if m in best_bof:
                mini_step3.append(1)
            else:
                mini_step3.append(0)

        step3.append(sum(mini_step3))

    better_results_df = pd.DataFrame({'Step 3 clusters': step3})
    compare_cluster_ga = pd.concat([grouped_clusters, better_results_df], axis=1)

    bd3_unique = []
    for i, row in compare_cluster_ga.iterrows():
        mini_unique = [(m, freq_df.loc[m]['Frequency']) for m in row['metab'].split(', ')]
        cluster = pd.DataFrame.from_records(mini_unique, columns=['metab', 'freq'])
        cluster.index = cluster.metab
        del cluster['metab']
        for m in cluster.index[cluster['freq'] == cluster.max()[0]]:
            bd3_unique.append(m)

    return bd3_unique


def _select_metabolites_method2(freq_df, cluster_matrix, best_bof):
    """
    Function that selects the metabolites to add to the BOF
    Returns a list of metabolites
    """
    # 1- Convert frequency dataframe to dictionary
    freq_dict = {row['Metab']: row['Frequency'] for i, row in freq_df.iterrows()}

    # 2- Sort to ensure the rows are in order of cluster
    cluster_matrix.sort_values('cluster', inplace=True)
    metab_frequency = [freq_dict.get(i) for i in cluster_matrix.index]
    cluster_matrix['metabolite frequency'] = metab_frequency

    # 3- Make a dictionary that maps cluster identifier number to their sum frequency
    cluster_frequency_dict = cluster_matrix.groupby('cluster').sum()['metabolite frequency'].to_dict()
    cluster_matrix['cluster frequency'] = [cluster_frequency_dict.get(i) for i in cluster_matrix.cluster]

    # 4- Remove the columns with distance, keeping only cluster and metabolite frequencies
    selection_matrix = cluster_matrix.iloc[:, -3:]

    # 5- Make a dictionary of max value and cluster identifier
    # Get the maximum metabolite frequency for each cluster
    grouped = selection_matrix.groupby('cluster').max()
    # max_metab_freq = [i for i in grouped['metabolite frequency']]
    cluster_frequency = [i for i in grouped['cluster frequency']]
    clusters = [i for i in grouped.index]
    # max_metab_dict = dict(zip(clusters,max_metab_freq))
    cluster_freq_dict = dict(zip(clusters, cluster_frequency))

    # Normalize cluster frequency
    cluster_freq = pd.DataFrame(cluster_freq_dict.items(), columns=['Cluster', 'Frequency'])

    from scipy.stats import zscore
    weight = []
    for z in zscore(cluster_freq['Frequency']):
        if z > 1:
            weight.append(3)
        elif 1 > z > 0:
            weight.append(2)
        else:
            weight.append(1)

    cluster_weight_dict = dict(zip([c for c in cluster_freq['Cluster']], weight))
    # 6- Add w metabolites from each cluster where w is the weight
    # determined based on the z-score of cluster frequency
    selected_metab = []
    for i, row in selection_matrix.iterrows():
        cluster = row['cluster']
        number_of_metab = cluster_weight_dict.get(cluster)

        for i in selection_matrix[selection_matrix['cluster'] == cluster].sort_values('metabolite frequency',
                                                                                      ascending=False)[
                 :number_of_metab].index:
            selected_metab.append(i)

    return set(selected_metab)

def cluster_metabolites(outpath,
                        path_to_model,
                        CONNECTIVITY_THRESHOLD,
                        THRESHOLD,
                        eps,
                        show_frequency,
                        show_matrix,
                        frequency_fig_name,
                        matrix_fig_name):
    """

    :param outpath:
    :param path_to_model:
    :param CONNECTIVITY_THRESHOLD:
    :param show_frequency:
    :return:
    """
    #Get inputs
    hofs = _get_hof_data(outpath)
    model = _import_model(path_to_model)

    #Generate distance matrix
    distance_matrix = dijkstra.generate_distance_matrix(model,CONNECTIVITY_THRESHOLD)

    #Analyze
    #1- Get the best biomass objective functions accross all evolutions
    best_bof = _filter_hofs(hofs,THRESHOLD)

    #2- Get the frequency of each metabolite
    freq_df  = _make_freq_df(best_bof)

    #3- Cluster the metabolites
    result_matrix = _generate_result_matrix(freq_df,distance_matrix)
    grouped_clusters, cluster_dictionary = _dbscan_clustering(result_matrix,eps)

    #4- Display
    if show_frequency:
        _display_occurence(freq_df,frequency_fig_name)

    if show_matrix:
        _display_result_matrix(result_matrix,
                               cluster_dictionary,
                               freq_df,
                               matrix_fig_name)
    #Select metabolites based on frequency of apparition in hall of fame
    list_of_metab = _select_metabolites(freq_df,grouped_clusters,best_bof)
    return list_of_metab