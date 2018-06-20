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

def _filter_hofs(hofs,BASELINE):
    #Select only the indviduals above baseline in each HOF
    select_hofs = []
    for hof in hofs:
        baseline_index = (hof['Fitness'] > BASELINE).idxmin()
        select_hofs.append(hof.iloc[:baseline_index,])
    #Pool the hall of fame together
    pooled_hof = pd.concat(select_hofs)
    sorted_pooled_hof = pooled_hof.sort_values('Fitness',ascending=False)
    #Take the 30% best individuals
    final_index = int(len(pooled_hof)*0.3)
    metrics_df = sorted_pooled_hof
    # Only the best hall of fames should be selected
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
    #Convert to DataFrame
    freq_df = pd.DataFrame.from_records(frequency_data,
                                        columns=['Metab', 'Frequency'],
                                        index=[m for m in set(best_bof)])
    return freq_df

def _generate_result_matrix(freq_df,dist_matrix):
    #Select the metabolites based on their apparition frequency in the HALL of FAME
    THRESHOLD = freq_df.mean()[0]
    selected_metab = [m for m in freq_df[freq_df.iloc[:,1] > THRESHOLD].iloc[:,0]]
    all_metab = set(selected_metab)
    #Generate the reduced distance matrix used for clustering of metabolic end goals
    best_bof_dist = dist_matrix[dist_matrix.index.isin(all_metab)]
    transposed_bbd = best_bof_dist.T
    result_matrix = transposed_bbd[transposed_bbd.index.isin(all_metab)]
    result_matrix.fillna(1000,inplace=True)
    return result_matrix

def _agglomerative_clustering(result_matrix):
    # Using agglomerative clustering
    k = int(float(len(result_matrix))/4)
    Hclustering = AgglomerativeClustering(n_clusters=k, affinity='euclidean', linkage='average')
    Hclustering.fit(result_matrix)
    clusters = []
    for i in range(len(Hclustering.labels_)):
        clusters.append((Hclustering.labels_[i], result_matrix.index[i]))

    cluster_df = pd.DataFrame.from_records(clusters, columns=['cluster', 'metab'])
    grouped = cluster_df.groupby('cluster')
    grouped_clusters = grouped.agg(lambda x: ', '.join(x))
    return grouped_clusters

def _display_occurence(freq_df,fig_name):
    from pylab import rcParams
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

def find(x,my_palette,grouped_clusters):
    #1- Find the cluster to which x belongs
    def find_cluster(x,grouped_clusters):
        for i in range(len(grouped_clusters)):
            if x in grouped_clusters.loc[i][0].split(', '):
                return grouped_clusters.loc[i].name

    cluster_id = find_cluster(x,grouped_clusters)
    return my_palette.get(cluster_id)



def _display_result_matrix(grouped_clusters,cluster_matrix,figname):
    from pylab import rcParams
    default_param = {'figsize': (15, 15),
                     'sb_style': 'whitegrid',
                     'fontsize': 12}

    my_palette = dict(zip(grouped_clusters.index, sb.color_palette('hls', len(grouped_clusters))))
    my_dark_palette = dict(zip(grouped_clusters.index, sb.hls_palette(len(grouped_clusters), h=0.01, l=0.5, s=0.6)))
    row_colors = [c for c in cluster_matrix.index.map(lambda x: find(x, my_palette, grouped_clusters))]

    rcParams['figure.figsize'] = default_param.get('figsize')
    rcParams['font.size'] = default_param.get('fontsize')
    g = sb.clustermap(cluster_matrix,
                      cmap='coolwarm',
                      method='centroid',
                      figsize=(15, 15),
                      linewidth=0.3,
                      row_colors=row_colors,
                      col_colors=row_colors,
                      xticklabels=False,
                      vmax=15.)

    for tick_label in g.ax_heatmap.get_yticklabels():
        tick_text = tick_label.get_text()
        # Get the color
        tick_label.set_alpha(1.)
        tick_label.set_color(find(tick_text, my_dark_palette, grouped_clusters))

    plt.savefig(figname)
    plt.show()

def _select_metabolites(freq_df,grouped_clusters,best_bof):
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
    print(better_results_df)
    compare_cluster_ga = pd.concat([grouped_clusters, better_results_df], axis=1)
    print(compare_cluster_ga)
    bd3_unique = []
    for i, row in compare_cluster_ga.iterrows():
        mini_unique = [(m, freq_df.loc[m]['Frequency']) for m in row['metab'].split(', ')]
        cluster = pd.DataFrame.from_records(mini_unique, columns=['metab', 'freq'])
        cluster.index = cluster.metab
        del cluster['metab']
        for m in cluster.index[cluster['freq'] == cluster.max()[0]]:
            bd3_unique.append(m)
    print(bd3_unique)
    return bd3_unique


def cluster_metabolites(outpath,
                        path_to_model,
                        CONNECTIVITY_THRESHOLD,
                        BASELINE,
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
    best_bof = _filter_hofs(hofs,BASELINE)
    freq_df  = _make_freq_df(best_bof)
    result_matrix = _generate_result_matrix(freq_df,distance_matrix)
    grouped_clusters = _agglomerative_clustering(result_matrix)
    #Display
    if show_frequency:
        _display_occurence(freq_df,frequency_fig_name)

    if show_matrix:
        _display_result_matrix(grouped_clusters,result_matrix,matrix_fig_name)
    #Select metabolites based on frequency of apparition in hall of fame
    list_of_metab = _select_metabolites(freq_df,grouped_clusters,best_bof)
    return list_of_metab


"""
DEPRECATED
Functions only useful if trying to compare with the original biomass which isnt the case here
def _generate_distance(dist_matrix,original_bof,compare_bof):
    non_shared_metab = [m for m in compare_bof if m not in original_bof]
    non_shared_dist_matrix = dist_matrix[dist_matrix.index.isin(non_shared_metab)]
    transposed = non_shared_dist_matrix.T
    transposed_cut  = transposed[transposed.index.isin(original_bof)]
    distance_result = pd.concat([transposed_cut.idxmin(),transposed_cut.min()],axis=1)
    return distance_result

def _generate_hof_dist_matrix(model, dist_matrix, hofs):
    # Incorporate the data and generate a matrix based on that
    #biomass = _get_biomass_objective_function(model) --> this is useful only to compare
    #original_bof = [str(m.id) for m in biomass.metabolites] --> this is useful only to compare
    data = []
    for i in range(len(hofs)):
        # 1- Make a list from the best individual of each all of fame
        HOF_bof = _make_list([m for m in hofs[i].loc[0:0]['Biomass composition']][0])
        #distance_df = _generate_distance(dist_matrix, original_bof, compare_bof) --> this is useful only to compare
        data.append(distance_metric(HOF_bof, hofs[i].loc[0:0]['Fitness'][0]))

    metrics_df = pd.DataFrame.from_records(data, columns=['Fitness', 'Sum', 'Mean', 'Median'])
    return metrics_df

def _make_list(m):
    #This function makes a list of metabolites from a hall of famer individual
    l = []
    for i in m.split(','):
        i = i.replace("'",'')
        i = i.replace("[",'')
        i = i.replace("]",'')
        i = i.replace(" ",'')
        l.append(i)
    return l

def _filter_hofs(metrics_df,hofs,BASELINE):
    # Only the best hall of fames should be selected
    best_bof = []
    for i in range(len(metrics_df)):
        if metrics_df.iloc[i, 0] > BASELINE:
            l = _make_list([m for m in hofs[i].loc[0:0]['Biomass composition']][0])
            for m in l:
                best_bof.append(m)
    return best_bof

def _make_freq_df(best_bof):
    def occurence(metab, best_bof):
        return float(len([m for m in best_bof if m == metab])) / len(set(best_bof))
    #Generate the frequency of each metabolite in hall of fames
    frequency_data = []
    for m in set(best_bof):
        f_m = occurence(m, best_bof)
        frequency_data.append((m, f_m))
    #Convert to DataFrame
    freq_df = pd.DataFrame.from_records(frequency_data,
                                        columns=['Metab', 'Frequency'],
                                        index=[m for m in set(best_bof)])
    return freq_df

This filter hofs function should be tweaked to account for the modelling
(maybe include that modelling function in the tool...)
"""
