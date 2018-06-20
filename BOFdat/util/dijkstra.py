import BOFdat as bd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

from collections import defaultdict
class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list)
        self.distances = {}

    def add_node(self, value):
        self.nodes.add(value)

    def add_edge(self, from_node, to_node, distance):
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.distances[(from_node, to_node)] = distance

def _dijsktra(graph, initial):
    visited = {initial: 0}
    path = {}

    nodes = set(graph.nodes)

    while nodes: 
        min_node = None
        for node in nodes:

            if node in visited:
                if min_node is None:
                    min_node = node
                elif visited[node] < visited[min_node]:
                    min_node = node

        if min_node is None:
            break

        nodes.remove(min_node)
        current_weight = visited[min_node]

        for edge in graph.edges[min_node]:
            weight = current_weight + graph.distances[(min_node, edge)]
            if edge not in visited or weight < visited[edge]:
                visited[edge] = weight
                path[edge] = min_node

    return visited, path

def generate_distance_matrix(model,THRESHOLD=15):
    #List of co-enzymes (highly connected metabolites)
    co_enzymes = []
    for m in model.metabolites:
        if len(m.reactions) > THRESHOLD:
            co_enzymes.append(m.id)

    metab_network = Graph()
    #Add nodes
    for m in model.metabolites:
        if m.id not in co_enzymes:
            metab_network.add_node(str(m.id))
    #Add edges
    for r in model.reactions:
        for reac in r.reactants:
            if reac.id not in co_enzymes:
                for prod in r.products:
                    if prod.id not in co_enzymes:
                        metab_network.add_edge(str(reac.id),str(prod.id),1)
        for prod in r.products:
            if  prod.id not in co_enzymes:
                for reac in r.reactants:
                    if reac.id not in co_enzymes:
                        metab_network.add_edge(str(prod.id),str(reac.id),1)

    dijsktra_matrix = []
    for m in metab_network.nodes:
        dijsktra_matrix.append(_dijsktra(graph=metab_network,initial=m))

    new_matrix = []
    for i in range(len(metab_network.nodes)):
        # Order the dijsktra matrix for this metabolite
        row = []
        for m in metab_network.nodes:
            row.append(dijsktra_matrix[i][0].get(m))

        new_matrix.append(row)

    dist_matrix = pd.DataFrame(new_matrix, columns=metab_network.nodes, index=metab_network.nodes)
    return dist_matrix
