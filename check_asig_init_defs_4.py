#!/usr/bin/python3

import subprocess
import pandas as pd
import numpy as np
import networkx as nx


def create_non_duplicated_asig_dist_matrix_4():

    with open('asig_dist_matrix_4.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", "asig_1.msh", "asig_1.msh", "-p", "24"], stdout=outfile)

    df = pd.read_csv("asig_dist_matrix_4.csv", sep='\t', header=None)
    df.columns = ['Source', 'Target', 'value', 'X1', 'X2']
    df1 = pd.DataFrame(np.sort(df[['Source', 'Target']], axis=1))
    df2 = df[~df1.duplicated()]

    return df2

def from_clique_dist_to_largest_clique_4(genome_i, asig_dist_matrix_4, distance4):
    #se podria coger el valor en vez de usar read
    tabla_asig = pd.read_csv("tabla_asig.csv")
    value_prev_clique_T1 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 1]
    value_prev_clique_T2 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 2]
    value_prev_clique_T3 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 3]

    #LIADA MAXIMA!!!! al decir que sea igual que el anterior cojo 2.1 y 3.1
    clique_prev = tabla_asig.loc[(tabla_asig['T1'] == value_prev_clique_T1) & (tabla_asig['T2'] == value_prev_clique_T2) & (tabla_asig['T3'] == value_prev_clique_T3)]
    #clique_prev = tabla_asig.loc[(tabla_asig['T3'] == value_prev_clique)]
    clique_prev_list = clique_prev['Genome'].tolist()

    #distance dataframe of clique T1 elements
    asig_1 = asig_dist_matrix_4[asig_dist_matrix_4['Target'].isin(clique_prev_list)]
    asig_2 = asig_1[asig_1['Source'].isin(clique_prev_list)]
    asig_4 = asig_2.drop_duplicates()

    clique_dataframe = asig_4[asig_4['value'] < distance4].dropna()
    G = nx.from_pandas_edgelist(clique_dataframe, 'Source', 'Target')
    largest_clique = max(nx.enumerate_all_cliques(G), key=len)

    return largest_clique, clique_dataframe

def add_first_clique_to_asig_data_and_asigmsh_4(largest_clique_input):
    "CLIQUE -> ASIG_DATA & ASIG.MSH | the first clique asign to asig_data and create of asig.msh"

    tabla_asig = pd.read_csv("tabla_asig.csv")
    #
    # for genome in largest_clique_input:
    #     tabla_asig.loc[tabla_asig.Genome == genome, 'T4'] = 1

    count1 = 0

    for genome_asig in largest_clique_input:
        tabla_asig.loc[tabla_asig.Genome == genome_asig, 'T4'] = 1

        count1 += 1
        if count1 == 1:
            input_i = genome_asig.replace('.fna', '.msh')
            subprocess.call(['cp', input_i, "asig_4.msh"])

        input_ii = genome_asig.replace('.fna', '.msh')
        subprocess.call(['mash', "paste", "tmp4", input_ii, "asig_4.msh"])
        subprocess.call(['mv', 'tmp4.msh', "asig_4.msh"])

    return tabla_asig



def search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage):

    asig_data1 = clique_dataframe[clique_dataframe['Target'].isin(largest_clique)]
    asig_data2 = clique_dataframe[clique_dataframe['Source'].isin(largest_clique)]
    asig_data4 = asig_data1.append(asig_data2)
    asig_data4 = asig_data4.drop_duplicates()

    H = nx.from_pandas_edgelist(asig_data4, 'Source', 'Target')
    a = list(nx.nodes(H))
    b = set(a) - set(largest_clique)
    satellites = []

    for x in b:
        g = (list(set(largest_clique).intersection(list(nx.all_neighbors(H, x)))))
        if len(g) > percentage * len(largest_clique):
            entrada7 = x.replace('.fna', '.msh')
            subprocess.call(['mash', "paste", "tmp4", entrada7, "asig_4.msh"])
            subprocess.call(['mv', 'tmp4.msh', "asig_4.msh"])
            satellites.append(x)

        else:
            print("satelite fallido")

    return satellites

def append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig_input):

    for satellite_asig in satellites:
        tabla_asig_input.loc[tabla_asig_input.Genome == satellite_asig, 'T4'] = 1


    return tabla_asig_input



