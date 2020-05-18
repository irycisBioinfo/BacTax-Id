#!/usr/bin/python3


import subprocess
import pandas as pd
import networkx as nx
import numpy as np

def mash_i_vs_asig_to_create_filtered_new_asigdatacsv_4(genome_i):
    "GENOME -> ASIG_DATA | compare genome to asig_data"

    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    with open('asig_data_4.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", input_i, "asig_4.msh", "-p", "24"], stdout=outfile)

    asig_data = pd.read_csv("asig_data_4.csv", sep='\t', header=None)
    asig_data.columns = ['Source', 'Target', 'value', 'X1', 'X2']

    tabla_asig = pd.read_csv("tabla_asig.csv")
    value_prev_clique_T1 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 1]
    value_prev_clique_T2 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 2]
    value_prev_clique_T3 = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 3]
    #value_prev_clique = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 3]
    clique_prev = tabla_asig.loc[(tabla_asig['T1'] == value_prev_clique_T1) & (tabla_asig['T2'] == value_prev_clique_T2) & (tabla_asig['T3'] == value_prev_clique_T3)]
    #clique_prev = tabla_asig.loc[(tabla_asig['T3'] == value_prev_clique)]
    clique_prev_list = clique_prev['Genome'].tolist()

    #distance dataframe of clique T1 elements
    asig_1 = asig_data[asig_data['Target'].isin(clique_prev_list)]
    asig_2 = asig_1[asig_1['Source'].isin(clique_prev_list)]
    asig_4 = asig_2.drop_duplicates()

    return asig_4, clique_prev_list, value_prev_clique_T3, clique_prev, tabla_asig

def list_of_connections_i_vs_possible_pseudoclique_4(asig_dist_4, distance_input, clique_prev_list, tabla_asig):

    lowest_value_genome_4 = asig_dist_4.loc[asig_dist_4['value'] == min(asig_dist_4['value'])].values[0, 1]
    lowest_value_clique_4 = tabla_asig.loc[tabla_asig['Genome'] == lowest_value_genome_4].values[0, 4]
    clique_4 = tabla_asig.loc[(tabla_asig['T3'] == lowest_value_clique_4)]
    clique_list_4 = clique_4['Genome'].tolist()


    # tabla_asig = pd.read_csv("tabla_asig.csv")
    # value_prev_clique = tabla_asig.loc[tabla_asig['Genome'] == genome_i].values[0, 1]
    # clique_prev = tabla_asig.loc[(tabla_asig['T1'] == value_prev_clique)]
    # clique_prev_list = clique_prev['Genome'].tolist()

    # distance dataframe of clique T1 elements
    # asig_1 = asig_dist_4[asig_dist_4['Target'].isin(clique_prev_list)]
    # asig_4 = asig_1[asig_1['Source'].isin(clique_prev_list)]
    # asig_4 = asig_4.drop_duplicates()

    satelite_posibles = asig_dist_4[asig_dist_4['value'] < distance_input].dropna()
    satelite_posibles = satelite_posibles[satelite_posibles['Target'].isin(clique_list_4)]
    satelite_posibles = satelite_posibles.drop_duplicates()
    satelite_list = satelite_posibles['Target'].tolist()

    return  satelite_list, clique_list_4, lowest_value_clique_4

def mash_paste_i_asig_4(genome_i):
    "GENOME -> ASIG.MSH | Create file  "
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['mash', "paste", "tmp4", input_i, "asig_4.msh"])
    subprocess.call(['mv', 'tmp4.msh', "asig_4.msh"])

def add_satellite_to_tabla_asig_4(genome_i, lowest_value_clique_input, tabla_asig_input):

    tabla_asig_input.loc[tabla_asig_input.Genome == genome_i, 'T4'] = lowest_value_clique_input

    return tabla_asig_input

def create_non_duplicated_asig_dist_matrix_4():

    with open('asig_dist_matrix_4.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", "asig_1.msh", "asig_1.msh", "-p", "24"], stdout=outfile)

    df = pd.read_csv("asig_dist_matrix_4.csv", sep='\t', header=None)
    df.columns = ['Source', 'Target', 'value', 'X1', 'X2']
    df1 = pd.DataFrame(np.sort(df[['Source', 'Target']], axis=1))
    df2 = df[~df1.duplicated()]

    return df2

def from_clique_dist_to_largest_clique_4(asig_dist_matrix_4,clique_prev_list, distance4, clique_prev):


    noasig_4 = clique_prev[clique_prev['T4'].isnull()]
    clique_prev_list_4 = noasig_4['Genome'].tolist()


    #distance dataframe of clique T1 elements and T4 nan
    asig_1 = asig_dist_matrix_4[asig_dist_matrix_4['Target'].isin(clique_prev_list_4)]
    asig_2 = asig_1[asig_1['Source'].isin(clique_prev_list_4)]
    asig_4 = asig_2.drop_duplicates()

    clique_dataframe = asig_4[asig_4['value'] < distance4].dropna()
    G = nx.from_pandas_edgelist(clique_dataframe, 'Source', 'Target')
    largest_clique = max(nx.enumerate_all_cliques(G), key=len)

    return largest_clique, clique_dataframe

def add_first_clique_to_asig_data_and_asigmsh_4(tabla_asig_input, largest_clique_input, clique_prev):
    "CLIQUE -> ASIG_DATA & ASIG.MSH | the first clique asign to asig_data and create of asig.msh"

    tabla_asig_max_value = (clique_prev.T4.max())
    if tabla_asig_max_value > 0:
        tabla_asig_max_value_new = tabla_asig_max_value + 1
    else:
        tabla_asig_max_value_new = 1

    count1 = 0

    for genome_asig in largest_clique_input:
        tabla_asig_input.loc[tabla_asig_input.Genome == genome_asig, 'T4'] = tabla_asig_max_value_new

        count1 += 1
        if count1 == 1:
            input_i = genome_asig.replace('.fna', '.msh')
            subprocess.call(['cp', input_i, "asig_4.msh"])

        input_ii = genome_asig.replace('.fna', '.msh')
        subprocess.call(['mash', "paste", "tmp4", input_ii, "asig_4.msh"])
        subprocess.call(['mv', 'tmp4.msh', "asig_4.msh"])

    return tabla_asig_input, tabla_asig_max_value_new


def search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage):

    #tengo qe filtrar clique_dataframe para quitar los ya asignados

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

def append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig_input, tabla_asig_max_value_new):

    for satellite_asig in satellites:
        tabla_asig_input.loc[tabla_asig_input.Genome == satellite_asig, 'T4'] = tabla_asig_max_value_new

    return tabla_asig_input