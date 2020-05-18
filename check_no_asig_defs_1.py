#!/usr/bin/python3

import subprocess
import pandas as pd
import networkx as nx


def mash_sketch_1(genome_i):
    "GENOMA -> SKETCH | Creates a mash sketch of the genome"
    output_i = genome_i.replace('.fna', '')
    subprocess.call(['mash', "sketch", genome_i, "-p", "24", "-o", output_i])

def mash_paste_i_noasig_1(genome_i):
    "GENOMA -> noasig_1.msh | mash paste of genome to noasig"
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['mash', "paste", "tmp", input_i, "noasig_1.msh"])
    subprocess.call(['mv', 'tmp.msh', "noasig_1.msh"])

def mash_i_vs_noasig_to_create_uniqcsv_and_new_noasig_data_1(genome_i ,noasig_data_input):
    "GENOME -> NEW NOASIG DATA | mash dist of i vs noasig and append to old noasig_data"
    mash_sketch_1(genome_i)
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    with open('uniq_1.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", input_i, "noasig_1.msh", "-p", "24"], stdout=outfile)
    mash_paste_i_noasig_1(genome_i)

    uniq_data = pd.read_csv("uniq_1.csv", sep='\t', header=None)
    uniq_data.columns = ['Source', 'Target', 'value', 'X1', 'X2']

    noasig_data1 = pd.concat([uniq_data, noasig_data_input], ignore_index=True, sort=True)
    #noasig_data1 = noasig_data_input.append(uniq_data)
    return  noasig_data1

def from_noasig_to_largest_clique_1(noasig_data_input, distance_input):
    "NOASIG_DATA -> LARGEST_CLIQUE | FINDS THE LARGEST CLIQUE IN NOSASIG DATA"
    noasig_data = noasig_data_input
    clique_data = noasig_data[noasig_data['value'] < distance_input].dropna()
    G = nx.from_pandas_edgelist(clique_data, 'Source', 'Target')
    largest_clique = max(nx.enumerate_all_cliques(G), key=len)
    return largest_clique, clique_data

def first_clique_to_asig_data_and_asigmsh_1(largest_clique_input, asig_data_input):
    "CLIQUE -> ASIG_DATA & ASIG.MSH | the first clique asign to asig_data and create of asig.msh"

    asig_data_input['Genome'] = largest_clique_input
    asig_data_input['T1'] = 1

    count1 = 0
    for genome_asig in largest_clique_input:
        count1 += 1
        if count1 == 1:

            input_i = genome_asig.replace('.fna', '.msh')
            subprocess.call(['cp', input_i, "asig_1.msh"])
        else:
            input_ii = genome_asig.replace('.fna', '.msh')
            subprocess.call(['mash', "paste", "tmp2", input_ii, "asig_1.msh"])
            subprocess.call(['mv', 'tmp2.msh', "asig_1.msh"])
    return asig_data_input

def search_for_satellites_from_new_clique_1(clique_data_input, largest_clique_input, noasig_data_input, percentage_input):
    "CLIQUE_DATA+LARGEST_CLIQUE -> noasig_data | find new satellites of clique"
    #clique_data = clique_data_input
    # esto es para coger las distacias de los del clique con respecto a todos para ver si hay satelites

    noasig_data1 = clique_data_input[clique_data_input['Target'].isin(largest_clique_input)]
    noasig_data2 = clique_data_input[clique_data_input['Source'].isin(largest_clique_input)]
    noasig_data3 = noasig_data1.append(noasig_data2)
    noasig_data4 = noasig_data3.drop_duplicates()

    H = nx.from_pandas_edgelist(noasig_data4, 'Source', 'Target')
    a = list(nx.nodes(H))
    # coger los satelites con respecto al clique
    b = set(a) - set(largest_clique_input)
    satelites = []
    for x in b:
        g = (list(set(largest_clique_input).intersection(list(nx.all_neighbors(H, x)))))
        if len(g) > percentage_input * len(largest_clique_input):
            entrada7 = x.replace('.fna', '.msh')
            subprocess.call(['mash', "paste", "tmp3", entrada7, "asig_1.msh"])
            subprocess.call(['mv', 'tmp3.msh', "asig_1.msh"])
            satelites.append(x)

        else:
            print("satelite fallido")

    genomes_to_erase = largest_clique_input + satelites

    # eliminar los del clique de noasigdata
    noasig_data1 = noasig_data_input[~noasig_data_input['Target'].isin(genomes_to_erase)]
    noasig_data1 = noasig_data1[~noasig_data1['Source'].isin(genomes_to_erase)]
    noasig_data1 = noasig_data1.drop_duplicates()

    new_noasig_source = noasig_data1['Source'].tolist()
    new_noasig_target = noasig_data1['Target'].tolist()
    new_noasig = new_noasig_source + new_noasig_target

    new_noasig = pd.DataFrame(new_noasig)
    new_noasig = new_noasig.dropna()
    new_noasig = new_noasig.drop_duplicates()

    #creo que esto lo tengo que separar en dos xk hay dos posibles return
    if len(new_noasig) > 0:

        new_noasig_list = new_noasig[0].tolist()
        count1 = 0
        for genome1 in new_noasig_list:
            count1 += 1
            if count1 == 1:
                entrada2 = genome1.replace('.fna', '.msh')
                subprocess.call(['cp', entrada2, "noasig_1.msh"])
            else:
                entrada3 = genome1.replace('.fna', '.msh')
                subprocess.call(['mash', "paste", "tmp33", entrada3, "noasig_1.msh"])
                subprocess.call(['mv', 'tmp33.msh', "noasig_1.msh"])
    else:
        subprocess.call('rm noasig_1.msh', shell=True)
        #esto falta en el original
        subprocess.call('rm noasig_1.csv', shell=True)

    return noasig_data1, satelites

def append_satellites_to_new_clique_in_tabla_asig_1(tabla_asig_input, satellites):

    tabla_asig_max_value = (tabla_asig_input.T1.max())
    #tabla_asig_max_value_new = tabla_asig_max_value + 1
    len_tmp = len(satellites)
    numbers = [tabla_asig_max_value]
    numbers = numbers * len_tmp
    # tmp_df = pd.DataFrame(np.zeros((len_tmp,2)))
    tmp_df = pd.DataFrame({"Genome": satellites,
                           "T1": numbers,
                           })
    tabla_asig1 = tabla_asig_input.append(tmp_df)
    return tabla_asig1

