#!/usr/bin/python3

from check_no_asig_defs_1 import *

import pandas as pd

def check_no_asig_1(genome_i, percentage, clique, distance):
    #noasig_data = pd.read_csv("noasig_1.csv", sep='\t', header=None)
    #noasig_data.columns = ['Source', 'Target', 'value', 'X1', 'X2']
    asig_data = pd.DataFrame(columns=['Genome', 'T1','T2','T3','T4'])

    noasig_data = pd.read_csv("noasig_1.csv")
    noasig_data = mash_i_vs_noasig_to_create_uniqcsv_and_new_noasig_data_1(genome_i, noasig_data)

    if min(noasig_data['value']) < distance:
        # clique_data es los datos filtrados por el valor del nivel
        largest_clique, clique_data = from_noasig_to_largest_clique_1(noasig_data, distance)
        if len(largest_clique) > clique:
            asig_data = first_clique_to_asig_data_and_asigmsh_1(largest_clique, asig_data)
            noasig_data, satellites  = search_for_satellites_from_new_clique_1(clique_data, largest_clique, noasig_data, percentage)
            asig_data = append_satellites_to_new_clique_in_tabla_asig_1(asig_data, satellites)
            asig_data.to_csv('tabla_asig.csv', index=False)
            noasig_data.to_csv('noasig_1.csv', index=False)

        noasig_data.to_csv('noasig_1.csv', index=False)

    noasig_data.to_csv('noasig_1.csv',index=False)