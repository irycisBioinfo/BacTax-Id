#!/usr/bin/python3

from check_asig_init_defs_2 import *

def check_asig_init_2(genome_i,clique, distance2, percentage):

    asig_dist_matrix_2 = create_non_duplicated_asig_dist_matrix_2()
    largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_2(genome_i, asig_dist_matrix_2, distance2)

    if len(largest_clique) > clique:

        tabla_asig = add_first_clique_to_asig_data_and_asigmsh_2(largest_clique)
        satellites = search_for_satellites_from_clique_2(clique_dataframe, largest_clique, percentage)
        tabla_asig = append_satellites_to_new_clique_tabla_asig_2(satellites, tabla_asig)
        tabla_asig.to_csv('tabla_asig.csv', index=False)



