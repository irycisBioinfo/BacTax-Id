#!/usr/bin/python3

from check_asig_init_defs_4 import *

def check_asig_init_4(genome_i,clique, distance4, percentage):

    asig_dist_matrix_4 = create_non_duplicated_asig_dist_matrix_4()
    largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_4(genome_i, asig_dist_matrix_4, distance4)

    if len(largest_clique) > clique:

        tabla_asig = add_first_clique_to_asig_data_and_asigmsh_4(largest_clique)
        satellites = search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage)
        tabla_asig = append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig)
        tabla_asig.to_csv('tabla_asig.csv', index=False)



