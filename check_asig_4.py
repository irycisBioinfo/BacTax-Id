#!/usr/bin/python3

from check_asig_defs_4 import *

def check_asig_4(genome_i,clique, distance4, percentage):

    asig_data_4, clique_prev_list, value_prev_clique,clique_prev, tabla_asig = mash_i_vs_asig_to_create_filtered_new_asigdatacsv_4(genome_i)

    if asig_data_4.empty == False:
        if min(asig_data_4['value']) < distance4:
            satelite_list, clique_list_4, lowest_value_clique_4 = list_of_connections_i_vs_possible_pseudoclique_4(asig_data_4,distance4, clique_prev_list, tabla_asig)
            # comparo la longitud de la lista de vecinos con la del clique
            if len(satelite_list) > percentage * len(clique_list_4):
                mash_paste_i_asig_4(genome_i)
                tabla_asig = add_satellite_to_tabla_asig_1(genome_i, lowest_value_clique_4, tabla_asig)

            else:
                asig_dist_matrix_4 = create_non_duplicated_asig_dist_matrix_4()
                largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_4(asig_dist_matrix_4, clique_prev_list, distance4, clique_prev)
                if len(largest_clique) > clique:
                    tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_4(tabla_asig, largest_clique, clique_prev)
                    satellites = search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage)
                    tabla_asig = append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig, tabla_asig_max_value_new)

        else:
            asig_dist_matrix_4 = create_non_duplicated_asig_dist_matrix_4()
            largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_4(asig_dist_matrix_4, clique_prev_list, distance4, clique_prev)
            if len(largest_clique) > clique:
                tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_4(tabla_asig,largest_clique,clique_prev)
                satellites = search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage)
                tabla_asig = append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig, tabla_asig_max_value_new)
    else:
        asig_dist_matrix_4 = create_non_duplicated_asig_dist_matrix_4()
        largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_4(asig_dist_matrix_4, clique_prev_list,
                                                                                distance4, clique_prev)
        if len(largest_clique) > clique:
            tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_4(tabla_asig,
                                                                                               largest_clique,
                                                                                               clique_prev)
            satellites = search_for_satellites_from_clique_4(clique_dataframe, largest_clique, percentage)
            tabla_asig = append_satellites_to_new_clique_tabla_asig_4(satellites, tabla_asig, tabla_asig_max_value_new)

    tabla_asig.to_csv('tabla_asig.csv', index=False)