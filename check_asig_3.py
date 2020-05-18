#!/usr/bin/python3

from check_asig_defs_3 import *

def check_asig_3(genome_i,clique, distance3, percentage):

    asig_data_3, clique_prev_list, value_prev_clique,clique_prev, tabla_asig = mash_i_vs_asig_to_create_filtered_new_asigdatacsv_3(genome_i)

    if asig_data_3.empty == False:
        if min(asig_data_3['value']) < distance3:
            satelite_list, clique_list_3, lowest_value_clique_3 = list_of_connections_i_vs_possible_pseudoclique_3(asig_data_3,distance3, clique_prev_list, tabla_asig)
            # comparo la longitud de la lista de vecinos con la del clique
            if len(satelite_list) > percentage * len(clique_list_3):
                mash_paste_i_asig_3(genome_i)
                tabla_asig = add_satellite_to_tabla_asig_1(genome_i, lowest_value_clique_3, tabla_asig)

            else:
                asig_dist_matrix_3 = create_non_duplicated_asig_dist_matrix_3()
                largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_3(asig_dist_matrix_3, clique_prev_list, distance3, clique_prev)
                if len(largest_clique) > clique:
                    tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_3(tabla_asig, largest_clique, clique_prev)
                    satellites = search_for_satellites_from_clique_3(clique_dataframe, largest_clique, percentage)
                    tabla_asig = append_satellites_to_new_clique_tabla_asig_3(satellites, tabla_asig, tabla_asig_max_value_new)

        else:
            asig_dist_matrix_3 = create_non_duplicated_asig_dist_matrix_3()
            largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_3(asig_dist_matrix_3, clique_prev_list, distance3, clique_prev)
            if len(largest_clique) > clique:
                tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_3(tabla_asig,largest_clique,clique_prev)
                satellites = search_for_satellites_from_clique_3(clique_dataframe, largest_clique, percentage)
                tabla_asig = append_satellites_to_new_clique_tabla_asig_3(satellites, tabla_asig, tabla_asig_max_value_new)
    else:
        asig_dist_matrix_3 = create_non_duplicated_asig_dist_matrix_3()
        largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_3(asig_dist_matrix_3, clique_prev_list,
                                                                                distance3, clique_prev)
        if len(largest_clique) > clique:
            tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_3(tabla_asig,
                                                                                               largest_clique,
                                                                                               clique_prev)
            satellites = search_for_satellites_from_clique_3(clique_dataframe, largest_clique, percentage)
            tabla_asig = append_satellites_to_new_clique_tabla_asig_3(satellites, tabla_asig, tabla_asig_max_value_new)

    tabla_asig.to_csv('tabla_asig.csv', index=False)