#!/usr/bin/python3

from check_asig_defs_2 import *

def check_asig_2(genome_i,clique, distance2, percentage):

    asig_data_2, clique_prev_list, value_prev_clique,clique_prev, tabla_asig = mash_i_vs_asig_to_create_filtered_new_asigdatacsv_2(genome_i)

    if asig_data_2.empty == False:
        if min(asig_data_2['value']) < distance2:
            satelite_list, clique_list_2, lowest_value_clique_2 = list_of_connections_i_vs_possible_pseudoclique_2(asig_data_2,distance2, clique_prev_list, tabla_asig)
            # comparo la longitud de la lista de vecinos con la del clique
            if len(satelite_list) > percentage * len(clique_list_2):
                mash_paste_i_asig_2(genome_i)
                tabla_asig = add_satellite_to_tabla_asig_1(genome_i, lowest_value_clique_2, tabla_asig)

            else:
                asig_dist_matrix_2 = create_non_duplicated_asig_dist_matrix_2()
                largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_2(asig_dist_matrix_2, clique_prev_list, distance2, clique_prev)
                if len(largest_clique) > clique:
                    tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_2(tabla_asig, largest_clique, clique_prev)
                    satellites = search_for_satellites_from_clique_2(clique_dataframe, largest_clique, percentage)
                    tabla_asig = append_satellites_to_new_clique_tabla_asig_2(satellites, tabla_asig, tabla_asig_max_value_new)

        else:
            asig_dist_matrix_2 = create_non_duplicated_asig_dist_matrix_2()
            largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_2(asig_dist_matrix_2, clique_prev_list, distance2, clique_prev)
            if len(largest_clique) > clique:
                tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_2(tabla_asig,largest_clique,clique_prev)
                satellites = search_for_satellites_from_clique_2(clique_dataframe, largest_clique, percentage)
                tabla_asig = append_satellites_to_new_clique_tabla_asig_2(satellites, tabla_asig, tabla_asig_max_value_new)
    else:
        print("ALFIN")
        asig_dist_matrix_2 = create_non_duplicated_asig_dist_matrix_2()
        largest_clique, clique_dataframe = from_clique_dist_to_largest_clique_2(asig_dist_matrix_2, clique_prev_list,
                                                                                distance2, clique_prev)
        if len(largest_clique) > clique:
            tabla_asig, tabla_asig_max_value_new = add_first_clique_to_asig_data_and_asigmsh_2(tabla_asig,
                                                                                               largest_clique,
                                                                                               clique_prev)
            satellites = search_for_satellites_from_clique_2(clique_dataframe, largest_clique, percentage)
            tabla_asig = append_satellites_to_new_clique_tabla_asig_2(satellites, tabla_asig, tabla_asig_max_value_new)

    tabla_asig.to_csv('tabla_asig.csv', index=False)