#!/usr/bin/python3

from check_asig_defs_1 import *

import pandas as pd

def check_asig_1(genome_i, percentage, clique, distance):
    noasig_data = pd.read_csv("noasig_1.csv")
    tabla_asig = pd.read_csv("tabla_asig.csv")

    asig_data = mash_i_vs_asig_to_create_new_asigdatacsv_1(genome_i)
    #print(asig_data['Target'])

    # chequeo si la distancia más baja es menor al threshold marcado para el nivel
    if min(asig_data['value']) < distance:
        print("A")
        clique_list, satelite_list, lowest_value_clique = list_of_connections_i_vs_possible_pseudoclique_1(asig_data, tabla_asig, distance)

        # comparo la longitud de la lista de vecinos con la del clique
        if len(satelite_list) > percentage * len(clique_list):
            print("B")
            mash_paste_i_asig_1(genome_i)
            tabla_asig = add_satellite_to_tabla_asig_1(genome_i, lowest_value_clique, tabla_asig)

        else:
            print("C")
            noasig_data = mash_i_vs_noasig_to_create_uniqcsv_and_new_noasig_data_1(genome_i, noasig_data)
            #noasig_data.to_csv('noasig_1.csv', index=False)
            # buscar clique
            if min(noasig_data['value']) < distance:
                print("D")
                largest_clique, clique_data = from_noasig_to_largest_clique_1(noasig_data, distance)

                if len(largest_clique) > clique:
                    print("E")
                    #!!!!!!!!!!!!
                    noasig_data, clique_and_satellites = search_for_satellites_from_new_clique_1(clique_data, largest_clique, noasig_data, percentage)

                    tabla_asig = append_new_clique_to_tabla_asig_1(tabla_asig, clique_and_satellites)

                    for genome_asig in largest_clique:
                        mash_paste_i_asig_1(genome_asig)




    else:
        print("F")

        noasig_data = mash_i_vs_noasig_to_create_uniqcsv_and_new_noasig_data_1(genome_i, noasig_data)

        # buscar clique
        if min(noasig_data['value']) < distance:
            print("G")
            largest_clique, clique_data = from_noasig_to_largest_clique_1(noasig_data, distance)

            if len(largest_clique) > clique:
                print("H")
                # !!!!!!!!!!!!
                noasig_data, clique_and_satellites = search_for_satellites_from_new_clique_1(clique_data, largest_clique,
                                                                                           noasig_data, percentage)
                tabla_asig = append_new_clique_to_tabla_asig_1(tabla_asig, clique_and_satellites)

                #no es elegante hacer por un lado el clique y por otro los satelites del paste
                for genome_asig in largest_clique:
                    mash_paste_i_asig_1(genome_asig)

    tabla_asig.to_csv('tabla_asig.csv', index=False)
    noasig_data.to_csv('noasig_1.csv', index=False)
