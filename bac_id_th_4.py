#!/usr/bin/python3

from check_asig_4 import *
from check_asig_init_4 import *

import os

def bac_id_th_4(genome_i, percentage_input, clique_input, distance_input):

    print("th4",genome_i)

    if os.path.isfile('tabla_asig.csv') == True:
        tabla_asig = pd.read_csv("tabla_asig.csv")

        if tabla_asig.T3.isnull().all() == False:
            #tabla_asig = pd.read_csv("tabla_asig.csv")

            if genome_i in tabla_asig["Genome"].values:
                fila = tabla_asig.loc[tabla_asig['Genome'] == genome_i]
                if pd.isnull(fila['T3'].iloc[0]) == False:
                #tengo que chequear si genome_i tiene valor en T3

                    #si la columna de tabla_asig esta vacia
                    if tabla_asig.T4.isnull().all():
                    #if os.path.isfile('tabla_asig_4.csv') == False:
                        check_asig_init_4(genome_i,clique_input, distance_input, percentage_input)
                    else:
                        check_asig_4(genome_i, clique_input, distance_input, percentage_input)