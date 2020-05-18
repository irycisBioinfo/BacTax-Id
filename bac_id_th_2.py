#!/usr/bin/python3

#from starting_defs_2 import *
from check_asig_2 import *
from check_asig_init_2 import *

import os

def bac_id_th_2(genome_i, percentage_input, clique_input, distance_input):

    print("th2",genome_i)

    # if os.path.isfile('noasig_2.msh') == False:
    #     create_noasig_msh_2(genome_i)
    #
    # elif os.path.isfile('noasig_2.csv') == False:
    #     create_new_asigdatacsv_2(genome_i)
    #
    # else:

    if os.path.isfile('tabla_asig.csv') == True:
        print("AAAAAAAAAAAAAA")
        tabla_asig = pd.read_csv("tabla_asig.csv")
        if genome_i in tabla_asig["Genome"].values:
            #si la columna de tabla_asig esta vacia
            if tabla_asig.T2.isnull().all():
            #if os.path.isfile('tabla_asig_2.csv') == False:
                check_asig_init_2(genome_i,clique_input, distance_input, percentage_input)
            else:
                check_asig_2(genome_i, clique_input, distance_input, percentage_input)
