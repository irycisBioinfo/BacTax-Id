#!/usr/bin/python3

from starting_defs_1 import *
from check_no_asig_1 import *
from check_asig_1 import *

import os

def bac_id_th_1(genome_i, percentage_input, clique_input, distance_input):
    print("th1",genome_i)

    if os.path.isfile('noasig_1.msh') == False:
        create_noasig_msh_1(genome_i)

    elif os.path.isfile('noasig_1.csv') == False:
        create_new_asigdatacsv_1(genome_i)

    else:

        if os.path.isfile('tabla_asig.csv') == False:
            print("NOASIG")
            check_no_asig_1(genome_i, percentage_input, clique_input, distance_input)

        else:
            print("ASIG")
            check_asig_1(genome_i, percentage_input, clique_input, distance_input)
