#!/usr/bin/python3

from bac_id_th_1 import *
from bac_id_th_2 import *
from bac_id_th_3 import *
from bac_id_th_4 import *



import os
import argparse
import random

parser = argparse.ArgumentParser(description= "¯\_(ツ)_/¯\t\t ")
parser.add_argument('-l','--list', nargs='+', help='<Required> Set flag', required=True, type=str)
parser.add_argument("-p","--percentage",default=0.8, help= 'percentage of similarity')
parser.add_argument("-c","--clique", default=5, help='clique size')
parser.add_argument("-d1", "--distance1", default=0.025, help='threshold distance 1')
parser.add_argument("-d2", "--distance2", default=0.02, help='threshold distance 2')
parser.add_argument("-d3", "--distance3", default=0.01, help='threshold distance 3')
parser.add_argument("-d4", "--distance4", default=0.005, help='threshold distance 4')

args = parser.parse_args()
genomes = args.list
random.shuffle(genomes)

percentage = args.percentage
percentage=float(percentage)

clique = args.clique
clique=int(clique)

distance1 = args.distance1
distance1=float(distance1)
distance2 = args.distance2
distance2=float(distance2)
distance3 = args.distance3
distance3=float(distance3)
distance4 = args.distance4
distance4=float(distance4)

#Folder of the genomes
current_directory = os.getcwd()
os.chdir(current_directory)
print(current_directory)

subprocess.call('rm *.msh', shell= True)
subprocess.call('rm *.csv', shell= True)

count = 0
for genome in genomes:

    print(count, genome)
    count += 1

    bac_id_th_1(genome, percentage, clique, distance1)
    bac_id_th_2(genome, percentage, clique, distance2)
    bac_id_th_3(genome, percentage, clique, distance3)
    bac_id_th_4(genome, percentage, clique, distance4)




