#!/usr/bin/python3

import subprocess
import pandas as pd



def mash_sketch_1(genome_i):
    "GENOMA -> SKETCH | Creates a mash sketch of the genome"
    output_i = genome_i.replace('.fna', '')
    subprocess.call(['mash', "sketch", genome_i, "-p", "24", "-o", output_i])

def mash_i_vs_noasig_to_create_noasigcsv_1(genome_i):
    "GENOME -> NOASIG_DATA | creates dataframe of noasig (genome vs allnoasig)"
    mash_sketch_1(genome_i)
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    with open('noasig_1.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", input_i, "noasig_1.msh", "-p", "24"], stdout=outfile)
    noasig_data = pd.read_csv("noasig_1.csv", sep='\t', header=None)
    noasig_data.columns = ['Source', 'Target', 'value', 'X1', 'X2']
    return noasig_data

def mash_paste_i_noasig_1(genome_i):
    "GENOMA -> NOASIG.MSH | mash paste of genome to noasig"
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['mash', "paste", "tmp", input_i, "noasig_1.msh"])
    subprocess.call(['mv', 'tmp.msh', "noasig_1.msh"])


def create_noasig_msh_1(genome_i):
    "GENOMA -> NOASIG.MSH | creates the file noasig.msh"

    mash_sketch_1(genome_i)
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['cp', input_i, 'noasig_1.msh'])

def create_new_asigdatacsv_1(genome_i):
    noasig_data = mash_i_vs_noasig_to_create_noasigcsv_1(genome_i)
    mash_paste_i_noasig_1(genome_i)
    noasig_data.to_csv('noasig_1.csv',index=False)
    return noasig_data
