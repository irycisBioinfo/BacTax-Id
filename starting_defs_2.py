#!/usr/bin/python3

import subprocess
import pandas as pd

def mash_i_vs_noasig_to_create_noasigcsv_2(genome_i):
    "GENOME -> NOASIG_DATA | creates dataframe of noasig (genome vs allnoasig)"
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    with open('noasig_2.csv', 'w') as outfile:
        subprocess.call(['mash', "dist", input_i, "noasig_2.msh", "-p", "24"], stdout=outfile)
    noasig_data = pd.read_csv("noasig_2.csv", sep='\t', header=None)
    noasig_data.columns = ['Source', 'Target', 'value', 'X1', 'X2']
    return noasig_data

def mash_paste_i_noasig_2(genome_i):
    "GENOMA -> NOASIG.MSH | mash paste of genome to noasig"
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['mash', "paste", "tmp", input_i, "noasig_2.msh"])
    subprocess.call(['mv', 'tmp.msh', "noasig_2.msh"])


def create_noasig_msh_2(genome_i):
    "GENOMA -> NOASIG.MSH | creates the file noasig.msh"
    output_i = genome_i.replace('.fna', '')
    input_i = output_i + ".msh"
    subprocess.call(['cp', input_i, 'noasig_2.msh'])

def create_new_asigdatacsv_2(genome_i):
    noasig_data = mash_i_vs_noasig_to_create_noasigcsv_2(genome_i)
    mash_paste_i_noasig_2(genome_i)
    noasig_data.to_csv('noasig_2.csv',index=False)
    return noasig_data
