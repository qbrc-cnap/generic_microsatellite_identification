import os
import json

def map_inputs(user, unmapped_data, id_list):
    '''
    This maps the genome string to the resources needed to run a WDL
    This is important because Cromwell localizes ALL files, so providing
    a Map[String, File] will pull all the files listed therein.

    unmapped_data is a string giving the genome
    id_list is a list of the WDL input names.  The order is given in the gui.json
    file.
    '''

    genome_choice = unmapped_data
    this_directory = os.path.dirname(os.path.abspath(__file__))
    resource_file = os.path.join(this_directory, 'genome_resources.json')
    j = json.load(open(resource_file))

    d = {}
    d[id_list[0]] = genome_choice
    d[id_list[1]] = j[genome_choice]['ref_fasta']
    d[id_list[2]] = j[genome_choice]['ref_fasta_index']
    d[id_list[3]] = j[genome_choice]['ref_index_1']
    d[id_list[4]] = j[genome_choice]['ref_index_2']
    d[id_list[5]] = j[genome_choice]['ref_index_3']
    d[id_list[6]] = j[genome_choice]['ref_index_4']
    d[id_list[7]] = j[genome_choice]['ref_rev_1']
    d[id_list[8]] = j[genome_choice]['ref_rev_2']
    return d
