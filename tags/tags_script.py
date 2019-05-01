
import argparse
import collections
import sys
import os


parser = argparse.ArgumentParser(description="LaTeX tagger - Python program for tags on LaTeX source code.")

parser.add_argument(
    '-id', '--input_dict',
    help='Specify the input dictionary(ies)',
    required=False,
    nargs='+'
    )

parser.add_argument(
    '-it', '--input_tex',
    help='Specify the input source LaTeX file(s)',
    required=False,
    nargs='+'
    )

parser.add_argument(
    '-od', '--output_dict',
    help='Specify to output a dictionary file. If a name is specified, it will be used as the file name; if not, the default \'tags.dict\' file name will be used.)',
    required=False,
    nargs='?',
    default=False,
    action='store'
    )

parser.add_argument(
    '-ot', '--output_tex',
    help='Specify to output a LaTeX file. If a name is specified, it will be used as the file name; if not, the default \'tags.tex\' file name will be used.)',
    required=False,
    nargs='?',
    default=False,
    action='store'
    )

argsdict=vars(parser.parse_args())

if argsdict['input_tex'] == None and argsdict['input_dict'] == None:
    print('Please, specify a tex or a dict file.')
    sys.exit(0)

if argsdict['output_tex'] == False and argsdict['output_dict'] == False:
    print('Please, specify an output file type (using the default or a custom name).')
    sys.exit(0)

if argsdict['output_dict'] == None:
    argsdict['output_dict'] = 'tags.dict'

if argsdict['output_tex'] == None:
    argsdict['output_tex'] = 'tags.tex'


#################################################
# ROUTES
# A: read_tex_tags   -> dict_list                                   => combine_dicts -> gen_dict => write_tex
# B: read_tex_tags   -> dict_list                                   => combine_dicts -> gen_dict => write_dict
# C: read_tex_tags   -> dict_list                                   => combine_dicts -> gen_dict => (write_tex + write_dict)
# D: read_dict_tags  -> dict_list                                   => combine_dicts -> gen_dict => write_tex
# E: read_dict_tags  -> dict_list                                   => combine_dicts -> gen_dict => write_dict
# F: read_dict_tags  -> dict_list                                   => combine_dicts -> gen_dict => (write_tex + write_dict)
# G: ((read_tex_tags -> dict_list) + (read_dict_tags -> dict_list)) => combine_dicts -> gen_dict => write_tex
# H: ((read_tex_tags -> dict_list) + (read_dict_tags -> dict_list)) => combine_dicts -> gen_dict => write_dict
# I: ((read_tex_tags -> dict_list) + (read_dict_tags -> dict_list)) => combine_dicts -> gen_dict => (write_tex + write_dict)
#################################################



def read_tex_tags(argsdict, r): #return part_dict
    # loading files
    file = open(argsdict['input_tex'][r], 'r')
    tags_dict_list  = []; i=-1
    ref = []
    for line in file:
        if line.find('\\label{day:') != -1:
            loca = line.find('\\label{')
            locb = line.find('}\n')
            ref.append('\\ref{' + str(line[loca+7:locb]) + '}')
            i+=1;
            tags_dict_list.append(dict())
            
        if line.find('\\tags{') != -1:
            locA = line.find('\\tags{')
            locB = line.find('}')
            tags_list = line[locA+6:locB]
            #print(tags_list)
            tags_list = tags_list.split(', ')
            #print(tags_dict_list)
            for j in range(len(tags_list)):
                tags_dict_list[i][tags_list[j]] = ref

    #Generation of a list with the different tags from the tex file
    tags_list_complete = []
    for i in range(len(tags_dict_list)):
        for j in range(len(list(tags_dict_list[i].keys()))):
            tags_list_complete.append(list(tags_dict_list[i].keys())[j])
    tags_list_complete_set = list(set(tags_list_complete))

    #Generation of the dictionary
    part_dict = dict()
    for i in range(len(tags_list_complete_set)):
        for j in range(len(tags_dict_list)):
            
            if tags_list_complete_set[i] in list(tags_dict_list[j].keys()):
                
                if tags_list_complete_set[i] in list(part_dict.keys()):
                
                    if len(part_dict[tags_list_complete_set[i]]) == 0:
                        part_dict[tags_list_complete_set[i]] = [ref[j]]
                    elif len(part_dict[tags_list_complete_set[i]]) != 0: 
                        part_dict[tags_list_complete_set[i]].append(ref[j])
                
                elif tags_list_complete_set[i] not in part_dict.keys():
                    part_dict[tags_list_complete_set[i]] = [ref[j]]
    
    part_dict = dict(collections.OrderedDict(sorted(part_dict.items())))

    return part_dict

def read_dict_tags(argsdict, s): #return part_dict
    file = open(argsdict['input_dict'][s], 'r').read()
    part_dict = eval(file)

    return part_dict


# CREATE COMBINE_DICTS AGAIN!!!!!!
def combine_dicts(dict_list): #return gen_dict
    tags_list_complete = []
    for i in range(len(dict_list)):
        for j in range(len(list(dict_list[i].keys()))):
            tags_list_complete.append(list(dict_list[i].keys())[j])
    tags_list_complete_set = list(set(tags_list_complete))

    gen_dict = dict()
    for i in range(len(tags_list_complete_set)):
        for j in range(len(dict_list)):
            if dict_list[j][str(tags_list_complete_set[i])]:
                if len(dict_list[j][str(tags_list_complete_set[i])]) != 0:
                    gen_dict[str(tags_list_complete_set[i])].append(dict_list[j][str(tags_list_complete_set[i])])
                    
                else :
                    gen_dict[str(tags_list_complete_set[i])] = []
                    gen_dict[str(tags_list_complete_set[i])].append(dict_list[j][str(tags_list_complete_set[i])])

    
    gen_dict = dict(collections.OrderedDict(sorted(gen_dict.items())))

    for k in range(len(list(gen_dict.keys()))):
        gen_dict[str(list(gen_dict.keys())[k])] = list(set(gen_dict[str(list(gen_dict.keys())[k])]))

    return gen_dict

def write_dict(argsdict,gen_dict):
    f = open(argsdict['output_dict'], 'w')

    f.write(str(gen_dict))

    f.close()
    #returns 0

def write_tex(argsdict, gen_dict):
    # printing on file
    out = open(argsdict['output_tex'], 'w')

    out.write('\\chapter{Table of Tags}\n\n')
    out.write('\\begin{itemize}\n')

    for key in gen_dict:
        line_ = '\t\\item \\textbf{' + str(key) + ':}'
        for i in range (len(gen_dict[key])):
            if i == 0:
                line_ = line_ + ' Section(s) ' + str(gen_dict[key][i])
            else :
                line_ = line_ + ', ' + str(gen_dict[key][i])
        out.write(str(line_) + '\n')
    out.write('\\end{itemize}\n')
    out.close()

    #return 0


dict_list = []

if argsdict['input_dict'] != None:
    
    for s in range(len(argsdict['input_dict'])):
        dict_list.append(read_dict_tags(argsdict, s))


if argsdict['input_tex'] != None:
    
    for r in range(len(argsdict['input_tex'])):
        dict_list.append(read_tex_tags(argsdict, r))


gen_dict = combine_dicts(dict_list)

if argsdict['output_dict'] != False:
     write_dict(argsdict,gen_dict)

if argsdict['output_tex'] != False:
     write_tex(argsdict,gen_dict)
     