#!/usr/bin/env python3
# coding: utf-8

from sys import argv, exit

#try :
if len(argv) == 1:
    print('Please, specify the input file as the first argument of the command. \nIt is also possible to activate the printing into a file by setting \'print\' as the second argument.')
    exit(0)

elif len(argv) > 1:
    file = argv[1]

    try :
        if argv[2].lower() == 'print':
                print_file = file.split('.')[0] + '_charges.txt'

    except IndexError:
        print_file = None

#except IndexError:
#    pass


def charges_parser(file):
    log = open(file).readlines()

    start = False
    charges = []

    for l in log:
        if start == False:
            if l == 'Summary of Natural Population Analysis:\n':
                start = True

        elif start == True:
            try :
                int(l.split()[0])
                charges.append([int(l.split()[0]), l.split()[1].capitalize(), float(l.split()[2])])

            except (ValueError, IndexError):
                if l.find('* Total *') != -1:
                    break

    return charges


def ask_indexes():
    quest = input('Specify the atoms to extract the charges (or type \'help\' or \'exit\'): ')

    if quest.lower() == 'help':
        helpmsg = '''
        There are three symbols that can be used:
        \t',' --> separates non-continuous atom's indexes
        \t'-' --> separates continuous atom's indexes (all atoms between the indexes will be selected)
        \t';' --> separates two different selections

        Example: 
        \t1, 3-5; 10-12 --> two selections will be created: Atoms 1, 3, 4 and 5 and Atoms 10, 11 and 12
        '''
        print(helpmsg)

        return 'Help', 'Help'

    elif quest.lower() == 'exit':
        print('Bye!')
        exit(1)

    sels = quest.split(';')

    indexes = []

    for sel in sels:
        noncont = sel.split(',')

        indexes_ = []
        for index in noncont:
            if index.find('-') == -1:
                try :
                    indexes_.append(int(index)-1)

                except ValueError:
                    return 'Error', 'Error'

            elif index.find('-') != -1:
                index_ = index.split('-')
                try :                  
                    if int(index_[0]) > int(index_[-1]):
                        ma = int(index_[0])
                        mi = int(index_[-1])

                    elif int(index_[-1]) > int(index_[0]):
                        ma = int(index_[-1])
                        mi = int(index_[0])

                    for i in range(mi, ma+1):
                        indexes_.append(i-1)

                except ValueError:
                    return 'Error', 'Error'

            print(indexes_)

        indexes.append(indexes_)

    return indexes, sels


def calculate_charges(charges_in, indexes):
    charges_out = []
    for indexes_ in indexes:

        charge = 0
        for index in indexes_:
            charge += charges_in[index][2]
            print(charges_in[index][0])

        charges_out.append(charge)

    return charges_out

def print_charges(sels, charges_out, print_file):
    for i in range(len(sels)):

        print(f'For the atoms {str(sels[i]).strip()}, the charge is: {round(charges_out[i], 2)}')

        if print_file != None:
            pf = open(print_file, 'a')
            pf.write(f'For the atoms {str(sels[i]).strip()}, the charge is: {round(charges_out[i], 2)}\n')
            pf.close()
            print(f'Charges printed to {print_file}')

def main(file, print_file):
    charges_in = charges_parser(file)

    while True:
        indexes, sels = ask_indexes()

        if indexes == 'Help':
            continue

        elif indexes == 'Error':
            print('Check the selection and retype it. Type \'help\' for help.')
            continue

        charges_out = calculate_charges(charges_in, indexes)

        print_charges(sels, charges_out, print_file)


if __name__ == '__main__':
    main(file, print_file)
