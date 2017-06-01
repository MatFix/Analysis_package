#! /usr/bin/env python3
'''Rename the .omf files in a directory to the format A_{ii}.ovf, where 
ii is an increasing counter

Args:
    -r <DIR>    :    renames the .ovf files in <DIR>
    no arg      :    renames all .ovf files in the current directory

Notes:
    If no <DIR> is set it defaults to the current working directory
'''

import os
import re
import sys
import fileinput

def renamer(folder):
    ''' Rename the .ovf files in a directory

    Args:
        <DIR>
    '''
    table = os.listdir(folder)
    table.sort()
    p = re.compile('\.ovf')
    i = 1
    for a in table:
        if p.search(a):
            os.rename(folder+'\\'+a,folder+'\\A_{}.ovf'.format(i))
            i += 1
            


def main():
    renamer(os.getcwd())

if __name__ == '__main__':
    try:
        folder = sys.argv[2]
    except IndexError:
        folder = os.getcwd()
    try:
        if sys.argv[1] == '-r':
            print('Renaming files in selected directory')
            renamer(folder)
        elif sys.argv[1] == '-s':
            try:
                if sys.argv[3].isdigit():
                    n = int(sys.argv[3])
                else:
                    raise(TypeError)
                if n < 0:
                    raise(ValueError)
            except (TypeError,ValueError):
                print('{} is not a positive integer, setting n = 0'.format(sys.argv[3]))
                n = 0
            except IndexError:
                n = 0
        elif sys.argv[1] == '-a':
            print('Processing files in selected directory')
            main(folder)
        else:
            print('Command unrecognized')
    except IndexError:
        print('Processing files in current directory')
        main()

