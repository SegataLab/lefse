#!/usr/bin/env python3

import sys
import os
import argparse

def read_params(args):
    parser = argparse.ArgumentParser(description='Convert LEfSe output to '
                        'Circlader input')
    parser.add_argument(    'inp_f', metavar='INPUT_FILE', nargs='?', 
                            default=None, type=str, 
                            help="the input file [stdin if not present]")    
    parser.add_argument(    'out_f', metavar='OUTPUT_FILE', nargs='?', 
                            default=None, type=str, 
                            help="the output file [stdout if not present]")
    parser.add_argument('-l', metavar='levels with label', default=0, type=int)

    return vars(parser.parse_args()) 

def lefse2circlader():
    par = read_params(sys.argv)
    finp,fout = bool(par['inp_f']), bool(par['out_f'])

    with open(par['inp_f']) if finp else sys.stdin as inpf:
        put_bm = (l.strip().split('\t') for l in inpf.readlines()) 
    biomarkers = [p for p in put_bm if len(p) > 2]

    circ = [    [   b[0],
                    "" if b[0].count('.') > par['l'] else b[0].split('.')[-1],
                    b[2],
                    b[2]+"_col" ] for b in biomarkers]

    with open(par['out_f'],'w') if fout else sys.stdout as out_file:
        for c in circ:
            out_file.write( "\t".join( c ) + "\n" )

if __name__ == '__main__':
    lefse2circlader()



