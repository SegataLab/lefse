#!/usr/bin/env python3

import sys

def read_params(args):
    import argparse as ap
    import textwrap

    p = ap.ArgumentParser( description= '''
    Script will convert QIIME TSV BIOM table for use with lefse.  It is imperative that this table has taxa metadata associated with it named 'Consensus Lineage', this can be down with e.g. the follow biom convert script: ----
    biom convert -i otu.biom -o otu.txt --to-tsv --header-key Taxonomy --output-metadata-id 'Consensus Lineage'
    ''')
    
    p.add_argument( '--in', metavar='INPUT_FILE', type=str, 
                    nargs='?', default=sys.stdin,
                    help=   "the Qiime OTU table file "
                            "[ stdin if not present ]" )
    p.add_argument( '--md', metavar='METADATA_FILE', type=str, 
                    nargs='?', default=None,
                    help=   "the Qiime OTU table file " 
                            "[ only OTU table without metadata if not present ]" )
    p.add_argument( '--out', metavar='OUTPUT_FILE', type=str, 
                    nargs = '?', default=sys.stdout,
                    help=   "the output file "
                            "[stdout if not present]")

    p.add_argument( '-c', metavar="class attribute", 
                    type=str,
                    help =  "the attribute to use as class"   )
    p.add_argument( '-s', metavar="subclass attribute", 
                    type=str,
                    help =  "the attribute to use as subclass"   )
    p.add_argument( '-u', metavar="subject attribute", 
                    type=str,
                    help =  "the attribute to use as subject"   )



    return vars(p.parse_args()) 



def qiime2lefse():
    pars = read_params( sys.argv )
    fin = pars['in']
    fmd = pars['md']
    fout = pars['out']
    all_md = not pars['c'] and not pars['s'] and not pars['u']
    sel_md = [pars['c'],pars['s'],pars['u']]
    with (fin if fin==sys.stdin else open(fin)) as inpf :
        lines = [list(ll) for ll in 
                    (zip(*[l.strip().split('\t') 
                        for l in inpf.readlines()[1:]]) ) ]
    for i,(l1,l2) in enumerate(zip( lines[0], lines[-1] )):
        if not l2 == 'Consensus Lineage':
            lines[-1][i] = l2+"|"+l1

    data = dict([(l[0],l[1:]) for l in lines[1:]])
    
    md = {}
    if fmd:
        with open(fmd) as inpf:
            mdlines = [l.strip().split('\t') for l in inpf.readlines()]
  
        mdf = mdlines[0][1:]

        for l in mdlines:
            mdd = dict(zip(mdf,l[1:]))
            md[l[0]] = mdd

    selected_md = list(md.values())[0].keys() if md else []

    if not all_md:
        selected_md = [s for s in sel_md if s]
    
    out_m = [   selected_md + 
                list([d.replace(";","|").replace("\"","") for d in data[ 'Consensus Lineage' ]])    ]
    for k,v in data.items():
        if k == 'Consensus Lineage':
            continue
        out_m.append( [md[k][kmd] for kmd in selected_md] + list(v) )

    with (fout if fout == sys.stdout else open( fout, "w" )) as outf:
        for l in zip(*out_m):
            outf.write( "\t".join(l) + "\n" )

if __name__ == '__main__':
    qiime2lefse()
