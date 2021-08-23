#!/usr/bin/env python3

import os,sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
from collections import defaultdict

from lefse.lefse import *
import argparse

colors = ['r','g','b','m','c','y','k','w']

def read_params(args):
    parser = argparse.ArgumentParser(description='Plot results')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="tab delimited input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str, help="the file for the output image")
    parser.add_argument('--feature_font_size', dest="feature_font_size", type=int, default=7, help="the file for the output image")
    parser.add_argument('--format', dest="format", choices=["png","svg","pdf"], default='png', type=str, help="the format for the output file")
    parser.add_argument('--dpi',dest="dpi", type=int, default=72)
    parser.add_argument('--title',dest="title", type=str, default="")
    parser.add_argument('--title_font_size',dest="title_font_size", type=str, default="12")
    parser.add_argument('--class_legend_font_size',dest="class_legend_font_size", type=str, default="10")
    parser.add_argument('--width',dest="width", type=float, default=7.0 )
    parser.add_argument('--height',dest="height", type=float, default=4.0, help="only for vertical histograms")
    parser.add_argument('--left_space',dest="ls", type=float, default=0.2 )
    parser.add_argument('--right_space',dest="rs", type=float, default=0.1 )
    parser.add_argument('--orientation',dest="orientation", type=str, choices=["h","v"], default="h" )
    parser.add_argument('--autoscale',dest="autoscale", type=int, choices=[0,1], default=1 )
    parser.add_argument('--background_color',dest="back_color", type=str, choices=["k","w"], default="w", help="set the color of the background")
    parser.add_argument('--subclades', dest="n_scl", type=int, default=1, help="number of label levels to be dislayed (starting from the leaves, -1 means all the levels, 1 is default )")
    parser.add_argument('--max_feature_len', dest="max_feature_len", type=int, default=60, help="Maximum length of feature strings (def 60)")
    parser.add_argument('--all_feats', dest="all_feats", type=str, default="")
    parser.add_argument('--otu_only', dest="otu_only", default=False, action='store_true', help="Plot only species resolved OTUs (as opposed to all levels)")
    parser.add_argument('--report_features', dest="report_features", default=False, action='store_true', help="Report important features to STDOUT")
    args = parser.parse_args()
    return vars(args)

def read_data(input_file,output_file,otu_only):
    with open(input_file, 'r') as inp:
        if not otu_only:
            rows = [line.strip().split()[:-1] for line in inp.readlines() if len(line.strip().split())>3]
        else:
            rows = [line.strip().split()[:-1] for line in inp.readlines() if len(line.strip().split())>3 and len(line.strip().split()[0].split('.'))==8] # a feature with length 8 will have an OTU id associated with it
    classes = list(set([v[2] for v in rows if len(v)>2]))
    if len(classes) < 1: 
        print("No differentially abundant features found in "+input_file)
        os.system("touch "+output_file)
        sys.exit()
    data = {}
    data['rows'] = rows
    data['cls'] = classes
    return data

def plot_histo_hor(path,params,data,bcl,report_features):
    cls2 = []
    if params['all_feats'] != "":
        cls2 = sorted(params['all_feats'].split(":"))
    cls = sorted(data['cls'])
    if bcl: data['rows'].sort(key=lambda ab: fabs(float(ab[3]))*(cls.index(ab[2])*2-1))
    else: 
        mmax = max([fabs(float(a)) for a in list(zip(*list(data['rows'])))[3]])
        data['rows'].sort(key=lambda ab: fabs(float(ab[3]))/mmax+(cls.index(ab[2])+1))
    pos = arange(len(data['rows']))
    head = 0.75
    tail = 0.5
    ht = head + tail
    ints = max(len(pos)*0.2,1.5)
    fig = plt.figure(figsize=(params['width'], ints + ht), edgecolor=params['back_color'],facecolor=params['back_color'])
    ax = fig.add_subplot(111,frame_on=False,facecolor=params['back_color'])
    ls, rs = params['ls'], 1.0-params['rs']
    plt.subplots_adjust(left=ls,right=rs,top=1-head*(1.0-ints/(ints+ht)), bottom=tail*(1.0-ints/(ints+ht)))

    fig.canvas.manager.set_window_title('LDA results')

    l_align = {'horizontalalignment':'left', 'verticalalignment':'baseline'}
    r_align = {'horizontalalignment':'right', 'verticalalignment':'baseline'}
    added = []
    m = 1 if data['rows'][0][2] == cls[0] else -1
    out_data = defaultdict(list) # keep track of which OTUs result in the plot
    for i,v in enumerate(data['rows']):
        if report_features:
            otu = v[0].split('.')[7].replace('_','.') # string replace retains format New.ReferenceOTUxxx
            score = v[3]
            otu_class = v[2]
            out_data[otu] = [score, otu_class]
        indcl = cls.index(v[2])
        lab = str(v[2]) if str(v[2]) not in added else ""
        added.append(str(v[2])) 
        col = colors[indcl%len(colors)] 
        if len(cls2) > 0: 
            col = colors[cls2.index(v[2])%len(colors)]
        vv = fabs(float(v[3])) * (m*(indcl*2-1)) if bcl else fabs(float(v[3]))
        ax.barh(pos[i],vv, align='center', color=col, label=lab, height=0.8, edgecolor=params['fore_color'])
    mv = max([abs(float(v[3])) for v in data['rows']])  
    if report_features:
        print('OTU\tLDA_score\tCLass')
        for i in out_data:
            print('%s\t%s\t%s' %(i, out_data[i][0], out_data[i][1]))
    for i,r in enumerate(data['rows']):
        indcl = cls.index(data['rows'][i][2])
        if params['n_scl'] < 0: rr = r[0]
        else: rr = ".".join(r[0].split(".")[-params['n_scl']:])
        if len(rr) > params['max_feature_len']: rr = rr[:params['max_feature_len']/2-2]+" [..]"+rr[-params['max_feature_len']/2+2:]
        if m*(indcl*2-1) < 0 and bcl: ax.text(mv/40.0,float(i)-0.3,rr, l_align, size=params['feature_font_size'],color=params['fore_color'])
        else: ax.text(-mv/40.0,float(i)-0.3,rr, r_align, size=params['feature_font_size'],color=params['fore_color'])
    ax.set_title(params['title'],size=params['title_font_size'],y=1.0+head*(1.0-ints/(ints+ht))*0.8,color=params['fore_color'])

    ax.set_yticks([])
    ax.set_xlabel("LDA SCORE (log 10)")
    ax.xaxis.grid(True)
    xlim = ax.get_xlim()
    if params['autoscale']: 
        ran = arange(0.0001,round(round((abs(xlim[0])+abs(xlim[1]))/10,4)*100,0)/100)
        if len(ran) > 1 and len(ran) < 100:
            ax.set_xticks(arange(xlim[0],xlim[1]+0.0001,min(xlim[1]+0.0001,round(round((abs(xlim[0])+abs(xlim[1]))/10,4)*100,0)/100)))
    ax.set_ylim((pos[0]-1,pos[-1]+1))
    leg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, borderaxespad=0., frameon=False,prop={'size':params['class_legend_font_size']})

    def get_col_attr(x):
                return hasattr(x, 'set_color') and not hasattr(x, 'set_facecolor')
    for o in leg.findobj(get_col_attr):
                o.set_color(params['fore_color'])
    for o in ax.findobj(get_col_attr):
                o.set_color(params['fore_color'])

    
    plt.savefig(path,format=params['format'],facecolor=params['back_color'],edgecolor=params['fore_color'],dpi=params['dpi'])
    plt.close()

def plot_histo_ver(path,params,data,report_features):
    cls = data['cls']
    mmax = max([fabs(float(a)) for a in zip(*data['rows'])[1]])
    data['rows'].sort(key=lambda ab: fabs(float(ab[3]))/mmax+(cls.index(ab[2])+1))
    pos = arange(len(data['rows'])) 
    if params['n_scl'] < 0: nam = [d[0] for d in data['rows']]
    else: nam = [d[0].split(".")[-min(d[0].count("."),params['n_scl'])] for d in data['rows']]
    fig = plt.figure(edgecolor=params['back_color'],facecolor=params['back_color'],figsize=(params['width'], params['height'])) 
    ax = fig.add_subplot(111,facecolor=params['back_color'])
    plt.subplots_adjust(top=0.9, left=params['ls'], right=params['rs'], bottom=0.3) 
    fig.canvas.manager.set_window_title('LDA results')   
    l_align = {'horizontalalignment':'left', 'verticalalignment':'baseline'}
    r_align = {'horizontalalignment':'right', 'verticalalignment':'baseline'} 
    added = []
    out_data = defaultdict(list) # keep track of which OTUs result in the plot
    for i,v in enumerate(data['rows']):
        if report_features:
            otu = v[0].split('.')[7].replace('_','.') # string replace retains format New.ReferenceOTUxxx
            score = v[3]
            otu_class = v[2]
            out_data[otu] = [score, otu_class]
        indcl = data['cls'].index(v[2])
        lab = str(v[2]) if str(v[2]) not in added else ""
        added.append(str(v[2])) 
        col = colors[indcl%len(colors)]
        vv = fabs(float(v[3])) 
        ax.bar(pos[i],vv, align='center', color=col, label=lab)
    if report_features:
        print('OTU\tLDA_score\tCLass')
        for i in out_data:
            print('%s\t%s\t%s' %(i, out_data[i][0], out_data[i][1]))
    xticks(pos,nam,rotation=-20, ha = 'left',size=params['feature_font_size'])  
    ax.set_title(params['title'],size=params['title_font_size'])
    ax.set_ylabel("LDA SCORE (log 10)")
    ax.yaxis.grid(True) 
    a,b = ax.get_xlim()
    dx = float(len(pos))/float((b-a))
    ax.set_xlim((0-dx,max(pos)+dx)) 
    plt.savefig(path,format=params['format'],facecolor=params['back_color'],edgecolor=params['fore_color'],dpi=params['dpi'])
    plt.close() 

def plot_res():
    params = read_params(sys.argv)
    params['fore_color'] = 'w' if params['back_color'] == 'k' else 'k'
    data = read_data(params['input_file'],params['output_file'],params['otu_only'])
    if params['orientation'] == 'v': plot_histo_ver(params['output_file'],params,data,params['report_features'])
    else: plot_histo_hor(params['output_file'],params,data,len(data['cls']) == 2,params['report_features'])


if __name__ == '__main__':
    plot_res()