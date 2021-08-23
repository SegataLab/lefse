#!/usr/bin/env python3

import os,sys,matplotlib,zipfile,argparse,string
matplotlib.use('Agg')
from pylab import *
from lefse.lefse import *
import random as rand

colors = ['r','g','b','m','c']

def read_params(args):
	parser = argparse.ArgumentParser(description='Cladoplot')
	parser.add_argument('input_file_1', metavar='INPUT_FILE', type=str, help="dataset files")
	parser.add_argument('input_file_2', metavar='INPUT_FILE', type=str, help="LEfSe output file")
	parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str, help="the file for the output (the zip file if an archive is required, the output directory otherwise)")
	parser.add_argument('--width',dest="width", type=float, default=10.0 )
	parser.add_argument('--height',dest="height", type=float, default=4.0) 
	parser.add_argument('--top',dest="top", type=float, default=-1.0, help="set maximum y limit (-1.0 means automatic limit)") 
	parser.add_argument('--bot',dest="bot", type=float, default=0.0, help="set minimum y limit (default 0.0, -1.0 means automatic limit)") 
	parser.add_argument('--title_font_size',dest="title_font_size", type=str, default="14")
	parser.add_argument('--class_font_size',dest="class_font_size", type=str, default="14")
	parser.add_argument('--class_label_pos',dest="class_label_pos", type=str, choices=["up","down"], default="up")
	parser.add_argument('--subcl_mean',dest="subcl_mean", type=str, choices=["y","n"], default="y")
	parser.add_argument('--subcl_median',dest="subcl_median", type=str, choices=["y","n"], default="y")
	parser.add_argument('--font_size',dest="font_size", type=str, default="10")
	parser.add_argument('-n',dest="unused", metavar="flt", type=float, default=-1.0,help="unused")
	parser.add_argument('--format', dest="format", default="png", choices=["png","pdf","svg"], type=str, help="the format for the output file")
	parser.add_argument('-f', dest="f", default="diff", choices=["all","diff","one"], type=str, help="wheter to plot all features (all), only those differentially abundant according to LEfSe or only one (the one given with --feature_name) ")
	parser.add_argument('--feature_name', dest="feature_name", default="", type=str, help="The name of the feature to plot (levels separated by .) ")
	parser.add_argument('--feature_num', dest="feature_num", default="-1", type=int, help="The number of the feature to plot ")
	parser.add_argument('--archive', dest="archive", default="none", choices=["zip","none"], type=str, help="")
	parser.add_argument('--background_color',dest="back_color", type=str, choices=["k","w"], default="w", help="set the color of the background")
	parser.add_argument('--dpi',dest="dpi", type=int, default=72)
		
	args = parser.parse_args()

	return vars(args)
	
def read_data(file_data,file_feats,params):
	with open(file_feats, 'r') as features:
		feats_to_plot = [(f.split()[:-1],len(f.split()) == 5) for f in features.readlines()]
	if not feats_to_plot:
		print("No features to plot\n")
		sys.exit(0)
	feats,cls,class_sl,subclass_sl,class_hierarchy,params['norm_v'] = load_data(file_data, True)	 	
	if params['feature_num'] > 0: 
		params['feature_name'] = [line.split()[0] for line in open(params['input_file_2'])][params['feature_num']-1]
	features = {}
	for f in feats_to_plot:
		if params['f'] == "diff" and not f[1]: continue
		if params['f'] == "one" and f[0][0] != params['feature_name']: continue
		features[f[0][0]] = {'dim':float(f[0][1]), 'abundances':feats[f[0][0]], 'sig':f[1], 'cls':cls, 'class_sl':class_sl, 'subclass_sl':subclass_sl, 'class_hierarchy':class_hierarchy} 
	if not features:
                print("No features to plot\n")
                sys.exit(0)
	return features

def plot(name,k_n,feat,params):
	fig = plt.figure(figsize=(params['width'], params['height']),edgecolor=params['fore_color'],facecolor=params['back_color'])
	ax = fig.add_subplot(111,facecolor=params['back_color']) 
	subplots_adjust(bottom=0.15)

	max_m = 0.0
	norm = 1.0 if float(params['norm_v']) < 0.0 else float(params['norm_v'])
	for v in feat['subclass_sl'].values():
		fr,to  = v[0], v[1]
		median = numpy.mean(feat['abundances'][fr:to])
		if median > max_m: max_m = median
	max_m /= norm
	max_v = max_m*3 if max_m*3 < max(feat['abundances'])*1.1/norm else max(feat['abundances'])/norm
	min_v = max(0.0,min(feat['abundances'])*0.9/norm)

	if params['top'] > 0.0: max_v = params['top'] 
	if params['bot'] >= 0.0: min_v = params['bot']
	
	if max_v == 0.0: max_v = 0.0001
	if max_v == min_v: max_v = min_v*1.1 

	cl_sep = max(int(sum([vv[1]/norm - vv[0]/norm for vv in feat['class_sl'].values()])/150.0),1)
	seps = []
	xtics = []
	x2tics = []
	last_fr = 0.0
	for i,cl in enumerate(sorted(feat['class_hierarchy'].keys())):
		for j,subcl in enumerate(feat['class_hierarchy'][cl]):
			fr = feat['subclass_sl'][subcl][0]
			to = feat['subclass_sl'][subcl][1]
			val = feat['abundances'][fr:to]
			fr += cl_sep*i
			to += cl_sep*i
			pos = arange(fr,to)
			max_x = to
			col = colors[j%len(colors)]
			vv = [v1/norm for v1 in val]
			median = numpy.median(vv)
			mean = numpy.mean(vv)
			valv = [max(min(v/norm,max_v),min_v) for v in val]
			ax.bar(pos,valv, align='center', color=col, edgecolor=col, linewidth=0.1 )
			if params['subcl_median'] == 'y': ax.plot([fr,to-1],[median,median],"k--",linewidth=1,color=params['fore_color'])
			if params['subcl_mean'] == 'y': ax.plot([fr,to-1],[mean,mean],"-",linewidth=1,color=params['fore_color'])
			nna = subcl if subcl.count("_") == 0 or not subcl.startswith(cl) else "_".join(subcl.split(cl)[1:])
			if nna == "subcl" or nna == "_subcl": nna = " "
			xtics.append(((fr+(to-fr)/2),nna))
		seps.append(float(to))
		x2tics.append(((last_fr+(to-last_fr)/2),cl))
		last_fr = to + float(cl_sep)
	for s in seps[:-1]:
		ax.plot([s,s],[min_v,max_v],"-",linewidth=5,color=params['fore_color'])	
	ax.set_title(k_n, size=params['title_font_size'])
	xticks([x[0] for x in xtics],[x[1] for x in xtics],rotation=-30, ha = 'left', fontsize=params['font_size'], color=params['fore_color'])
	yticks(fontsize=params['font_size'])

	ylabel('Relative abundance')
	ax.set_ylim((min_v,max_v))
	a,b = ax.get_xlim()
	ax.set_xlim((0-float(last_fr)/float(b-a),max_x))		
	ax.yaxis.grid(True)
	
	def get_col_attr(x):
                return hasattr(x, 'set_color') and not hasattr(x, 'set_facecolor')
	def get_edgecol_attr(x):
                return hasattr(x, 'set_edgecolor') 


	for o in fig.findobj(get_col_attr):
                o.set_color(params['fore_color'])	
	for o in fig.findobj(get_edgecol_attr):
		if o.get_edgecolor() == params['back_color']:
                	o.set_edgecolor(params['fore_color'])

	for t in x2tics:
		m = ax.get_ylim()[1]*0.97 if params['class_label_pos']=='up' else 0.07
		plt.text(t[0],m, "class: "+t[1], ha ="center", size=params['class_font_size'], va="top", bbox = dict(boxstyle="round", ec='k', fc='y'))
	

	plt.savefig(name,format=params['format'],facecolor=params['back_color'],edgecolor=params['fore_color'],dpi=params['dpi'])
	plt.close()
	return name 

def plot_features():
	params = read_params(sys.argv)
	params['fore_color'] = 'w' if params['back_color'] == 'k' else 'k'
	features = read_data(params['input_file_1'],params['input_file_2'],params)
	if params['archive'] == "zip": file = zipfile.ZipFile(params['output_file'], "w")
	for k,f in features.items():
		print("Exporting ", k)
		if params['archive'] == "zip":
			of = plot("/tmp/"+str(int(f['sig']))+"_"+"-".join(k.split("."))+"."+params['format'],k,f,params)
			file.write(of, os.path.basename(of), zipfile.ZIP_DEFLATED)
		else:
			if params['f'] == 'one': plot(params['output_file'],k,f,params)
			else: plot(params['output_file']+str(int(f['sig']))+"_"+"-".join(k.split("."))+"."+params['format'],k,f,params) 
	if params['archive'] == "zip": file.close()


if __name__ == '__main__':
	plot_features()