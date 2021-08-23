#!/usr/bin/env python3

import os,sys,math,pickle
from lefse.lefse import *

def read_params(args):
    parser = argparse.ArgumentParser(description='LEfSe 1.1.01')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="the input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str,
                help="the output file containing the data for the visualization module")
    parser.add_argument('-o',dest="out_text_file", metavar='str', type=str, default="",
                help="set the file for exporting the result (only concise textual form)")
    parser.add_argument('-a',dest="anova_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Anova test (default 0.05)")
    parser.add_argument('-w',dest="wilcoxon_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Wilcoxon test (default 0.05)")
    parser.add_argument('-l',dest="lda_abs_th", metavar='float', type=float, default=2.0,
                help="set the threshold on the absolute value of the logarithmic LDA score (default 2.0)")
    parser.add_argument('--nlogs',dest="nlogs", metavar='int', type=int, default=3,
        help="max log ingluence of LDA coeff")
    parser.add_argument('--verbose',dest="verbose", metavar='int', choices=[0,1], type=int, default=0,
        help="verbose execution (default 0)")
    parser.add_argument('--wilc',dest="wilc", metavar='int', choices=[0,1], type=int, default=1,
        help="wheter to perform the Wicoxon step (default 1)")
    parser.add_argument('-r',dest="rank_tec", metavar='str', choices=['lda','svm'], type=str, default='lda',
        help="select LDA or SVM for effect size (default LDA)")
    parser.add_argument('--svm_norm',dest="svm_norm", metavar='int', choices=[0,1], type=int, default=1,
        help="whether to normalize the data in [0,1] for SVM feature waiting (default 1 strongly suggested)")
    parser.add_argument('-b',dest="n_boots", metavar='int', type=int, default=30,
                help="set the number of bootstrap iteration for LDA (default 30)")
    parser.add_argument('-e',dest="only_same_subcl", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test only among the subclasses with the same name (default 0)")
    parser.add_argument('-c',dest="curv", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test ing the Curtis's approach [BETA VERSION] (default 0)")
    parser.add_argument('-f',dest="f_boots", metavar='float', type=float, default=0.67,
                help="set the subsampling fraction value for each bootstrap iteration (default 0.66666)")
    parser.add_argument('-s',dest="strict", choices=[0,1,2], type=int, default=0,
                help="set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for dependent comparison")
#       parser.add_argument('-m',dest="m_boots", type=int, default=5,
#               help="minimum cardinality of classes in each bootstrapping iteration (default 5)")
    parser.add_argument('--min_c',dest="min_c", metavar='int', type=int, default=10,
                help="minimum number of samples per subclass for performing wilcoxon test (default 10)")
    parser.add_argument('-t',dest="title", metavar='str', type=str, default="",
                help="set the title of the analysis (default input file without extension)")
    parser.add_argument('-y',dest="multiclass_strat", choices=[0,1], type=int, default=0,
                help="(for multiclass tasks) set whether the test is performed in a one-against-one ( 1 - more strict!) or in a one-against-all setting ( 0 - less strict) (default 0)")
    args = parser.parse_args()

    params = vars(args)
    if params['title'] == "":
        params['title'] = params['input_file'].split("/")[-1].split('.')[0]

    return params


def lefse_run():
    init()
    params = read_params(sys.argv)
    feats,cls,class_sl,subclass_sl,class_hierarchy = load_data(params['input_file'])
    kord,cls_means = get_class_means(class_sl,feats)
    wilcoxon_res = {}
    kw_n_ok = 0
    nf = 0
    for feat_name,feat_values in list(feats.items()):
        if params['verbose']:
            print("Testing feature",str(nf),": ",feat_name)
            nf += 1
        kw_ok,pv = test_kw_r(cls,feat_values,params['anova_alpha'],sorted(cls.keys()))
        if not kw_ok:
            if params['verbose']: print("\tkw ko")
            del feats[feat_name]
            wilcoxon_res[feat_name] = "-"
            continue
        if params['verbose']: print("\tkw ok\t")

        if not params['wilc']: continue
        kw_n_ok += 1
        res_wilcoxon_rep = test_rep_wilcoxon_r(subclass_sl,class_hierarchy,feat_values,params['wilcoxon_alpha'],params['multiclass_strat'],params['strict'],feat_name,params['min_c'],params['only_same_subcl'],params['curv'])
        wilcoxon_res[feat_name] = str(pv) if res_wilcoxon_rep else "-"
        if not res_wilcoxon_rep:
            if params['verbose']: print("wilc ko")
            del feats[feat_name]
        elif params['verbose']: print("wilc ok\t")

    if len(feats) > 0:
        print("Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon")
        if params['lda_abs_th'] < 0.0:
            lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
        else:
            if params['rank_tec'] == 'lda': lda_res,lda_res_th = test_lda_r(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0000000001,params['nlogs'])
            elif params['rank_tec'] == 'svm': lda_res,lda_res_th = test_svm(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0,params['svm_norm'])
            else: lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
    else:
        print("Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon")
        print("No features with significant differences between the two classes")
        lda_res,lda_res_th = {},{}
    outres = {}
    outres['lda_res_th'] = lda_res_th
    outres['lda_res'] = lda_res
    outres['cls_means'] = cls_means
    outres['cls_means_kord'] = kord
    outres['wilcox_res'] = wilcoxon_res
    print("Number of discriminative features with abs LDA score >",params['lda_abs_th'],":",len(lda_res_th))
    save_res(outres,params["output_file"])


if __name__ == '__main__':
    lefse_run()