import argparse
from SigProfilerExtractor import sigpro as sig
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1" # ulimit fix

def main():
    """ Performs NMF decomposition with provided matrix """
    args = parse_arguments()
    sig.sigProfilerExtractor(input_type = 'matrix', \
                             output = args.output, \
                             input_data = args.input_data, \
                             reference_genome = args.ref, \
                             opportunity_genome = args.ref, \
                             context_type = args.context_type, \
                             exome = args.exome, \
                             minimum_signatures = args.minimum_signatures, \
                             maximum_signatures = args.maximum_signatures, \
                             nmf_replicates = args.nmf_replicates, \
                             matrix_normalization = args.matrix_normalization, \
                             nmf_init = args.nmf_init, \
                             cpu = args.cpu, \
                             make_decomposition_plots = False, \
                             resample = args.resample, \
                             precision = args.precision, \
                             min_nmf_iterations = args.min_nmf_iterations, \
                             max_nmf_iterations = args.max_nmf_iterations, \
                             nmf_test_conv = args.nmf_test_conv, \
                             nmf_tolerance = args.nmf_tolerance, \
                             stability = args.stability, \
                             min_stability = args.min_stability, \
                             combined_stability = args.combined_stability, \
                             cosmic_version = args.cosmic_version, \
                             de_novo_fit_penalty = args.de_novo_fit_penalty, \
                             nnls_add_penalty = args.nnls_add_penalty, \
                             nnls_remove_penalty = args.nnls_remove_penalty, \
                             initial_remove_penalty = args.initial_remove_penalty, \
                             refit_denovo_signatures = args.refit_denovo_signatures
                             )
    return()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_data', help = 'Path to matrix that will undergo NMF.')
    parser.add_argument('output', help = 'Path to output directory.')
    parser.add_argument('ref', help = 'Reference build of matrix.')
    parser.add_argument('context_type', default = '96,DINUC,ID', help = 'A string of mutaion context name/names separated by comma (","). The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUC,ID", where "96" is the SBS96 context, "DINUC" is the DINUCLEOTIDE context and ID is INDEL context.')
    parser.add_argument('exome', default = False, type = bool, help = 'Defines if the exomes will be extracted. The default value is "%(default)s".')
    parser.add_argument('minimum_signatures', default = 1, type = int, help = "The minimum number of signatures to be extracted. The default value is %(default)s.")
    parser.add_argument('maximum_signatures', default = 25, type = int, help = "The maximum number of signatures to be extracted. The default value is %(default)s.")
    parser.add_argument('cpu', default = -1, type = int, help = 'Number of threads to use. Default use all threads.')
    parser.add_argument('--nmf_replicates', default = 100 , type = int, help = "The number of iteration to be performed to extract each number signature. The default value is %(default)s.")
    parser.add_argument('--matrix_normalization', default = 'gmm', choices = ['gmm','log2','custom','none'], help = 'Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "%(default)s". Other options are, "log2", "custom" or "none".')
    parser.add_argument('--nmf_init', default = 'nndsvd_min', choices = ['nndsvd_min','random','nndsvd','nndsvda','nndsvdar'], help = "The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'. Default is '%(default)s'.")
    parser.add_argument('--resample', default = True, type = bool, help = 'Default is %(default)s. If True, add poisson noise to samples by resampling.')
    parser.add_argument('--precision', default = 'single', choices = ['single','double'], help = 'Values should be single or double. Default is %(default)s.')
    parser.add_argument('--min_nmf_iterations', default = 10000, type = int, help = 'Value defines the minimum number of iterations to be completed before NMF converges. Default is %(default)s.')
    parser.add_argument('--max_nmf_iterations', default = 1000000, type = int, help = 'Value defines the maximum number of iterations to be completed before NMF converges. Default is %(default)s.')
    parser.add_argument('--nmf_test_conv', default = 10000, type = int, help = 'Value defines the number of iterations to done between checking next convergence. Default is %(default)s.')
    parser.add_argument('--nmf_tolerance', default = 1e-15, type = float, help = 'Value defines the tolerance to achieve to converge. Default is %(default)s.')
    parser.add_argument('--stability', default = 0.8, type = float, help = 'Default is %(default)s. The cutoff thresh-hold of the average stability. Solutions with average stabilities below this thresh-hold will not be considered.')
    parser.add_argument('--min_stability', default = 0.2, type = float, help = 'Default is %(default)s. The cutoff thresh-hold of the minimum stability. Solutions with minimum stabilities below this thresh-hold will not be considered.')
    parser.add_argument('--combined_stability', default = 1, type = float, help = 'Default is %(default)s. The cutoff thresh-hold of the combined stability (sum of average and minimum stability). Solutions with combined stabilities below this thresh-hold will not be considered.')
    parser.add_argument('--cosmic_version', default = 3.1, type = float, choices = [1,2,3,3.1,3.2], help = 'Takes a positive float among 1, 2, 3, 3.1, 3.2. Default is %(default)s. Defines the version of COSMIC reference signatures.')
    parser.add_argument('--de_novo_fit_penalty', default = 0.02, type = float, help = 'Takes any positive float. Default is %(default)s. Defines the weak (remove) thresh-hold cutoff to assign denovo signatures to a sample.')
    parser.add_argument('--nnls_add_penalty', default = 0.05, type = float, help = 'Takes any positive float. Default is %(default)s. Defines the strong (add) thresh-hold cutoff to assign COSMIC signatures to a sample.')
    parser.add_argument('--nnls_remove_penalty', default = 0.01, type = float, help = 'Takes any positive float. Default is %(default)s. Defines the weak (remove) thresh-hold cutoff to assign COSMIC signatures to a sample.')
    parser.add_argument('--initial_remove_penalty', default = 0.05, type = float, help = 'Takes any positive float. Default is %(default)s. Defines the initial weak (remove) thresh-hold cutoff to COSMIC assign signatures to a sample.')
    parser.add_argument('--refit_denovo_signatures', default = True, type = bool, help = 'Default is %(default)s. If True, then refit the denovo signatures with nnls.')
    
    return(parser.parse_args())

if __name__ == '__main__':
    main()
