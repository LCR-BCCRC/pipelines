import argparse
import os.path
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as mg

def main():
    """ Generates matrix from MAF file """
    args = parse_arguments()
    outdir = os.path.dirname(args.maf) # output made in MAF directory

    # Generates matrix from MAF found in outdir
    matrices = mg.SigProfilerMatrixGeneratorFunc(args.project, args.ref, \
                                                 outdir, tsb_stat = True)
    return()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('project', help = 'Project name. Will be file prefixes.')
    parser.add_argument('ref', help = 'Genome reference of input file')
    parser.add_argument('maf', help = 'Input MAF')
    return(parser.parse_args())

if __name__ == '__main__':
    main()
