import subprocess
import os, sys, getopt, time
import argparse
import pandas as pd
from Bio.Seq import Seq

def runMenthu(args):
    #subprocess.run(["Rscript", "ExampleRScript.R"])
    subprocess.run(["Rscript", "Menthu-cmd.R", args.Outfile, args.Crispr, args.Pam, args.Dist2dsb, args.Overhang,
                       args.TalenOption, args.TalenScheme, args.GenInputType, args.GenInput, args.ScoreThreshold,
                        args.T7opt, args.Verbose, args.Validate])
    time_to_wait = 20
    time_counter = 0
    while not os.path.exists(args.Outfile):
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            break
    runLindel(args.Outfile)


def runLindel(outfile):
    df = pd.read_csv(outfile)
    i = 0
    for j in df['Tool_Type']:
        #print(j)
        if j == 'NGG' and df['Strand'][i] == 'forward':    #Change tool type to NGG!
            dna_string = df['Context'][i][10::]
            dna_string = dna_string[:-10]
            #print(dna_string)
        elif j == 'NGG' and df['Strand'][i] == 'complement':   #Change tool type to NGG
            dna_string = Seq(df['Context'][i]).reverse_complement()[10::]
            dna_string = dna_string[:-10]
            #print(dna_string)
        subprocess.run(["python", "Lindel_prediction.py", dna_string, "test_out"])
        i += 1

    #something here

# Rscript menthu.R [outFile] [CRISPR Option] [PAM Sequence] [Distance to DSB] [Overhang] [TALEN Option] [TALEN
# scheme] [Gen Input Type] [Gen Input] [Score Threshold] [T7 opt] [verbose] [validate]

def main():
    parser = argparse.ArgumentParser(description='Inputs arguments for input to Menthu')
    parser.add_argument('-o', '--Outfile', type=str, default='menthu_outfile_2.csv',
                        help='character string file name The name of the file to output your results to. If using a '
                             'fasta file with multiple sequences, multiple files will be created, using this as a '
                             'prefix')
    parser.add_argument('-c', '--Crispr', type=str, default='T',
                        help='T or F Flags the system to use CRISPR nuclease processing. If this option is T (true), '
                             '"TALEN Option" must be F (false)')
    parser.add_argument('-p', '--Pam', type=str, default='NGG',
                        help='A PAM sequence The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are '
                             'allowed. Using N will scan every possible cut site in the target sequence. This '
                             'parameter must be present, but is not used, if "CRISPR Option" is false (i.e., '
                             'you can put a 0 or NA in this spot.)')
    parser.add_argument('-d', '--Dist2dsb', type=str, default='-3',
                        help='Integer The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, '
                             'use a negative value (e.g., -3 for SpCa9); for downstream, use a positive value (e.g., '
                             '18 for Cas12a.) This parameter must be present, but is not used, if "CRISPR Option" is '
                             'false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-oh', '--Overhang', type=str, default='0',
                        help='Integer >= 0 The length of 5\' overhang produced by the nuclease (e.g., 5 for Cas12a). '
                             'Use 0 for blunt-cutting nucleases. This parameter must be present, but is not used, '
                             'if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-to', '--TalenOption', type=str, default='F',
                        help='T or F Flags the system to use TALEN processing. If this option is T (true), "CRISPR '
                             'Option" must be F (false)')
    parser.add_argument('-ts', '--TalenScheme', type=str, default='0',
                        help='15-18/14 or 16/15-18 The left arm length, spacer length, and right arm length to use '
                             'when searching for TALEN locations. E.g., for a TALEN with arms 15 nt long, with spacer '
                             '14 nt, use 15/14/15. TALEN arms can be 15-18 nt in length; the spacer should be 14 OR '
                             '16 nt in length (15 is not allowed for the spacer) This parameter must be present, '
                             'but is not used, if "TALEN Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-g', '--GenInputType', type=str, default='gb',
                        help='gb ens seq file Flags the system to get a GenBank/RefSeq ID (gb), Ensembl ID (ens), '
                             'DNA sequence (seq), or to expect a FASTA file (file)')
    parser.add_argument('-i', '--GenInput', type=str, default='AY214391.1',
                        help='See explanation Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for '
                             '"file" option, or DNA sequence for "seq". If the file name has spaces in it, '
                             'put this parameter in quotes.')
    parser.add_argument('-st', '--ScoreThreshold', type=str, default='1',
                        help='Positive number Only output results with MENTHU score above this threshold. Default is '
                             '1.0. We recommend to only use sites with score >= 1.5')
    parser.add_argument('-t7', '--T7opt', type=str, default='F',
                        help='T or F If T (true), only displays results where the gRNA is compatible with T7-cloning.')
    parser.add_argument('-v', '--Verbose', type=str, default='F',
                        help='T or F If T (true), outputs progress messages to the console.')
    parser.add_argument('-va', '--Validate', type=str, default='F',
                        help='T or F If T (true), checks the command line arguments to make sure they are all valid ('
                             'this may take some time); if F, skip validation checks')
    args = parser.parse_args()
    runMenthu(args)


if __name__ == '__main__':
    main()

#subprocess.run(["Rscript", "ExampleRScript.R"])
