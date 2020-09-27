import subprocess
import os, sys, getopt, time, shutil, glob
import argparse
import pandas as pd
from Bio.Seq import Seq

def runMenthu(args):
    #subprocess.run(["Rscript", "ExampleRScript.R"])
    menthudir = ''
    lindeldir = ''
    for i in os.listdir():
        if "MENTHU" in i:
            menthudir = i
        elif "Lindel" in i:
            lindeldir = i
    os.chdir(menthudir)
    subprocess.run(["Rscript", "menthu.R", args.Outfile, args.crispr, args.pam, args.dist2dsb, args.overhang,
                       args.talenOption, args.talenScheme, args.genInputType, args.genInput, args.scoreThreshold,
                        args.T7opt, args.verbose, args.validate])
    shutil.move(args.Outfile, "..\\"+args.Outfile)
    os.chdir("..")
    time_to_wait = 20
    time_counter = 0
    while not os.path.exists(args.Outfile):
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            break
    runLindel(args.Outfile, lindeldir)


def runLindel(outfile, lindeldir):
    shutil.move(outfile, lindeldir+"\\"+outfile)
    os.chdir(lindeldir)
    df = pd.read_csv(outfile)
    i = 0
    for j in df['Tool_Type']:
        #print(df.iloc[i])
        #print(j)
        if j == 'NGG' and df['Strand'][i] == 'forward':
            dna_string = df['Context'][i][10::]
            dna_string = dna_string[:-10]
            #print(dna_string)
            #print(df['MENTHU_Score'][i])
        elif j == 'NGG' and df['Strand'][i] == 'complement':
            dna_string = Seq(df['Context'][i]).reverse_complement()[10::]
            dna_string = dna_string[:-10]
            #print(dna_string)
            #print(df['MENTHU_Score'][i])
        subprocess.run(["python", "Lindel_prediction.py", dna_string, "test_out"])
        txt = glob.glob(os.getcwd()+"\\test_out*.txt")
        for txtfile in txt:
            shutil.move(txtfile, str(outfile[:-4])+"_"+str(i+1)+".tsv")
            fh = open(str(outfile[:-4])+"_"+str(i+1)+".tsv", "a")
            fh.seek(0)
            fh.write(str(df.iloc[i]))
            fh.close()
        shutil.move(str(outfile[:-4])+"_"+str(i+1)+".tsv", "..\\"+str(outfile[:-4])+"_"+str(i+1)+".tsv")
        i += 1
    shutil.move(outfile, "..\\"+outfile)


def main():
    parser = argparse.ArgumentParser(description='Inputs arguments for input to Menthu')
    parser.add_argument('-o', '--Outfile', type=str, default='MENdel_outfile.csv',
                        help='character string file name The name of the file to output your results to. If using a '
                             'fasta file with multiple sequences, multiple files will be created, using this as a '
                             'prefix')
    parser.add_argument('-c', '--crispr', type=str, default='T',
                        help='T or F Flags the system to use CRISPR nuclease processing. If this option is T (true), '
                             '"TALEN Option" must be F (false)')
    parser.add_argument('-p', '--pam', type=str, default='NGG',
                        help='A PAM sequence The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are '
                             'allowed. Using N will scan every possible cut site in the target sequence. This '
                             'parameter must be present, but is not used, if "CRISPR Option" is false (i.e., '
                             'you can put a 0 or NA in this spot.)')
    parser.add_argument('-d', '--dist2dsb', type=str, default='-3',
                        help='Integer The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, '
                             'use a negative value (e.g., -3 for SpCa9); for downstream, use a positive value (e.g., '
                             '18 for Cas12a.) This parameter must be present, but is not used, if "CRISPR Option" is '
                             'false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-oh', '--overhang', type=str, default='0',
                        help='Integer >= 0 The length of 5\' overhang produced by the nuclease (e.g., 5 for Cas12a). '
                             'Use 0 for blunt-cutting nucleases. This parameter must be present, but is not used, '
                             'if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-to', '--talenOption', type=str, default='F',
                        help='T or F Flags the system to use TALEN processing. If this option is T (true), "CRISPR '
                             'Option" must be F (false)')
    parser.add_argument('-ts', '--talenScheme', type=str, default='0',
                        help='15-18/14 or 16/15-18 The left arm length, spacer length, and right arm length to use '
                             'when searching for TALEN locations. E.g., for a TALEN with arms 15 nt long, with spacer '
                             '14 nt, use 15/14/15. TALEN arms can be 15-18 nt in length; the spacer should be 14 OR '
                             '16 nt in length (15 is not allowed for the spacer) This parameter must be present, '
                             'but is not used, if "TALEN Option" is false (i.e., you can put a 0 or NA in this spot.)')
    parser.add_argument('-g', '--genInputType', type=str, default='gb',
                        help='gb ens seq file Flags the system to get a GenBank/RefSeq ID (gb), Ensembl ID (ens), '
                             'DNA sequence (seq), or to expect a FASTA file (file)')
    parser.add_argument('-i', '--genInput', type=str, default='AY214391.1',
                        help='See explanation Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for '
                             '"file" option, or DNA sequence for "seq". If the file name has spaces in it, '
                             'put this parameter in quotes.')
    parser.add_argument('-st', '--scoreThreshold', type=str, default='1',
                        help='Positive number Only output results with MENTHU score above this threshold. Default is '
                             '1.0. We recommend to only use sites with score >= 1.5')
    parser.add_argument('-t7', '--T7opt', type=str, default='F',
                        help='T or F If T (true), only displays results where the gRNA is compatible with T7-cloning.')
    parser.add_argument('-v', '--verbose', type=str, default='F',
                        help='T or F If T (true), outputs progress messages to the console.')
    parser.add_argument('-va', '--validate', type=str, default='F',
                        help='T or F If T (true), checks the command line arguments to make sure they are all valid ('
                             'this may take some time); if F, skip validation checks')
    args = parser.parse_args()
    runMenthu(args)


if __name__ == '__main__':
    main()

