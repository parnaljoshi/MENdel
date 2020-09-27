# MENdel

This repository contains python scripts for running MENTHU and Lindel tools via commandline interface, to generate the MENdel output. MENdel is compatible with Linux, Mac, Windows (using Cygwin), and Windows CMD. In order to run MENdel, it is required that the user has MENTHU-command-line and Lindel cloned locally. Follow the instructions given below to install each of these tools.

```
mkdir MENdel_root
cd MENdel_root
```

1. **MENTHU-command-line installation**
   
   Follow instructions at https://github.com/parnaljoshi/MENTHU-command-line to install the commandline version of MENTHU within the previously created directory MENdel_root. The README file contains instructions to download and install R and necessary R packages to run MENTHU, if not already present.

2. **Lindel installation**

   Lindel requires Python installation.
   
   (Include instructions to install Python)
   
   Follow instructions at https://github.com/shendurelab/Lindel to install Lindel within the directory MENdel_root

3. **MENdel installation**

   Navigate to the previously created directory MENdel_root and follow these instructions to install scripts to run MENdel.
   
   ```
   git clone https://github.com/parnaljoshi/MENdel/
   cd MENdel
   ```
   
   Alternatively, use the green "Code" button in the upper right corner and choose "Download ZIP"

4. **Running MENdel**

   MENdel-command-line can be run from Unix-like command lines (Linux, Mac, Cygwin on Windows) or Windows CMD. However, due to end-of-line conversion differences between Unix and DOS for MENTHU, **you must use MENdelScript.py on Linux, Mac, and Cygwin, and Windows_MENdelScript.py in Windows CMD**. 

   Check to make sure you're using the correct command for your system!

   MENTHU-command-line can be run using the following syntax:
   
   - For Windows users
   
   ```
   python Windows_MENdelScript.py [outFile] [CRISPR Option] [PAM Sequence] [Distance to DSB] [Overhang] [TALEN Option] [TALEN scheme] [Gen Input Type] [Gen Input] [Score Threshold] [T7 opt] [verbose] [validate]
   ```
   
   - For Linux/Mac/Windows using WSL or Cygwin users
   
   ```
   python MENdelScript.py [-o outFile] [-c CRISPR Option] [-p PAM Sequence] [-d Distance to DSB] [-oh Overhang] [TALEN Option] [TALEN scheme] [Gen Input Type] [Gen Input] [Score Threshold] [T7 opt] [verbose] [validate]
   ```
   
   "python" tells the system to use ```python``` to execute ```MENdelScript.py``` (or ```Windows_MENdelScript.py``` for Windows CMD).
   
5. **Parameter explanation**

   The parameters are explained below. Each parameter is delimited by a space. Parameter values should not have spaces; if you want to put spaces in the output file name, the name should be in quotes, e.g. "output File.csv". Parameter values (including strings) do not have to be in quotes, except for the output file name exception. All parameters values can be kept to default with the exception of output file name, gen input type and gen input
   
   |**Parameter        |Accepted Values|Default|Explanation** |
   |-----------------  |---------------|-------|--------------|
   |**outFile**        |Character string followed by extension|MENdel_outfile.csv|The name of the file to output your results to. If using a fasta file with multiple sequences, multiple files will be created, using this as a prefix|
   |**CRSIPR Option**  |T or F         |T      |T or F Flags the system to use CRISPR nuclease processing. If this option is T (true) "TALEN Option" must be F (false)|
   |**PAM Sequence**   |A PAM sequence |NGG    |The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are allowed. Using N will scan every possible cut site in the target sequence. This parameter must be present, but is not used, if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)|
   |**Distance to DSB**|Integer      |-3     |The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, use a negative value (e.g., -3 for SpCa9); for downstream, use a positive value (e.g., 18 for Cas12a.) This parameter must be present, but is not used, if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)|

