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
   python MENdelScript.py [outFile] [CRISPR Option] [PAM Sequence] [Distance to DSB] [Overhang] [TALEN Option] [TALEN scheme] [Gen Input Type] [Gen Input] [Score Threshold] [T7 opt] [verbose] [validate]
   ```
   
"Rscript" tells the system to use ```R``` to execute ```menthu.R``` (or ```menthu-cmd.R``` for Windows CMD).
