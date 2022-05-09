# CBMF4761-BWT-Final-Project
CBMF4761 Final project

Welcome to the repository for my final project for CBMF4761!
The project is Burrows Wheeler Transform for Protein and DNA alignment. See ahead for instructions on using the aligner tool.

Please find included the following files:
1. aligner.py
2. FinalReport.pdf
3. generate_reads.py
4. model.ipynb
5. output.txt
6. paladin_testing.ipynb
7. reads.txt
8. reference.txt
9. sample_reads.txt
10. sample_reference.txt

# aligner.py
The aligner tool runs from the command line in python and takes 2 arguments: the first argument is the read dataset and the second is the reference dataset (where the files are in the same directory as aligner.py). There is an optional 3rd argument to specify the output location (the default is output.txt). <br>

An example of usage is: $python aligner.txt reads.txt reference.txt

# final report.pdf
The final report is a paper in the style of an Oxford University press Bioinformatics applications note. The report provides further information on the project and testing of the aligner tool. 

# generate_reads.py
This is a script used to create randomly generated reads and reference (such as used in the testing). The script is run from the command line with 3 arguments: (1) the number of reads, (2) the read length, (3) the reference length. <br>

An example of usage is: $python generate_reads.py 100 20 1000<br>

which would generate 100 reads each of length 20 from a reference of length 1000. 
The output of the script for the reads genereated is to reads.txt and for the reference is reference.txt

# model.ipynb
The notebook is where most of my testing and adapting occurred on Google Colab, with explanatory comments. 

# output.txt
This is a sample output file for the aligner tool. Each row corresponds to a single read of the form and gives the alignment score followed by the location of the alignment in the reference. (See ahead for usage).

# paladin_testing.ipynb
While this did not come to fruition, I spent a significant time exploring the Paladin framework in attempts to merge the aligner into the platform. The notebook allows for installing and making Paladin in Google Colab (which was necessary for me as Paladin had errors with Apple's M1 chip). 

# reads.txt
This is the output file of the script generate_reads.py giving the generated reads, one read per line. 

# reference.txt
This is the output file of the script generate_reads.py giving the generated reference. 

# sample_reads.txt
This is a sample input file of randomly generated reads using the script generate_reads.py (see ahead for usage).

## sample input and output
The provided sample_reads.txt and sample_reference.txt are sample input files for use with aligner.py to create output.txt.<br>
In future work and practical applications, the reads and references used would be provided by the studies themselves or taken from public resources such as NCBI or ECBI.<br>

Usage for the sample input and output is from the command line:<br>

$python aligner.txt sample_reads.txt sample_reference.txt
