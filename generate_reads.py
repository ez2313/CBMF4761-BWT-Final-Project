'''Script to randomly generate a reference of specified length "reference.txt",
draw specified number of reads and length output to "reads.txt".'''
import numpy as np
import pandas as pd
import random as rand
import sys


def simulate_read(read_length,ref_length,number_of_reads):
  ref=''
  for i in range(ref_length):
    toss=rand.random()
    if toss<.25:
      ref=ref+"A"
    if .25<=toss and toss<.5:
      ref=ref+"C"
    if .5<=toss and toss<.75:
      ref=ref+"G"
    else:
      ref=ref+"T"
  reads=[]
  for i in range(number_of_reads):
      guess=rand.randint(0, ref_length-read_length)
      read=ref[guess:guess+read_length]
      reads.append(read)
  return reads, ref

def main():
    '''randomly generate reference, and reads of specified length and number '''
    #read in input information
    number_of_reads=int(sys.argv[1])
    read_length=int(sys.argv[2])
    reference_length=int(sys.argv[3])
    reads,ref=simulate_read(read_length,reference_length,number_of_reads)


    #output results to csv, text
    reads_df=pd.DataFrame(reads)
    reads_df.to_csv('reads.txt',header=None,index=None)
    reference_output=open('reference.txt','w')
    reference_output.write(ref)
    reference_output.close()



if __name__ == "__main__":
    main()
