#! /usr/bin/env python
import re, os
import getopt, sys, subprocess
import shutil


# defaults
settings=dict(nproc=2, memory='500MB')

try:
    opts, args = getopt.getopt(sys.argv[1:], "hvn:st:", ["help", "nproc=", "saveonly", "template="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

    # OK, now we get to the part where I know what's going on.
    # for each file supplied at the command line, do the following

searchExpression1=re.compile(" ! ")
searchExpression2=re.compile(" !")
for filename in args:

    (fileroot,filextension) = os.path.splitext(filename)
    
    print fileroot
            
    fin=open(filename,'r')
    com=fin.read()
    fin.close

    replaceme = ('DUPLICATE') #RAS07 ING/BOZ03

    
    com = com.replace(replaceme, '  //  DUPLICATE ')
    
    name = filename
    fout=open(name,'w')
    fout.write(com)
    fout.close
