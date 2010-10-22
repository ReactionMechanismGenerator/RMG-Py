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
    
    com = com.replace('=', ' = ')
    com = com.replace('+', ' + ')
    com = com.replace('( + ', ' (+')
    name = 'new_file.txt'
    fout=open(name,'w')
    fout.write(com)
    fout.close


    fin=open(name,'r+')
    f = open(name+'.new', 'w')
    for line in fin:
        match1 = searchExpression1.search(line)
        match2 = searchExpression2.search(line)
        if match1:
            (data, comment) = line.split(' ! ')
            f.write('// '+comment+data+'  0.0 0.0 0.0  '+'\n')
        elif match2:
            (data, comment) = line.split(' !')
            f.write(data+'\n')
        else:    
            f.write(line)    
    f.close        
            
 #   fin=open(name+'.new','r')
 #   com=fin.read()
 #   fin.close
    
 #   com = com.replace('!', '//')
#    com = com.replace('DUPLICATE', '   //     DUPLICATE')
 #   fout=open(name,'w')
 #   fout.write(com)
 #   fout.close
 #   os.remove(name+'.new')
