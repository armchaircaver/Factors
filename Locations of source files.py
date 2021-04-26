from __future__ import print_function

# Import the os module, for the os.walk function
import os
 
# Set the directory you want to start from
rootDir = 'C:\\Users\\Arthur\\source\\repos\\Factors'
for dirName, subdirList, fileList in os.walk(rootDir):
    #print('Found directory: %s' % dirName)
    for fname in fileList:
        if  fname.endswith('.cpp') or fname.endswith('.cc') or fname.endswith('.h') :
          print( '%s' % fname, '\t',
                 '%s' % dirName.replace('C:\\Users\\Arthur\\source\\repos\\Factors\\',''))