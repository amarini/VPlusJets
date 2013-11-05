#!/usr/bin/python
import sys,os
from glob import glob 
from subprocess import call
from optparse import OptionParser, make_option


#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------    

parser = OptionParser(option_list=[
    ## doesn't work
    make_option("-d", "--dir",
                action="store", type="string", dest="dir",
                default="",
                help="Directory"
                ),
    make_option("-v","--verbose",action="store_true",default=False,dest="verbose",help="verbose"),
    make_option("-x","--delete",action="store_true",default=False,dest="delete",help="Delete not just print command"),
    make_option("-m","--merge",action="store_true",default=False,dest="merge",help="Do Hadd"),
    ])

(options, args) = parser.parse_args()

fileList=[]
for fileCand in glob(options.dir+"/*.root"):
    	fileList.append(fileCand)

offset=len(options.dir.split('_') )

Database={}
#PATZJetsExpress_100_4_vP4.root
for file in fileList:
	fileParts=file.split('_')
	try:
		Database[int(fileParts[offset])]+=1
	except KeyError:
		Database[int(fileParts[offset])] = 1
		
	if Database[int(fileParts[offset])] > 1 :
		print "rm "+file		
		if options.delete:
			call(["rm","-v",file])

if options.merge:
	dir=options.dir
	if dir[-1]=="/":
		target=dir[:-1]+".root"
	else:
		target=dir+".root"
	cmd=["hadd",target]
	for i in glob(options.dir+"/*.root"):
		cmd.append(i)
	call(cmd)

	if( len(cmd) > 1000): print "POSSIBLE ERROR - to many files"

if options.verbose:	
	print "------DATABASE------"
	for i in Database:
		print "-- " + str(i) + " : " + str(Database[i])
	print "------OFFSET--------"
	print offset
	print "--------DIR---------"
	print options.dir
	print "--------END---------"
	
