#!/usr/bin/python


import sys,os,array
import subprocess
from glob import glob
from optparse import OptionParser
import time
import random

#global
preCommand=[]
submitOptions=[]

random.seed()
suffix="."+str(random.randint(0,10000))
print "-- suffix is "+suffix+" --"

def CrabStatus(dirName):
	out=subprocess.call("crab -status -c "+dirName+" > /tmp/amarini/monitor.status"+suffix,shell=True)

def CrabGet(dirName):
	out=subprocess.call("crab -get -c "+dirName+" > /tmp/amarini/monitor.get"+suffix,shell=True)

def Status(dirName):
	Result={}
	#out=subprocess.check_output(["crab","-status","-c",dirName],shell=True)
	#out=subprocess.check_output("crab -status -c "+dirName,shell=True)
	#out,err = p.communicate()
	CrabStatus(dirName)
	return ReadStatus("/tmp/amarini/monitor.status"+suffix)

def FullStatus(dirName):
	CrabStatus(dirName)
	CrabGet(dirName)
	CrabStatus(dirName)
	return ReadStatus("/tmp/amarini/monitor.status"+suffix)
	

def ReadStatus(txtName):
	Result={}
	out=open(txtName,'r')

	for line in out:
		R={}
		parts=line.split()
		try:
			R['jobId']=int(parts[0])
			R['end']=parts[1]
			R['status']=parts[2]
			R['action']=parts[3]
			try:
				R['exeExit']=int(parts[4])
				R['jobExit']=int(parts[5])
				R['host']=parts[6]
			except (ValueError,IndexError):
				R['exeExit']=int(-1)
				R['jobExit']=int(-1)
				try:
					R['host']=parts[4]
				except (ValueError,IndexError):
					R['host']=""
			Result[ R['jobId'] ] = R
		except (ValueError,IndexError): pass
	return Result

def CallSubmit(jobList,dirName,Resubmit=False):
	cmd=[]+preCommand
	if Resubmit:
		cmd += ['crab','-resubmit']
	else:
		cmd += ['crab','-submit']
	jobs=''
	for iJob in jobList:
		jobs=jobs+str(iJob)+','
	#delete last
	jobs=jobs[:-1]
	cmd=cmd+[jobs,"-c",dirName]
	cmd=cmd+submitOptions
	subprocess.call(cmd)

def Submit(jobList,dirName,Resubmit=False):
	for i in range(0,len(jobList),450):
		CallSubmit(jobList[i:i+450],dirName,Resubmit);

def SubmitAll(Result,dirName,Resubmit=True):
	jobList=[]
	for key in Result:
		jobList.append(key)
	Submit(jobList,dirName,Resubmit)	

def ReSubmitES(exitN,Result,dirName):
	jobList=[]
	for key in Result:
		if Result[key]['exeExit']==exitN or Result[key]['jobExit']==exitN:
			jobList.append(key)
	Submit(jobList,dirName,True)

def ReSubmitCA(Result,dirName):
	jobList=[]
	for key in Result:
		if Result[key]['status']=="Cancelled" or Result[key]['status']=='Aborted':
			jobList.append(key)
	Submit(jobList,dirName,True)

def PrintStatus(Result,dirName=""):
	DONE=0
	NOTDONE=0
	Cancelled=0
	Aborted=0
	Submitted=0
	Running=0
	Exit={}
	for key in Result:
		try:
			Exit[ Result[key]['jobExit'] ] += 1
		except KeyError:
			Exit[ Result[key]['jobExit'] ] = 1
		if Result[key]['status'] =="Cancelled":
			Cancelled += 1
		if Result[key]['status'] =="Aborted":
			Aborted += 1
		if Result[key]['status'] =="Submitted":
			Submitted += 1
		if Result[key]['status'] =="Running":
			Running += 1
		if Result[key]['end'] == "Y" and Result[key]['jobExit'] ==0:
			DONE += 1
		else:
			NOTDONE += 1
	
	if DONE+NOTDONE == 0: NOTDONE=-1
	
	print "---------------------------"
	print "| " +dirName + " |"
	print "---------------------------"
	print "\033[01;32mDONE    = %4d     (%.2f%%)"%(DONE,float(DONE)/(DONE+NOTDONE)*100)
	print "\033[01;31mNOTDONE = %4d     (%.2f%%)"%(NOTDONE,float(NOTDONE)/(DONE+NOTDONE)*100)
	print "\033[0m---------------------------"
	line="\033[01;34m["
	tot=25
	dash=int(tot*float(DONE)/(DONE+NOTDONE) )
	space=tot-dash-1
	for i in range(0,dash):line += '='
	if dash<tot: line+='>'
	for i in range(0,space):line+= ' '
	line +=']\033[0m'
	print line
	print "---------------------------"
	print "Cancelled: %d"%(Cancelled)
	print "Aborted: %d"%(Aborted)
	print "Submitted: %d"%(Submitted)
	print "Running: %d"%(Running)
	for key in Exit:
		print "ExitStatus %d: %d"%(key,Exit[key])
	print "---------------------------"

def PrintDatabase(Results):
	for key in Results:
		print str(key)+":"
		print Results[key]


if __name__ == "__main__":
	parser=OptionParser(usage=" %prog [options] dir1 dir2 ...\nProvide help to run/submit/resubmit crab jobs")
	parser.add_option("-l","--loop",help="Run in loop",dest="loop",default=False,action='store_true')
	parser.add_option("-n","--dryrun",help="Dry Run: Print command instead of re/submit",dest="dryrun",default=False,action='store_true')
	parser.add_option("-s","--submit",help="Submit jobs",dest='submit',default=False,action='store_true')
	parser.add_option("-x","--extra",help="Extra Options for submit. E.g. '-GRID.se_black_list=T2;-GRID.ce_black..'",dest='extra',default='',type='string')
	parser.add_option("-e","--exitNumber",help="Resubmit Exit 1,2,3 ...\nDefault=8020,8021,60317,60307,50664,50115,10034,8001,8022",dest='exit',default='8020,8021,60317,60307,50664,50115,10034,8001,8022,8028',type='string') 
	(options,arg) = parser.parse_args()

	if options.dryrun:
		preCommand=["echo"]

	if options.extra != '':
		submitOptions=options.extra.split(';')
	
	if options.exit != "":
		numList=[]
		for i in options.exit.split(','):
			numList.append(int(i))
	else:
		numList=[]
	
	#main loop
	while True:		
		for dir in arg:
			#CrabStatus(dir)
			#Results=ReadStatus("/tmp/amarini/status")
			Results=FullStatus(dir)
			#PrintDatabase(Results)
			PrintStatus(Results,dir)
			if options.submit:
				SubmitAll(dir)
			print "-- Resubmit Cancelled & Aborted --"		
			ReSubmitCA(Results,dir)
			for exitStatus in numList:
				print "-- Resubmit Exit "+str(exitStatus)+" --"
				ReSubmitES(exitStatus,Results,dir)
		if options.loop : time.sleep(15*60)
		else: break;

