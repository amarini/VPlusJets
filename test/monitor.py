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

user=os.environ['USER']

import signal
def signal_handler(signal, frame):
        print 'You pressed Ctrl+C!'
        sys.exit(0)

def CrabStatus(dirName):
	out=subprocess.call("crab -status -c "+dirName+" > /tmp/"+user+"/monitor.status"+suffix,shell=True)

def CrabGet(dirName):
	out=subprocess.call("crab -get -c "+dirName+" > /tmp/"+user+"/monitor.get"+suffix,shell=True)

def Status(dirName):
	Result={}
	#out=subprocess.check_output(["crab","-status","-c",dirName],shell=True)
	#out=subprocess.check_output("crab -status -c "+dirName,shell=True)
	#out,err = p.communicate()
	CrabStatus(dirName)
	return ReadStatus("/tmp/"+user+"/monitor.status"+suffix)

def FullStatus(dirName):
	CrabStatus(dirName)
	CrabGet(dirName)
	CrabStatus(dirName)
	return ReadStatus("/tmp/"+user+"/monitor.status"+suffix)
	

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

def CallForceResubmit(jobList,dirName):
	cmd=[]+preCommand
	cmd += ['crab','-forceResubmit']
	jobs=''
	for iJob in jobList:
		jobs=jobs+str(iJob)+','
	#delete last
	jobs=jobs[:-1]
	cmd=cmd+[jobs,"-c",dirName]
	cmd=cmd+submitOptions
	subprocess.call(cmd)
	print " -- GOING TO SLEEP FOR 60s --"
	time.sleep(60)
	print " -- WAKE UP --"

def ForceResubmit(jobList,dirName):
	for i in range(0,len(jobList),450):
		CallForceResubmit(jobList[i:i+450],dirName);

def SubmitAll(Result,dirName,Resubmit=False):
	jobList=[]
	for key in Result:
		jobList.append(key)
	Submit(jobList,dirName,Resubmit)	

def ReSubmitES(exitN,Result,dirName):
	jobList=[]
	sites={}
	for key in Result:
		if Result[key]['exeExit']==exitN or Result[key]['jobExit']==exitN:
			jobList.append(key)
			if exitN == 8020: ## site error: ban T2
				SE=subprocess.check_output("cat "+dirName+"/res/CMSSW_"+str(key)+".stdout | grep '^SyncSite='  | tail -n 1 | cut -d '=' -f2", shell=True)	
				try:
					sites[ SE ] += 1
				except KeyError:
					try: sites[ SE ] = 1
					except KeyError: pass
	for se in sites:
		if sites[ se ] > 10 and len(se) >4 : #avoid malformed strings
			print " BANNING SE " + se
			found=False
			for i in range(0,len(submitOptions)):
				if "se_black_list" in submitOptions[i]:
					if not se in submitOptions[i]:submitOptions[i]+=","+se
					found=True
			if not Found:
				submitOptions.append("-GRID.se_black_list="+se)
			
	Submit(jobList,dirName,True)

def ReSubmitCA(Result,dirName):
	jobList=[]
	for key in Result:
		if Result[key]['status']=="Cancelled" or Result[key]['status']=='Aborted':
			jobList.append(key)
	Submit(jobList,dirName,True)

def ForceReSubmitSubmitted(Result,dirName):
	jobList=[]
	for key in Result:
		if Result[key]['status']=="Submitted":
			jobList.append(key)
	ForceResubmit(jobList,dirName)

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
	parser.add_option("","--forceResubmitSubmittedJobs",help="Force Resubmition of jobs in status submitted",dest="forceResubmit",default=False,action='store_true')
	(options,arg) = parser.parse_args()

	if options.dryrun:
		preCommand=["echo"]
	
	if options.forceResubmit and options.loop:
		print "---> CANT LOOP ON FORCE RESUBMIT <---"
		sys.exit(0)
	
	if options.forceResubmit:
		pwd=random.randint(0,10000)
		line=raw_input('You are going to force resubmit jobs: Enter %d to confirm: '%(pwd))
		if int(line) != pwd:
			print "you have entered '"+line+"' instead of '"+str(pwd)+"': Going to abort execution."
			sys.exit(0)

	if options.extra != '':
		submitOptions=options.extra.split(';')
	
	if options.exit != "":
		numList=[]
		for i in options.exit.split(','):
			numList.append(int(i))
	else:
		numList=[]
	
	#main loop
	if options.loop:
		signal.signal(signal.SIGINT, signal_handler)
		
	if options.loop and options.submit:
		print "-- cannot be loop & submit --"
		sys.exit(0)

	while True:		
		for dir in arg:
			if options.submit:
				Results=Status(dir)
				SubmitAll(Results,dir)
			Results=FullStatus(dir)
			PrintStatus(Results,dir)
			print "-- Resubmit Cancelled & Aborted --"		
			ReSubmitCA(Results,dir)
			for exitStatus in numList:
				print "-- Resubmit Exit "+str(exitStatus)+" --"
				ReSubmitES(exitStatus,Results,dir)

			if options.forceResubmit:
				ForceReSubmitSubmitted(Results,dir)
				sys.exit(0)

		if options.loop : 
			print "--- GOING TO SLEEP ----"
			time.sleep(15*60)
			print "--- WAKE UP ----"
		else: break;

