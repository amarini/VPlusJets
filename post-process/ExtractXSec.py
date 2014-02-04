import ROOT 
import sys,os,array
from optparse import OptionParser, make_option

ROOT.gROOT.SetBatch()

parser=OptionParser()
parser.add_option("-t","--type",type="string",dest="type",help="type %default",default="INT")

(options,args)=parser.parse_args()


for file in args:
	print "FILE IS --- "+file+" ---"
	f=ROOT.TFile.Open(file)
	h=f.Get("accepted/XSec")
	if( options.type.lower() == 'int' or options.type.lower() == 'all'): print " Int xSec   = " + str(h.GetBinContent(h.FindBin(0))/h.GetBinContent(h.FindBin(1)))
	if( options.type.lower() == 'ext' or options.type.lower() == 'all'): print " Ext xSec   = " + str(h.GetBinContent(h.FindBin(2))/h.GetBinContent(h.FindBin(3)))
	if( options.type.lower() == 'nlo' or options.type.lower() == 'all'): 	print " ExtNLOxSec = " + str(h.GetBinContent(h.FindBin(4))/h.GetBinContent(h.FindBin(5)))
