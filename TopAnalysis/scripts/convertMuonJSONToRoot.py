#!/usr/bin/env python

import sys, math

def getRange(tkn):
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+",tkn)]

if len(sys.argv)<2:
  print('ERROR: missing input file.')
  quit()
  

import json
import ROOT
import re
from array import array

url=sys.argv[1]
fOut=ROOT.TFile.Open(url.replace('.json','.root'),'RECREATE')

with open(url,'r') as f:
    results = json.load(f)
    params='abseta_pt'
    for key in results.keys():

        effVals=[]
        for etaTkn in results[key][params]:
            if 'binning' in etaTkn: continue
            ieffVals=[]
            for ptTkn in results[key][params][etaTkn]:
                data=results[key][params][etaTkn][ptTkn]			
                value=data["value"]
                error=0.
                for _k in ['error','stat','syst']:
                  if _k in data.keys(): error+=data[_k]**2
                error=math.sqrt(error);
                ieffVals.append( getRange(ptTkn)+[value,error] )
            ieffVals.sort(key=lambda x : x[0])
            effVals.append( (getRange(etaTkn),ieffVals) )

        effVals.sort(key=lambda x : x[0][0])

        etaVals=[]
        for i in xrange(0,len(effVals)):
            etaVals.append( effVals[i][0][0] )
        etaVals.append(effVals[-1][0][1])

        ptVals=[]
        for i in xrange(0,len(effVals[0][1])):
            ptVals.append( effVals[0][1][i][0] )
        ptVals.append( effVals[0][1][-1][1] )

        fOut.cd()
        sfH=ROOT.TH2F(key+'_'+params,key,len(etaVals)-1,array('d',etaVals),len(ptVals)-1,array('d',ptVals))
        for i in xrange(0,len(effVals)):
            for j in xrange(0,len(effVals[i][1])):
                sf,sfUnc=effVals[i][1][j][2:4]
                sfH.SetBinContent(i+1,j+1,sf)
                sfH.SetBinError(i+1,j+1,sfUnc)
        sfH.Write()

print('Writes %s'%fOut.GetName())
fOut.Close()
        


