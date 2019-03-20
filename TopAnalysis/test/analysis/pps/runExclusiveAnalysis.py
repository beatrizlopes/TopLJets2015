#! /usr/bin/env python
import ROOT
import copy
import optparse
import json
import sys
import os
import numpy
import array
import json
import random
import pickle
import array
import re
from random import shuffle
from collections import OrderedDict
from TopLJets2015.TopAnalysis.HistoTool import *
from mixedEvent import *

VALIDLHCXANGLES=[120,130,140,150]

def isValidRunLumi(run,lumi,runLumiList):

    """checks if run is available and lumi section was certified"""

    #no run/lumi to select, all is good by default
    if not runLumiList:
        return True

    #run is available
    if not run in runLumiList: 
        return True

    for lran in runLumiList[run]:
        if lumi>=lran[0] and lumi<=lran[1]:
            return False

    #reached this far, nothing found in list
    return True

def getTracksPerRomanPot(tree,mcTruth=False,usefarRP=True):

    """loops over the availabe tracks in the event and groups them by roman pot id"""

    tkPos=[]
    tkNeg=[]
    for itk in xrange(0,tree.nRPtk):
        rpid=tree.RPid[itk]
        try:
            csi=tree.RPfarcsi[itk] if usefarRP else tree.RPnearcsi[itk]
        except:            
            if notMCTruth:
                csi=tree.RPfarx[itk] if useFarRP else tree.RPnearx[itk]
            else:
                csi=tree.RPtruecsi[itk] 
        if rpid==23  and csi>0: tkPos.append(csi)
        if rpid==123 and csi>0: tkNeg.append(csi)
        
    return (tkPos,tkNeg)

def buildDiproton(rptks,sqrts=13000.):

    """build a diproton system from to tracks in the arms of the roman pots"""

    if not rptks : return None
    if len(rptks[0])==0 or len(rptks[1])==0 : return None
    beamP=0.5*sqrts
    csiL,csiR=rptks[0][0],rptks[1][0]
    diproton=ROOT.TLorentzVector(0.,0.,beamP*(csiL-csiR),beamP*(csiL+csiR))
    return diproton

def getRandomEra():

    """generates a random era according to the integrated luminosity in each one"""

    r=random.random()
    if r<0.115   : return '2017B'
    elif r<0.348 : return '2017C'
    elif r<0.451 : return '2017D'
    elif r<0.671 : return '2017E'
    return '2017F'

def runExclusiveAnalysis(inFile,outFileName,runLumiList,mixFile):
    
    """event loop"""

    isData=True if 'Data' in inFile else False
    isSignal=True if 'MC13TeV_2017_ppxz_' in inFile else False

    #bind main tree with pileup discrimination tree, if failed return
    tree=ROOT.TChain('analysis/data' if isSignal else 'tree')
    tree.AddFile(inFile)
    #try:
    #    pudiscr_tree=ROOT.TChain('pudiscr')
    #    baseName=os.path.basename(inFile)
    #    baseDir=os.path.dirname(inFile)
    #    pudiscr_file=os.path.join(baseDir,'pudiscr',baseName)
    #    if not os.path.isfile(pudiscr_file):
    #        raise ValueError(pudiscr_file+' does not exist')
    #    pudiscr_tree.AddFile(pudiscr_file)
    #    tree.AddFriend(pudiscr_tree)
    #    print 'Added pu tree for',inFile
    #except:
    #    #print 'Failed to add pu discrimination tree as friend for',inFile
    #    return


    #identify data-taking era
    era=None
    if isData:
        era=os.path.basename(inFile).split('_')[1]
    
    #check if it is signal and load     
    signalPt=[]
    mcEff={}
    if isSignal: 
        signalPt=[float(x) for x in re.findall(r'\d+', os.path.basename(inFile) )[2:]]
        for ch in ['ee','mm']:
            effIn=ROOT.TFile.Open('$CMSSW_BASE/src/TopLJets2015/TopAnalysis/plots/effsummary_%sz_ptboson.root'%ch)
            mcEff[ch]=effIn.Get('gen%sz2trec_ptboson_ZH#rightarrowllbb_eff'%ch)
            effIn.Close()

    #filter events to mix according to tag if needed
    mixedRP=None
    xangleRelFracs={}
    try:
        print 'Analysing mixfile',mixFile
        with open(mixFile,'r') as cachefile:
            mixedRP=pickle.load(cachefile)        

        #build the list of probabilities for the crossing angles in each era
        for key in mixedRP:
            era,xangle=key
            if not xangle in VALIDLHCXANGLES: continue
            n=len(mixedRP[key])
            if not era in xangleRelFracs:
                xangleRelFracs[era]=ROOT.TH1F('xanglefrac_%s'%era,'',len(VALIDLHCXANGLES),0,len(VALIDLHCXANGLES))
            xbin=(xangle-120)/10+1
            xangleRelFracs[era].SetBinContent(xbin,n)
        for era in xangleRelFracs:
            xangleRelFracs[era].Scale(1./xangleRelFracs[era].Integral())
       
    except Exception as e:
        if mixFile : print e
        pass
    
    #start histograms
    ht=HistoTool()   

    #main analysis histograms
    ht.add(ROOT.TH1F('mll',';Dilepton invariant mass [GeV];Events',50,20,250))
    ht.add(ROOT.TH1F('yll',';Dilepton rapidity;Events',50,0,5))
    ht.add(ROOT.TH1F('ptll',';Dilepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l1eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l1pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l2eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l2pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('acopl',';A=1-|#Delta#phi|/#pi;Events',50,0,1))
    ht.add(ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',4,120,160))
    ht.add(ROOT.TH1F('mpp',';Di-proton invariant mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ypp',';Di-proton rapidity;Events',50,0,2))
    ht.add(ROOT.TH2F('mpp2d',';Far di-proton invariant mass [GeV];Near di-proton invariant mass [GeV];Events',50,0,3000,50,0,3000))
    ht.add(ROOT.TH2F('ypp2d',';Far di-proton rapidity;Near di-proton rapidity;Events',50,0,2,50,0,2))
    ht.add(ROOT.TH1F('mmass',';Missing mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('ppcount',';pp candidates;Events',3,0,3))
    ht.add(ROOT.TH1F('csi',';#xi;Events',50,0,0.3))
    ht.add(ROOT.TH2F('csi2d',';#xi(far);#xhi(near);Events',50,0,0.3,50,0,0.3))

    #pileup control
    ht.add(ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100))
    ht.add(ROOT.TH1F('rho',';Fastjet #rho;Events',50,0,30))
    #ht.add(ROOT.TH1F('rfc',';Random forest classifier probability;Events',50,0,1))
    for d in ['HF','HE','EE','EB']:
        ht.add(ROOT.TH1F(d+'PFMult',';PF multiplicity (%s);Events'%d,50,0,1000))
        ht.add(ROOT.TH1F(d+'PFHt',';PF HT (%s) [GeV];Events'%d,50,0,1000))
        ht.add(ROOT.TH1F(d+'PFHt',';PF P_{z} (%s) [TeV];Events'%d,50,0,40))
    ht.add(ROOT.TH1F('met',';Missing transverse energy [GeV];Events',50,0,200))
    ht.add(ROOT.TH1F('metbits',';MET filters;Events',124,0,124))
    ht.add(ROOT.TH1F('njets',';Jet multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('nch', ';Charged particle multiplicity;Events',50,0,50))
    ht.add(ROOT.TH1F('nextramu',';Additional muons ;Events',10,0,10))
    ht.add(ROOT.TH1F('extramupt',';Additional muon p_{T} [GeV] ;Events',10,0,50))
    ht.add(ROOT.TH1F('extramueta',';Additional muon pseudo-rapidty ;Events',10,0,2.5))
        
    nEntries=tree.GetEntries()          
    print '....analysing',nEntries,'in',inFile,', with output @',outFileName
    if mixedRP : print '    events mixed with',mixFile

    #loop over events
    rpData={}
    selEvents=[]
    summaryVars='cat:wgt:nvtx:rho:PFMultSumHF:PFHtSumHF:PFPFPzSumHF:nch:xangle:l1pt:l1eta:l2pt:l2eta:acopl:bosonpt:csi1:csi2:mpp:mmiss:mixcsi1:mixcsi2:mixmpp:mixmmiss'.split(':')
    for i in xrange(0,nEntries):

        tree.GetEntry(i)
        
        if i%1000==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
    
        #base event selection
        if tree.evcat==11*11   : evcat='ee'
        elif tree.evcat==11*13 : evcat='em'
        elif tree.evcat==13*13 : evcat='mm'
        elif tree.evcat==22 and tree.hasLowPtATrigger  and tree.bosonpt<200  : evcat="lpta"
        elif tree.evcat==22 and tree.hasHighPtATrigger and tree.bosonpt>=200 : evcat="hpta"
        elif tree.evcat==0  and tree.hasZBTrigger : evcat=='zbias'
        else : continue

        #filter same-sign events
        if tree.evcat!=22 and tree.evcat!=0 and tree.isSS : continue        

        #check RPs are in
        if isData and not isValidRunLumi(tree.run,tree.lumi,runLumiList): continue

        wgt=tree.evwgt
        nvtx=tree.nvtx
        nch=tree.nchPV
        rho=tree.rho
        met=tree.met_pt
        njets=0 if isSignal else tree.nj 

        #lepton specific
        l1p4=ROOT.TLorentzVector(0,0,0,0)
        l2p4=ROOT.TLorentzVector(0,0,0,0)        
        acopl=0
        if tree.evcat!=22 and tree.evcat!=0:
            acopl=1.0-abs(ROOT.TVector2.Phi_mpi_pi(tree.l1phi-tree.l2phi))/ROOT.TMath.Pi()
            l1p4.SetPtEtaPhiM(tree.l1pt,tree.l1eta,tree.l1phi,tree.ml1)
            l2p4.SetPtEtaPhiM(tree.l2pt,tree.l2eta,tree.l2phi,tree.ml2)
            if l1p4.Pt()<l2p4.Pt(): l1p4,l2p4=l2p4,l1p4
            if abs(l1p4.Eta())>2.1 : continue

        #boson kinematics
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        isZ=tree.isZ
        isA=tree.isA
        isHighPt=(boson.Pt()>50)
        isLowPt=(boson.Pt()<10)

        #Nicola's initial discriminator
        extra_muons=[]
        for im in range(tree.nrawmu):
            mup4=ROOT.TLorentzVector(0,0,0,0)
            mup4.SetPtEtaPhiM(tree.rawmu_pt[im],tree.rawmu_eta[im]/10.,tree.rawmu_phi[im]/10.,0.105)
            if mup4.DeltaR(l1p4)<0.05 : continue
            if mup4.DeltaR(l2p4)<0.05 : continue
            extra_muons.append( ROOT.TLorentzVector(mup4) )
            
        #proton tracks (standard and mixed)
        far_rptks, near_rptks = None, None
        beamXangle     = int(tree.beamXangle)
        if isSignal or isData: 
            far_rptks  = getTracksPerRomanPot(tree)
            near_rptks = getTracksPerRomanPot(tree,False,False)
            if not beamXangle in VALIDLHCXANGLES : continue

        mixed_far_rptks  = None
        mixed_far_1rptk  = None
        if not isData:
            evEra=getRandomEra()
            if not isSignal: 
                xbin=ROOT.TMath.FloorNint(xangleRelFracs[era].GetRandom())
                beamXangle=VALIDLHCXANGLES[xbin]
        else:
            evEra=era
            if not mixedRP:
                #if isZ and tree.evcat==13*13 and tree.bosonpt<10 and beamXangle in VALIDLHCXANGLES:
                if evcat=='em':
                    rfc=0 #getattr(tree,'rfc_%d'%beamXangle)
                    rpDataKey=(era,beamXangle)
                    if not rpDataKey in rpData: rpData[rpDataKey]=[]
                    rpData[rpDataKey].append( MixedEvent(beamXangle,
                                                         [len(extra_muons),nvtx,rho,rfc],
                                                         far_rptks,
                                                         near_rptks
                                                         ) )
                continue

        try:
            mixedEv=random.choice( mixedRP[(evEra,beamXangle)] )
            mixed_far_rptks=mixedEv.far_rptks
            if far_rptks and mixed_far_rptks:
                if random.random()<0.5:
                    mixed_far_1rptk = (mixed_far_rptks[0],far_rptks[1])
                else:
                    mixed_far_1rptk = (far_rptks[0],mixed_far_rptks[1])
        except Exception as e:
            print e
            pass

        #build combinations of histograms
        goldenSel=None
        mon_tasks=[]
        if isData:
            mon_tasks.append( (far_rptks,near_rptks,'') )
            mon_tasks.append( (mixed_far_1rptk,None,'_mix1') )
            mon_tasks.append( (mixed_far_rptks,None,'_mix2') )
        else:
            if isSignal:
                #embed pileup to signal
                tksPos=mixed_far_rptks[0]+far_rptks[0]
                shuffle(tksPos)
                tksNeg=mixed_far_rptks[1]+far_rptks[1]
                shuffle(tksNeg)
                mixed_far_rptks=(tksPos,tksNeg)                
                mon_tasks.append( (far_rptks,near_rptks,'nopu') )
            mon_tasks.append( (mixed_far_rptks,None,'') )

        for protons,near_protons, pfix in mon_tasks:

            #no calibration, not worth it...
            if not beamXangle in VALIDLHCXANGLES: continue
    
            #high purity selection for proton tracks
            highPur = True if protons and len(protons[0])==1 and len(protons[1])==1 else False  
            if highPur:
                if protons[0][0]<0.015 or protons[0][0]>0.20 : highPur=False
                if protons[1][0]<0.015 or protons[1][0]>0.20 : highPur=False

            #check if Z and di-proton combination is consistent with elastic scattering
            pp=buildDiproton(protons)
            near_pp=buildDiproton(near_protons) if near_protons else None
            isElasticLike=False
            mmass=0
            if pp:
                isElasticLike=(13000.-boson.E()-pp.E()>0)
                inPP=ROOT.TLorentzVector(0,0,0,13000.)
                if isElasticLike: 
                    mmass=(pp-boson).M()

            #categories to fill
            cats=[]            
            cats.append(evcat)
            if isZ : 
                cats.append(evcat+'Z')
                if isHighPt : cats.append(evcat+'hptZ')
                if isLowPt  : cats.append(evcat+'lptZ')
            if isElasticLike and highPur :
                ppCats=[c+'hpur' for c in cats]
                cats += ppCats
            beamAngleCats=[c+'%d'%beamXangle for c in cats]
            cats += beamAngleCats
                
            if (isData and 'mix' in pfix) or (isSignal and pfix==''):
                if (isZ or isA) and isElasticLike and highPur:
                    if goldenSel:
                        goldenSel += [protons[0][0],protons[1][0],pp.M(),mmass]
                    else:
                        goldenSel=[tree.evcat,wgt,
                                   nvtx,tree.rho,tree.PFMultSumHF,tree.PFHtSumHF,tree.PFPzSumHF,0,
                                   nch,beamXangle,l1p4.Pt(),l1p4.Eta(),l2p4.Pt(),l2p4.Eta(),acopl,boson.Pt(),
                                   protons[0][0],protons[1][0],pp.M(),mmass]

            #final plots (for signal correct wgt by efficiency curve and duplicate for mm channel)
            finalPlots=[[wgt,cats]]
            if isSignal:
                #signal has been pre-selected in the fiducial phase space so nEntries is 
                #the effective number of events (or sum of weights) generated
                finalPlots=[ [wgt*mcEff['ee'].Eval(boson.Pt())/nEntries, cats],
                             [wgt*mcEff['mm'].Eval(boson.Pt())/nEntries, [c.replace(evcat,'mm') for c in cats if c[0:2]=='ee']] ]

            for pwgt,pcats in finalPlots:
                
                ht.fill((nvtx,pwgt),                  'nvtx',  pcats,pfix)
                ht.fill((rho,pwgt),                   'rho',   pcats,pfix)
                ht.fill((met,pwgt),                   'met',   pcats,pfix)
                ht.fill((tree.metfilters,pwgt),       'metbits',   pcats,pfix)
                ht.fill((njets,pwgt),                 'njets',  pcats,pfix)
                ht.fill((nch,pwgt),                   'nch',    pcats,pfix)                                
                ht.fill((l1p4.Pt(),pwgt),             'l1pt',   pcats,pfix)
                ht.fill((l2p4.Pt(),pwgt),             'l2pt',   pcats,pfix)
                ht.fill((abs(l1p4.Eta()),pwgt),       'l1eta',  pcats,pfix)
                ht.fill((abs(l2p4.Eta()),pwgt),       'l2eta',  pcats,pfix)
                ht.fill((acopl,pwgt),                 'acopl',  pcats,pfix) 
                ht.fill((boson.M(),pwgt),             'mll',    pcats,pfix)
                ht.fill((abs(boson.Rapidity()),pwgt), 'yll',    pcats,pfix)
                ht.fill((boson.Pt(),pwgt),            'ptll',   pcats,pfix)
                ht.fill((beamXangle,pwgt),                'xangle', pcats,pfix)
                #ht.fill((getattr(tree,'rfc_%d'%beamXangle),pwgt), 'rfc',         pcats,pfix)
                for sd in ['HF','HE','EE','EB']:
                    ht.fill((getattr(tree,'PFMultSum'+sd),pwgt),     'PFMult'+sd,    pcats,pfix)
                    ht.fill((getattr(tree,'PFHtSum'+sd),pwgt),       'PFHt'+sd,      pcats,pfix)
                    ht.fill((getattr(tree,'PFPzSum'+sd)/1.e3,pwgt),  'PFPZ'+sd,      pcats,pfix)
                
                ht.fill((len(extra_muons),pwgt), 'nextramu', pcats, pfix)
                for mp4 in extra_muons:
                    ht.fill((mp4.Pt(),pwgt), 'extramupt', pcats,pfix)
                    ht.fill((abs(mp4.Eta()),pwgt), 'extramueta', pcats,pfix)

                for irp,rpside in [(0,'pos'),(1,'neg')]:
                    ht.fill((len(protons[irp]),pwgt), 'ntk',    pcats,rpside+pfix)

                ht.fill((0,pwgt), 'ppcount', pcats,pfix)
                if not pp: continue        

                ht.fill((1,pwgt), 'ppcount', pcats,pfix)
                ht.fill((pp.M(),pwgt),             'mpp',       pcats,pfix)
                ht.fill((abs(pp.Rapidity()),pwgt), 'ypp',       pcats,pfix)
                
                if near_pp:
                    ht.fill((2,pwgt), 'ppcount', pcats,pfix)
                    ht.fill((pp.M(),near_pp.M(),pwgt),                         'mpp2d',  pcats,pfix)
                    ht.fill((abs(pp.Rapidity()),abs(near_pp.Rapidity()),pwgt), 'ypp2d',  pcats,pfix)

                ht.fill((mmass,pwgt),              'mmass',     pcats,pfix)
                for irp,rpside in [(0,'pos'),(1,'neg')]:

                    ht.fill((len(protons[irp]),pwgt), 'ntk',    pcats,rpside+pfix)

                    for csi in protons[irp]:
                        ht.fill((csi,pwgt), 'csi',             pcats,rpside+pfix)
                        
                        if not near_protons : continue
                        if len(near_protons[irp])==0 : continue
                        csi_near=near_protons[irp][0]
                        ht.fill((csi,csi_near,pwgt), 'csi2d', pcats,rpside+pfix)
                    



        #select events
        if goldenSel:
            nVarsMissed=len(summaryVars)-len(goldenSel)
            if nVarsMissed>0: goldenSel += [0.]*nVarsMissed

            if isSignal:
                                
                #add a copy for ee
                eeGoldenSel=copy.copy(goldenSel)
                eeGoldenSel[0]=11*11
                eeGoldenSel[1]=goldenSel[1]*mcEff['ee'].Eval(boson.Pt())/nEntries
                selEvents.append(eeGoldenSel)
                
                #add a copy for mm
                mmGoldenSel=copy.copy(goldenSel)
                mmGoldenSel[0]=13*13
                mmGoldenSel[1]=goldenSel[1]*mcEff['mm'].Eval(boson.Pt())/nEntries
                selEvents.append(mmGoldenSel)

            else:

                selEvents.append(goldenSel)
            

    #dump events for the mixing
    nSelRPData=sum([len(rpData[x]) for x in rpData])
    if nSelRPData>0:
        rpDataOut=outFileName.replace('.root','.pck')        
        print 'Saving',nSelRPData,'events for mixing in',rpDataOut
        with open(rpDataOut,'w') as cachefile:
            pickle.dump(rpData,cachefile, pickle.HIGHEST_PROTOCOL)        

        #at this point no need to save anything else
        return

    #save results
    ht.writeToFile(outFileName)

    #dump events for fitting
    nSelEvents=len(selEvents)
    if nSelEvents>0:
        print 'Adding',nSelEvents,'selected events to',outFileName
        fOut=ROOT.TFile.Open(outFileName,'UPDATE')
        fOut.cd()
        t=ROOT.TNtuple('data','data',':'.join(summaryVars))
        for v in selEvents : t.Fill(array.array("f",v))
        t.Write()
        fOut.Close()
 
def runExclusiveAnalysisPacked(args):

    """wrapper for parallel execution"""

    try:
        runExclusiveAnalysis(*args)
    except Exception as e:
        print 50*'<'
        print "  Problem with", args[1], "continuing without"
        print e
        print 50*'<'
        return False
    
def runAnalysisTasks(opt):

    """create analysis tasks"""

    #read samples to process
    task_dict={}    
    with open(opt.json,'r') as cachefile:
        samples=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        for x in samples:
            task_dict[x[0]]=[]

    #group chunks matching the same name
    for file_path in os.listdir(opt.input):
        
        if opt.only and file_path!=os.path.basename(opt.only): continue

        file_name,ext=os.path.splitext(file_path)
        if ext != '.root' : continue
        
        #check if file tag is already in the list of samples to process
        lastTkn=file_name.rfind('_')
        tag=file_name[0:lastTkn]

        if not tag in task_dict: continue
        task_dict[tag].append( os.path.join(opt.input,file_path) )

    #parse json file with list of run/lumi sections
    with open(opt.RPout,'r') as cachefile:        
        runLumi=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        runLumi={int(x[0]):x[1] for x in runLumi}

    #create the tasks and submit them
    import multiprocessing as MP
    pool = MP.Pool(opt.jobs)
    task_list=[]
    for x in task_dict.keys():
        
        isData = True if 'Data' in x else False
        if opt.step==0 and not isData : continue
        runLumiList=runLumi if isData else None
        for f in task_dict[x]:
            fOut='%s/Chunks/%s'%(opt.output,os.path.basename(f))
            task_list.append( (f,fOut,runLumiList,opt.mix) )

    pool.map(runExclusiveAnalysisPacked,task_list)


def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/Chunks',
                      help='input directory with the files [default: %default]')
    parser.add_option('--jobs',
                      dest='jobs', 
                      default=8,
                      type=int,
                      help='# of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--json',
                      dest='json', 
                      default='pps_samples.json',
                      type='string',
                      help='json with the files to process')
    parser.add_option('--only',
                      dest='only', 
                      default=None,
                      type='string',
                      help='only this file')
    parser.add_option('--RPout',
                      dest='RPout', 
                      default='golden_noRP.json',
                      type='string',
                      help='json with the runs/lumi sections in which RP are out')
    parser.add_option('--step',
                      dest='step', 
                      default=0,
                      type=int,
                      help='analysis step: 0 - prepare event mixing bank; 1 - analysis')
    parser.add_option('--mix',
                      dest='mix',
                      default=None,
                      type='string',
                      help='bank of events to use for the mixing')
    parser.add_option('-o', '--output',
                      dest='output', 
                      default='analysis',
                      help='Output directory [default: %default]')
    (opt, args) = parser.parse_args()
    
    os.system('mkdir -p %s/Chunks' % opt.output)
    runAnalysisTasks(opt)
        

if __name__ == "__main__":
    sys.exit(main())
