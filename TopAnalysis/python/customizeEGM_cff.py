import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# EGM corrections : https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#Recipe_for_running_Scales_and_sm

def customizeEGM(process,era,runOnAOD=False):

    if '2016' in era: 
        egmEra='2016-Legacy'
        runEnergyCorrections=True
    if '2017' in era:
        egmEra='2017-UL'   
        runEnergyCorrections=True
    if '2018' in era: 
        egmEra="2018-Prompt"
        runEnergyCorrections=True

    setupEgammaPostRecoSeq(process,
                           isMiniAOD=(not runOnAOD),
                           era=egmEra,
                           runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           runEnergyCorrections=runEnergyCorrections)

    if runOnAOD:
        print 'Adapting e/g sources to AOD'
        process.electronMVAValueMapProducer.src = cms.InputTag("")
        process.photonMVAValueMapProducer.src = cms.InputTag("")
        #process.photonIDValueMapProducer.src = cms.InputTag("")

    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)