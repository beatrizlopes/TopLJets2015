import FWCore.ParameterSet.Config as cms


def setupProtonSim(process,xangle,withPU=False):

  # update settings of beam-smearing module
  process.beamDivergenceVtxGenerator.src = cms.InputTag("")
  
  # input collections
  if withPU:
    process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(
      #cms.InputTag("genPUProtons","genPUProtons"), # works with step2_premix modifier
      cms.InputTag("genPUProtons",""),
      #cms.InputTag("prunedGenParticles"), # when ~premix_stage2 signal protons proporate to genPUProtons
    )
  else:
    process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(cms.InputTag("prunedGenParticles"))
  
  # do not apply vertex smearing again
  process.ctppsBeamParametersESSource.vtxStddevX = 0
  process.ctppsBeamParametersESSource.vtxStddevY = 0
  process.ctppsBeamParametersESSource.vtxStddevZ = 0
  
  #undo CMS vertex shift (example)
  process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
  process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
  process.ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1
  
  # set up Xangle  
  process.ctppsLHCInfoESSource.xangle = cms.double(xangle)
  process.ctppsBeamParametersESSource.halfXangleX45 = xangle * 1E-6
  process.ctppsBeamParametersESSource.halfXangleX56 = xangle * 1E-6
  
  # for multiRP fit, set if you want to use x* and y* as free parameters or set them to zero
  process.ctppsProtons.fitVtxY = False

  # processing path 
  process.pps_fastsim = cms.Path(process.beamDivergenceVtxGenerator
    * process.ctppsDirectProtonSimulation
    * process.reco_local
    * process.ctppsProtons
 )

def ctppsCustom(process):

  # define conditions
  SetConditions(process)
  
def setupProtonReco(process, reMiniAOD = False):
  
  if reMiniAOD:
    process.load("CondCore.CondDB.CondDB_cfi")
    from CondCore.CondDB.CondDB_cfi import CondDB
    # override alignment
    process.CondDBAlignment = CondDB.clone( connect = "frontier://FrontierProd/CMS_CONDITIONS" )
    process.PoolDBESSourceAlignment = cms.ESSource("PoolDBESSource", process.CondDBAlignment,
      toGet = cms.VPSet(cms.PSet(record = cms.string("RPRealAlignmentRecord"),tag = cms.string("CTPPSRPAlignment_real_offline_v8")))
    )
    process.esPreferDBFileAlignment = cms.ESPrefer("PoolDBESSource", "PoolDBESSourceAlignment")
    
    # override optics
    process.CondDBOptics = CondDB.clone( connect = "frontier://FrontierProd/CMS_CONDITIONS" )
    process.PoolDBESSourceOptics = cms.ESSource("PoolDBESSource",process.CondDBOptics,DumpStat = cms.untracked.bool(False),
      toGet = cms.VPSet(cms.PSet(record = cms.string('CTPPSOpticsRcd'),tag = cms.string("PPSOpticalFunctions_offline_v7")))
    )
    process.esPreferDBFileOptics = cms.ESPrefer("PoolDBESSource", "PoolDBESSourceOptics")
  else:
    process.load("CalibPPS.ESProducers.ctppsAlignment_cff")

  # reconstruction (load standard sequence)
  process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
  
  # reconstruction (remove modules which do not need re-running)
  process.recoCTPPSTask.remove(process.totemRPClusterProducer)
  process.recoCTPPSTask.remove(process.totemRPRecHitProducer)
  process.recoCTPPSTask.remove(process.ctppsPixelClusters)
  process.recoCTPPSTask.remove(process.ctppsPixelRecHits)

	  

def SetConditions(process):
  
  # chose LHCInfo
  global lhcInfoDefined
  lhcInfoDefined = True

  # chose alignment
  global alignmentDefined
  alignmentDefined = True

  # chose optics
  global opticsDefined
  opticsDefined = True
  
  # check choices
  if not lhcInfoDefined:
    raise ValueError("LHCInfo not defined")

  if not alignmentDefined:
    raise ValueError("alignment not defined")

  if not opticsDefined:
    raise ValueError("optics not defined")
	
  # processing path 
  

