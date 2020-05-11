import FWCore.ParameterSet.Config as cms


def setupProtonSim(process,xangle,withPU=False):

  # update settings of beam-smearing module
  process.beamDivergenceVtxGenerator.src = cms.InputTag("")
  
  # input collections
  if withPU:
    process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(cms.InputTag("genPUProtons","genPUProtons"),
    cms.InputTag("prunedGenParticles"))
  else:
    process.beamDivergenceVtxGenerator.srcGenParticle = cms.VInputTag(cms.InputTag("prunedGenParticles"))
  
  # do not apply vertex smearing again
  process.ctppsBeamParametersESSource.vtxStddevX = 0
  process.ctppsBeamParametersESSource.vtxStddevY = 0
  process.ctppsBeamParametersESSource.vtxStddevZ = 0
  
  #undo CMS vertex shift
  process.ctppsBeamParametersESSource.vtxOffsetX45 = +0.2475 * 1E-1
  process.ctppsBeamParametersESSource.vtxOffsetY45 = -0.6924 * 1E-1
  process.ctppsBeamParametersESSource.vtxOffsetZ45 = -8.1100 * 1E-1
  
  # set up Xangle  
  process.ctppsLHCInfoESSource.xangle = cms.double(xangle)
  process.ctppsBeamParametersESSource.halfXangleX45 = xangle * 1E-6
  process.ctppsBeamParametersESSource.halfXangleX56 = xangle * 1E-6

  # processing path 
  process.pps_fastsim = cms.Path(process.beamDivergenceVtxGenerator
    * process.ctppsDirectProtonSimulation
    * process.reco_local
    * process.ctppsProtons
 )

def ctppsCustom(process):

  # define conditions
  SetConditions(process)
  
	

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
  

