import FWCore.ParameterSet.Config as cms


def setupProtonSim(process,xangle,withPU=False):

  # load settings
  #process.load("direct_simu_reco_cff")
   
  # update settings of beam-smearing module
  process.load("IOMC.EventVertexGenerators.beamDivergenceVtxGenerator_cfi")
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

  # processing path 
  process.pps_fastsim = cms.Path(process.beamDivergenceVtxGenerator
    * process.ctppsDirectProtonSimulation
    * process.reco_local
    * process.ctppsProtons
 )

