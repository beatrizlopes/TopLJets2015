import FWCore.ParameterSet.Config as cms

def customizeJetTools(process,jecDB,jecTag,jerDB,jerTag,baseJetCollection,runOnData):

	#general configurations
	jecTag += '_DATA' if runOnData else '_MC'
        jerTag += '_DATA' if runOnData else '_MC'
	payload='AK4PFchs'
	jecLevels=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']
	print '[customizeJetTools]',jecDB,jecTag,payload,jecLevels

	#setup the source for JEC 
	process.load('CondCore.CondDB.CondDB_cfi')
	from CondCore.CondDB.CondDB_cfi import CondDB
	process.jec = cms.ESSource("PoolDBESSource",
				   DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
				   timetype = cms.string('runnumber'),
				   toGet = cms.VPSet( cms.PSet(record = cms.string('JetCorrectionsRecord'),
							       tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFchs'%jecTag),
							       label  = cms.untracked.string('AK4PFchs')
							       )						      
						      ), 
				   connect = cms.string('sqlite_file:%s'%jecDB)
				   )

	## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
	process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        
        process.jerDB = cms.ESSource('PoolDBESSource', CondDBSetup,
                                     connect = cms.string('sqlite_file:%s'%jerDB),
                                     toGet = cms.VPSet(cms.PSet(record = cms.string('JetResolutionRcd'),
                                                                tag = cms.string('JR_%s_PtResolution_AK4PFchs'%jerTag),
                                                                label = cms.untracked.string('AK4PFchs_pt')
                                                                ),
                                                       #cms.PSet(record = cms.string('JetResolutionRcd'),
                                                       #         tag = cms.string('JR_%s_PhiResolution_AK4PFchs'%jerTag),
                                                       #         label = cms.untracked.string('AK4PFchs_phi')
                                                       #         ),
                                                       cms.PSet(record = cms.string('JetResolutionScaleFactorRcd'),
                                                                tag = cms.string('JR_%s_SF_AK4PFchs'%jerTag),
                                                                label = cms.untracked.string('AK4PFchs')
                                                                ),
                                                       )
                                     )        
        process.jerDBPreference = cms.ESPrefer('PoolDBESSource', 'jerDB')

        process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
                                                CondDBSetup,
                                                toGet = cms.VPSet(cms.PSet(record = cms.string('QGLikelihoodRcd'),
                                                                           tag = cms.string('QGLikelihoodObject_v1_AK4'),
                                                                           label = cms.untracked.string('QGL_AK4PFchs')
                                                                           ),
                                                                  ),
                                                connect = cms.string('sqlite:qg_db.db')
                                                )
        process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')

        process.load('RecoJets.JetProducers.QGTagger_cfi')
        process.QGTagger.srcJets = baseJetCollection
        process.QGTagger.srcVertexCollection = 'offlineSlimmedPrimaryVertices'
        process.QGTagger.useQualityCuts = cms.bool(False)
        	
	from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
	updateJetCollection(
		process,
		labelName='UpdatedJECBTag',
		jetSource = cms.InputTag(baseJetCollection),
		jetCorrections = (payload, cms.vstring(jecLevels), 'None')
		)
    
	#add MET Muon filters
	#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#How_to_run_the_Bad_Charged_Hadro
	process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
	process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
	process.BadPFMuonFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices")
	process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

	from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
	process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
	  muons = cms.InputTag("slimmedMuons"),
	  vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
	  PFCandidates = cms.InputTag("packedPFCandidates"),
	  minDzBestTrack = cms.double(0.5),
	  taggingMode    = cms.bool(True)
	)

