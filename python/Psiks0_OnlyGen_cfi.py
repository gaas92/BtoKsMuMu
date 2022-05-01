import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0_onlyGen',
                          dimuons = cms.InputTag("slimmedMuons"),
                          GenParticles = cms.InputTag("genParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          Lost_Tracks = cms.InputTag("lostTracks"),
                          pPFC_Tracks = cms.InputTag("packedPFCandidates"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          isRes = cms.bool(False),     # if reco MC resonant  
                          OnlyGen = cms.bool(True),
                          )
 