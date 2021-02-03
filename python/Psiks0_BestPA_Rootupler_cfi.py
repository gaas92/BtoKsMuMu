import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0_PVpa',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          GenParticles = cms.InputTag("genParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          isRes = cms.bool(True),     # if reco MC resonant  
                          OnlyGen = cms.bool(False),
                          )
