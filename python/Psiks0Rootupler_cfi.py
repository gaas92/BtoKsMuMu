import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          GenParticles = cms.InputTag("genParticles"),
                          #packedGenParticles = cms.InputTag("packedGenParticles"),
                          packedGenParticles = cms.InputTag("genParticlesPlusSim"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(False),
                          isRes = cms.bool(False),     # if reco MC resonant  
                          OnlyGen = cms.bool(False),
                          )
