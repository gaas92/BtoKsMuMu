# MiniAODBphysicsUltraLegacyRun2

## Intructions to run the examples.
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_10_6_12
cd CMSSW_10_6_12/src/
cmsenv
voms-proxy-init -voms cms -valid 192:00
git clone https://github.com/jmejiagu/MiniAODBphysicsUltraLegacyRun2.git myAnalyzers/BtoKsMuMu
cd myAnalyzers/BtoKsMuMu/
git checkout master
cd ../..
scram b

```

Run: (use your favorite input sample. You will see examples in the confi files)


```
cd myAnalyzers/BtoKsMuMu/test/
cmsRun PsikaonRootupler.py
```

This example is for Bu hadron. In test directory you can find other examples to run Bs and Bd hadrons too.