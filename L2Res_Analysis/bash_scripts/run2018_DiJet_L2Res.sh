path_Gl_MC=/nfs/dust/cms/user/karavdia/DiJetJERC/RunII_102X_v1/2018/Autumn18V13h_Closure/
path_Gl_DATA=/nfs/dust/cms/user/karavdia/DiJetJERC/RunII_102X_v1/2018/Autumn18V13h_Closure/

cd ../src/

./main --run ABC -FP --inputMC $path_Gl_MC/uhh2.AnalysisModuleRunner.MC.QCDHT.root --input $path_Gl_DATA/uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

./main --run ABC -aAP --inputMC $path_Gl_MC/uhh2.AnalysisModuleRunner.MC.QCDHT.root --input $path_Gl_DATA/uhh2.AnalysisModuleRunner.DATA.DATA_RunABC_2018.root --outSuffix _NoJERsf_V13hL2L3Res_MadGraph --Generator madgraph

cd ../bash_scripts/