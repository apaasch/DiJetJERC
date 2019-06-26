LIBRARY := SUHH2DiJetJERC
#DICT := include/Dijet_LinkDef.h
#DICT := include/ReconstructionHypothesis.h include/SUHH2common_LinkDef.h
#TEST := 1
USERLDFLAGS := -lSUHH2core -lSUHH2common
USERCXXFLAGS := -Wno-reorder -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-parameter
PAR := 1
include ../Makefile.common


