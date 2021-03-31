import glob
print "year", "etabin", "sm_unc/sm_cor", "fe_unc/fe_cor", "sm_unc/fe_unc", "sm_cor/fe_cor"

for etabin in range(2,15):
  for year in ["UL16preVFP","UL16postVFP","UL17","UL18","Legacy"]:
    sm = open(glob.glob("file/eta_JER/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_ST.txt")[0]).readlines()
    fe = open(glob.glob("file/eta_JER/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_FE.txt")[0]).readlines()
    sm_eta = str(float(sm[etabin].split()[0])-float(sm[etabin].split()[1])) +" " +str(float(sm[etabin].split()[0])+float(sm[etabin].split()[1]))
    fe_eta = str(float(fe[etabin-1].split()[0])-float(fe[etabin-1].split()[1])) +" " +str(float(fe[etabin-1].split()[0])+float(fe[etabin-1].split()[1]))
    sm_unc = float(sm[etabin].split()[2])
    sm_cor = float(sm[etabin].split()[4])
    fe_unc = float(fe[etabin-1].split()[2])
    fe_cor = float(fe[etabin-1].split()[4])
    print etabin, year, (" "*(15-len(year))), sm_eta,fe_eta,(" "*(30-len(sm_eta)-len(fe_eta))), round(abs((sm_unc/sm_cor)-1)*100,1), "\t", round(abs((fe_unc/fe_cor)-1)*100,1), "\t", round(abs((sm_unc/fe_unc)-1)*100,1), "\t", round(abs((sm_cor/fe_cor)-1)*100,1)


print "year", "etabin", "sm_unc/sm_cor", "fe_unc/fe_cor", "sm_unc/fe_unc", "sm_cor/fe_cor"

for etabin in range(2,3):
  for year in ["UL16preVFP","UL16postVFP","UL17","UL18","Legacy"]:
    sm = open(glob.glob("file/eta_simple/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_ST.txt")[0]).readlines()
    fe = open(glob.glob("file/eta_simple/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_FE.txt")[0]).readlines()
    # print glob.glob("file/eta_simple/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_ST.txt")[0], glob.glob("file/eta_simple/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_FE.txt")[0]
    sm_eta = str(float(sm[etabin].split()[0])-float(sm[etabin].split()[1])) +" " +str(float(sm[etabin].split()[0])+float(sm[etabin].split()[1]))
    fe_eta = str(float(fe[etabin-1].split()[0])-float(fe[etabin-1].split()[1])) +" " +str(float(fe[etabin-1].split()[0])+float(fe[etabin-1].split()[1]))
    sm_unc = float(sm[etabin].split()[2])
    sm_cor = float(sm[etabin].split()[4])
    fe_unc = float(fe[etabin-1].split()[2])
    fe_cor = float(fe[etabin-1].split()[4])
    print etabin, year, (" "*(15-len(year))), sm_eta,fe_eta, (" "*(30-len(sm_eta)-len(fe_eta))), round(abs((sm_unc/sm_cor)-1)*100,1), "\t", round(abs((fe_unc/fe_cor)-1)*100,1), "\t", round(abs((sm_unc/fe_unc)-1)*100,1), "\t", round(abs((sm_cor/fe_cor)-1)*100,1)


print "year", "etabin", "simple_unc/simple_cor", "jer_unc/jer_cor", "simple_unc/jer_unc", "simple_cor/jer_cor"

for etabin in range(1,9):
  etabin_simple = 2 if etabin>3 else 1
  for year in ["UL16preVFP","UL16postVFP","UL17","UL18","Legacy"]:
    jer = open(glob.glob("file/eta_JER/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_FE.txt")[0]).readlines()
    simple = open(glob.glob("file/eta_simple/"+year+"/*/AK4CHS/standard/QCDHT/*/output/scalefactors_ST.txt")[0]).readlines()
    jer_eta = str(float(jer[etabin].split()[0])-float(jer[etabin].split()[1])) +" " +str(float(jer[etabin].split()[0])+float(jer[etabin].split()[1]))
    simple_eta = str(float(simple[etabin_simple].split()[0])-float(simple[etabin_simple].split()[1])) +" " +str(float(simple[etabin_simple].split()[0])+float(simple[etabin_simple].split()[1]))
    jer_unc = float(jer[etabin].split()[2])
    jer_cor = float(jer[etabin].split()[4])
    simple_unc = float(simple[etabin_simple].split()[2])
    simple_cor = float(simple[etabin_simple].split()[4])
    print etabin, year, (" "*(15-len(year))), jer_eta, simple_eta, (" "*(30-len(jer_eta)-len(simple_eta))), round(abs((simple_unc/simple_cor)-1)*100,1), "\t", round(abs((jer_unc/jer_cor)-1)*100,1), "\t", round(abs((simple_unc/jer_unc)-1)*100,1), "\t", round(abs((simple_cor/jer_cor)-1)*100,1)
