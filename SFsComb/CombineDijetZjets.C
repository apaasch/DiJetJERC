// Combine Dijet and Z+jets JER SFs, stored as 2D histograms
// author: A.Karavdina
// date: 14.07.2019
// Run it with following command:
// root -l -b -q CombineDijetZjets.C\(\"Comb_ZjetV16_DijetV16\"\)

//exract x values for one y bin of TH2 histogram and store them as TGraphErrors (Z+jets)
TGraphErrors* extractGraph(TH2F *h1, double y_val){ 
  int n_y = h1->GetNbinsY();
  int n_x = h1->GetNbinsX();
  int y_bin =-1;
  for(int iy=1;iy<n_y;iy++){
    if(h1->GetYaxis()->GetBinCenter(iy)<=y_val && h1->GetYaxis()->GetBinCenter(iy+1)>y_val) 
      y_bin = iy;
  }
  double x_val[n_x], x_val_err[n_x], SF_val[n_x], SF_val_err[n_x];
  for(int ix=1;ix<n_x+1;ix++){
    x_val[ix-1] = h1->GetXaxis()->GetBinCenter(ix);
    x_val_err[ix-1] = 1e-2;
    SF_val[ix-1] = h1->GetBinContent(ix,y_bin); 
    SF_val_err[ix-1] = h1->GetBinError(ix,y_bin);
  }
  TGraphErrors* gr = new TGraphErrors(n_x,x_val,SF_val,x_val_err,SF_val_err);
  return gr;
}


//exract x values for one y bin of TH2Poly histogram and store them as TGraphErrors (Dijets)
TGraphErrors* extractGraph(TH2Poly *h1, double y_val){ 
  TList* binlist = h1->GetBins();
  vector<double> xval, xvalerr,SFval,SFvalerr; 
  for(const auto&& bin: *binlist){
    if(y_val>=((TH2PolyBin*)bin)->GetYMin() && y_val<((TH2PolyBin*)bin)->GetYMax()){
      SFval.push_back(h1->GetBinContent(((TH2PolyBin*)bin)->GetBinNumber()));
      SFvalerr.push_back(h1->GetBinError(((TH2PolyBin*)bin)->GetBinNumber()));
      xval.push_back(0.5*(((TH2PolyBin*)bin)->GetXMin()+((TH2PolyBin*)bin)->GetXMax()));
      xvalerr.push_back(1e-2);
    }
  }
  int n_x = SFval.size();
  double x_val[n_x], x_val_err[n_x], SF_val[n_x], SF_val_err[n_x];
  for(int i=0;i<n_x;i++){
    x_val[i] = xval[i];
    x_val_err[i] = xvalerr[i];
    SF_val[i] = SFval[i];
    SF_val_err[i] = SFvalerr[i];
  }

  TGraphErrors* gr = new TGraphErrors(n_x,x_val,SF_val,x_val_err,SF_val_err);
  return gr;
} 

//Clean points not filled due to low statistics
TGraphErrors* CleanEmptyPoints(TGraphErrors* input){
  double *Yval = input->GetY();
  double *YvalError = input->GetEY();
  double *Xval = input->GetX();
  double *XvalError = input->GetEX();
  int count=0;
  vector<double> Xnew,Ynew,Xerrornew,Yerrornew;
  for(int i=0;i<input->GetN();i++){
    //    if(YvalError[i]<1e-4 || Yval[i]==0) continue;
    //    if(Yval[i]==0) continue;
    if(Yval[i]==0 or Yval[i]<1) continue;//remove empty points and points with SFs<1
    count++;
    Xnew.push_back(Xval[i]);       
    Ynew.push_back(Yval[i]);
    Xerrornew.push_back(XvalError[i]);       
    Yerrornew.push_back(YvalError[i]);    
  }

  const int NnewSize =  count;
  double Xnew_m[NnewSize],Ynew_m[NnewSize],Xerrornew_m[NnewSize],Yerrornew_m[NnewSize]; //because silly ROOT doesn't know how to treat vectors
  for(int i=0;i<NnewSize;i++){
    Xnew_m[i] = Xnew[i];
    Ynew_m[i] = Ynew[i];
    Xerrornew_m[i] = Xerrornew[i];
    Yerrornew_m[i] = Yerrornew[i];
  }

  TGraphErrors* output = new TGraphErrors(count,Xnew_m,Ynew_m,Xerrornew_m,Yerrornew_m);
  if(input->GetN()!=output->GetN()) cout<<"Number of points in input: "<<input->GetN()<<" in output: "<<output->GetN()<<endl;
  return output;
}


//add Dijet systematics, hardcoded values are based on total JER SF uncertainty estimate V1 
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#2018_data
void addDijetSystemtics(TGraphErrors* gr, double eta_val){
  double systm_err = 0;
  // if(fabs(eta_val)<1.3)  systm_err = 0.10;
  // else if(fabs(eta_val)<2.85)  systm_err = 0.20;
  //   else if(fabs(eta_val)<2.96 && fabs(eta_val)>2.85) systm_err = 0.90;
  // else if(fabs(eta_val)<5.2) systm_err = 0.20;
  // //  cout<<"eta = "<<eta_val<<"  systm_err = "<<systm_err<<endl;
  if(fabs(eta_val)<1.3)  systm_err = 0.07;
  else if(fabs(eta_val)<2.5)  systm_err = 0.10;
    else if(fabs(eta_val)<2.96 && fabs(eta_val)>2.5) systm_err = 0.30;
  else if(fabs(eta_val)<5.2) systm_err = 0.15;
  for(int i;i<gr->GetN();i++)
    gr->SetPointError(i,1e-2,TMath::Hypot(gr->GetErrorY(i),systm_err));
}

//Step-like function (stolen from trigger thresholds fit)
Double_t SmoothFit(Double_t *v, Double_t *par){
  Double_t fitval  = 0.;
  if(par[2] != 0.){
    //    fitval = 0.5 * par[2] * (1. + TMath::Erf((v[0]-par[0]) / (TMath::Power(2, 0.5) * par[1] ) ) );
    fitval = 0.5 * par[2] * (par[3] + TMath::Erf((v[0]-par[0]) / (TMath::Power(2, 0.5) * par[1] ) ) );
  }
  return fitval;
}

Double_t SmoothFitECAL(Double_t *v, Double_t *par){
  Double_t fitval  = 0.;
  if(par[2] != 0.){
    fitval = sqrt(pow(par[0],2)/(v[0]*v[0])+pow(sqrt(par[3])*par[1],2)/v[0]+pow(par[3]*par[2],2)) / sqrt(pow(par[3]*par[0],2)/(v[0]*v[0])+pow(par[3]*par[1],2)/v[0]+pow(par[3]*par[2],2)) * par[3];
  }
  return fitval;
}

//fit 2 graphs, store result as plot
int doFit(TGraphErrors* gr1, TGraphErrors*gr2, TString lab1, TString lab2, TString leg_title, TString name, double chi2_thr, TString tag_name, bool skipZjets){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);
  double w = 600;
  double h = 600;
  gr1->SetMarkerColor(kRed);
  gr2->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(22);
  gr1->SetMarkerSize(1.2);
  gr2->SetMarkerSize(1.2);

  TMultiGraph *mg = new TMultiGraph();
  if(!skipZjets) mg->Add(gr1,"p");
  mg->Add(gr2,"p");
  TLegend* legend = new TLegend(0.15,0.75,0.38,0.95);
  legend->SetHeader(leg_title,"C"); // option "C" allows to center the header
  legend->AddEntry(gr1,lab1,"lep");
  legend->AddEntry(gr2,lab2,"lep");

  TLatex *tex2_rel = new TLatex();
  tex2_rel->SetNDC();
  tex2_rel->SetTextSize(0.035); 

  TCanvas* c0 = new TCanvas("SFCombined","Dijet and Z+jets JER SFs",w,h);
  mg->Draw("a");

  //Try with const fit first
  TF1* fit_func = new TF1("poly0","pol0");
  TFitResultPtr fit_SF = mg->Fit(fit_func,"FQN");
  float chi2 = fit_func->GetChisquare()/fit_func->GetNDF();
  TString chi2_val_str; 
  chi2_val_str.Form("%3.1f",chi2);
  TString chi2_rel = "#chi^{2}/n.d.f = ";
  chi2_rel +=chi2_val_str;
  TString par0_val_str; TString par0_err_str;
  float SF_const = fit_func->GetParameter(0);
  float SF_const_err = fit_func->GetParError(0);
  par0_val_str.Form("%3.2f",SF_const);
  par0_err_str.Form("%3.2f",SF_const_err);
  TString prt_SF_const = "Const Fit: ";
  prt_SF_const +=par0_val_str;  prt_SF_const += " +/- "; prt_SF_const += par0_err_str;
  tex2_rel->DrawLatex(0.54,0.85,prt_SF_const);
  tex2_rel->DrawLatex(0.54,0.80,chi2_rel);

  //If const fit does not look good, try step-like function
  if(chi2>chi2_thr){
    // //    TF1* fit_func_custom = new TF1("SF_custom",SmoothFit,40,600,3);
    // TF1* fit_func_custom = new TF1("SF_custom",SmoothFit,40,600,4);
    // //   fit_func_custom->SetParameters(200, 20., 2.5);
    // //    fit_func_custom->SetParameters(210, 10., 0.5, 4.0);
    // if(skipZjets) fit_func_custom->SetParameters(195, 30., 0.5, 4.0);
    // else fit_func_custom->SetParameters(180, 50., 1.8, 4.0);
    // fit_func_custom->SetParNames("p0", "p1", "N","p3");


    TF1* fit_func_custom = new TF1("SF_custom",SmoothFitECAL,40,600,4);
    //    fit_func_custom->SetParameters(5*2, 0, 0.05*sqrt(2), 2.020);
    fit_func_custom->SetParameters(10, 0.1, 0.05*sqrt(2), 1000);
    fit_func_custom->SetParNames("p0", "p1", "p2","p3");
    TFitResultPtr fit_SF2 = mg->Fit(fit_func_custom,"F");

    //    cout<<"chi2 fit: "<<fit_func_custom->GetChisquare()/fit_func_custom->GetNDF()<<endl;
    /*Create a TGraphErrors to hold the confidence intervals*/
   TGraphErrors *grint1 = new TGraphErrors(gr1->GetN());
   for (int i=0; i<gr1->GetN(); i++)
     grint1->SetPoint(i, gr1->GetX()[i], 0);
   /*Compute the confidence intervals at the x points of the created graph*/
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint1);
   for (int i=0; i<gr1->GetN(); i++)
     grint1->SetPoint(i, gr1->GetX()[i], fit_func_custom->Eval(gr1->GetX()[i]));
   TGraphErrors *grint2 = new TGraphErrors(gr2->GetN());
   for (int i=0; i<gr2->GetN(); i++)
     grint2->SetPoint(i, gr2->GetX()[i], 0);
   /*Compute the confidence intervals at the x points of the created graph*/
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint2);
   for (int i=0; i<gr2->GetN(); i++)
     grint2->SetPoint(i, gr2->GetX()[i], fit_func_custom->Eval(gr2->GetX()[i]));
   //Now the "grint" graph contains function values as its y-coordinates
   //and confidence intervals as the errors on these coordinates
   grint1->SetLineColor(kRed-1);
   grint1->SetFillColorAlpha(kRed-1,0.3);
   grint1->Draw("e3 same");
   grint2->SetLineColor(kBlue-1);
   grint2->SetFillColorAlpha(kBlue-1,0.3);
   grint2->Draw("e3 same");

   chi2 = fit_func_custom->GetChisquare()/fit_func_custom->GetNDF();
   chi2_val_str.Form("%3.1f",chi2);
   TString chi2_rel = "#chi^{2}/n.d.f = ";
   chi2_rel +=chi2_val_str;
   TString step_par0_val_str; TString step_par0_err_str;
   TString step_par1_val_str; TString step_par1_err_str;
   TString step_par2_val_str; TString step_par2_err_str;
   TString step_par3_val_str; TString step_par3_err_str;
   float Step_par0 = fit_func_custom->GetParameter(0);
   float Step_par0_err = fit_func_custom->GetParError(0);
   float Step_par1 = fit_func_custom->GetParameter(1);
   float Step_par1_err = fit_func_custom->GetParError(1);
   float Step_par2 = fit_func_custom->GetParameter(2);
   float Step_par2_err = fit_func_custom->GetParError(2);
   float Step_par3 = fit_func_custom->GetParameter(3);
   float Step_par3_err = fit_func_custom->GetParError(3);

   step_par0_val_str.Form("%3.2f",Step_par0);
   step_par0_err_str.Form("%3.2f",Step_par0_err);
   step_par1_val_str.Form("%3.2f",Step_par1);
   step_par1_err_str.Form("%3.2f",Step_par1_err);
   step_par2_val_str.Form("%3.2f",Step_par2);
   step_par2_err_str.Form("%3.2f",Step_par2_err);
   step_par3_val_str.Form("%3.2f",Step_par3);
   step_par3_err_str.Form("%3.2f",Step_par3_err);

   TString prt0_SF_step = "par0 = ";
   prt0_SF_step +=step_par0_val_str;  prt0_SF_step += " +/- "; prt0_SF_step += step_par0_err_str;
   TString prt1_SF_step = "par1 = ";
   prt1_SF_step +=step_par1_val_str;  prt1_SF_step += " +/- "; prt1_SF_step += step_par1_err_str;
   TString prt2_SF_step = "par2 = ";
   prt2_SF_step +=step_par2_val_str;  prt2_SF_step += " +/- "; prt2_SF_step += step_par2_err_str;
   TString prt3_SF_step = "par3 = ";
   prt3_SF_step +=step_par3_val_str;  prt3_SF_step += " +/- "; prt3_SF_step += step_par3_err_str;

   tex2_rel->DrawLatex(0.54,0.33,"Step Fit:");
   tex2_rel->DrawLatex(0.54,0.29,prt0_SF_step);
   tex2_rel->DrawLatex(0.54,0.25,prt1_SF_step);
   tex2_rel->DrawLatex(0.54,0.21,prt2_SF_step);
   tex2_rel->DrawLatex(0.54,0.17,prt3_SF_step);
   tex2_rel->DrawLatex(0.54,0.13,chi2_rel);
  }
  else{
    TFitResultPtr fit_SF = mg->Fit(fit_func,"F");

    /*Create a TGraphErrors to hold the confidence intervals*/
   TGraphErrors *grint1 = new TGraphErrors(gr1->GetN());
   for (int i=0; i<gr1->GetN(); i++)
     grint1->SetPoint(i, gr1->GetX()[i], 0);
   /*Compute the confidence intervals at the x points of the created graph*/
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint1);
   for (int i=0; i<gr1->GetN(); i++)
     grint1->SetPoint(i, gr1->GetX()[i], fit_func->Eval(gr1->GetX()[i]));
   TGraphErrors *grint2 = new TGraphErrors(gr2->GetN());
   for (int i=0; i<gr2->GetN(); i++)
     grint2->SetPoint(i, gr2->GetX()[i], 0);
   /*Compute the confidence intervals at the x points of the created graph*/
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint2);
   for (int i=0; i<gr2->GetN(); i++)
     grint2->SetPoint(i, gr2->GetX()[i], fit_func->Eval(gr2->GetX()[i]));
   //Now the "grint" graph contains function values as its y-coordinates
   //and confidence intervals as the errors on these coordinates
   grint1->SetLineColor(kRed-1);
   grint1->SetFillColorAlpha(kRed-1,0.3);
   grint1->Draw("e3 same");
   grint2->SetLineColor(kBlue-1);
   grint2->SetFillColorAlpha(kBlue-1,0.3);
   grint2->Draw("e3 same");
  }

  mg->GetXaxis()->SetTitle("p_{T}");
  mg->GetYaxis()->SetTitle("JER SF");
  mg->GetYaxis()->SetRangeUser(0,3.0);
  legend->Draw();
  c0->SaveAs(name+"_"+tag_name+".pdf");
  return 0;

}


void CombineDijetZjets(TString tag_name="Ngenjet"){
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045,"x");  
  gStyle->SetTitleSize(0.045,"y");
  gStyle->SetTitleYOffset(0.9);

 double w = 600;
 double h = 600;


 //Files after dijet and Z+jets standalone analysis
 // TString path_dijet = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/MergeL2Res/Autumn18_V16h/AK4CHS/standard/QCDHT/RunABCD/output/"; //V16h
 TString path_dijet = "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v1/CMSSW_10_2_10/src/UHH2/DiJetJERC/JERSF_Analysis/JER/wide_eta_binning/file/MergeL2Res/Autumn18_V16/AK4CHS/standard/QCDHT/RunABCD/output/"; //V16
 // TString path_dijet = "/nfs/dust/cms/user/karavdia/CMSSW_10_2_10/src/UHH2/JERSF/Analysis/JER/wide_eta_bin_moreECbins_morePtbins_EnergyEtaCut_fixJetSorting/file/StandardPtBins/Autumn18_V15/AK4CHS/standard/QCDHT/RunABCD/output/"; //V15
 // TString path_dijet = "/nfs/dust/cms/user/karavdia/CMSSW_10_2_10/src/UHH2/JERSF/Analysis/JER/wide_eta_bin_moreECbins_morePtbins_EnergyEtaCut_fixJetSorting_positiveEta/file/StandardPtBins/Autumn18_V15/AK4CHS/standard/QCDHT/RunABCD/output/";
 TString name_dijet = "DijetJERSF2D.root";

 // TString path_zjet = "/afs/desy.de/user/k/karavdia/www/JEC_plots/ZJets_JER_combination/RunABCD/";
 // TString name_zjet = "h2_jer_sf.root";

 // TString path_zjet = "/afs/desy.de/user/k/karavdia/www/JEC_plots/ZJets_JER_combination/August19/";
 TString path_zjet = "/afs/desy.de/user/k/karavdia/www/JEC_plots/ZJets_JER_combination/August19/abs_eta/";
 TString name_zjet = "results.root";

 TString gl_label_dijet = "Dijet";
 TString gl_label_zjet = "Zjets";

 TFile *input_dijet = TFile::Open(path_dijet+name_dijet);
 TFile *input_zjet = TFile::Open(path_zjet+name_zjet);

 //input Z+jets root files contain TCanvas with hists
 // TCanvas *c1 = (TCanvas*)input_zjet->Get("c4");
 // TH2F *h_zjets = (TH2F*)c1->GetPrimitive("jersf");
 TH2F *h_zjets = (TH2F*)input_zjet->Get("JER_SF");

 //input dijet root files contain hists stored directly
 TH2Poly *h_dijets = (TH2Poly*)input_dijet->Get("2D_SF_FE");


 // eta bins:
 int n_eta_bins = 15;
 double eta_bins[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 5.191};
 TString eta_bins_str[] = {"0", "0.522", "0.783", "1.131", "1.305", "1.740", "1.930", "2.043", "2.322", "2.5", "2.65", "2.853", "2.964", "3.139", "5.191"};
 TString eta_bins_name [] = {"0", "0_522", "0_783", "1_131", "1_305", "1_740", "1_930", "2_043", "2_322", "2_5", "2_65", "2_853", "2_964", "3_139", "5_191"};
 TGraphErrors* gr_zjets; TGraphErrors* gr_dijets;

 for(int ieta=0;ieta< n_eta_bins-1;ieta++){
   double eta_val =  eta_bins[ieta]+1e-2;
   cout<<"eta_val = "<<eta_val<<endl;
   bool skipZjets = false;
   // if(eta_val>2.65 && eta_val<2.85) skipZjets = true;//FixME: adjust Z+jets binning
   gr_zjets = extractGraph(h_zjets, eta_val);
   addDijetSystemtics(gr_zjets,eta_val);//Add total systematics as for dijet
   gr_dijets = extractGraph(h_dijets, eta_val);
   gr_dijets = CleanEmptyPoints(gr_dijets);
   addDijetSystemtics(gr_dijets,eta_val);//Add total systematics for dijet
   TString legend_title = eta_bins_str[ieta];
   legend_title+="<|#eta|<"; legend_title+=eta_bins_str[ieta+1];
   //   int fitres = doFit(gr_zjets, gr_dijets, "Z+jets", "Dijets",legend_title,"JERSF_Eta_"+eta_bins_name[ieta],2, tag_name, skipZjets);
   int fitres = doFit(gr_zjets, gr_dijets, "Z+jets", "Dijets",legend_title,"JERSF_Eta_"+eta_bins_name[ieta],3.0, tag_name, skipZjets);
 }
 

}
