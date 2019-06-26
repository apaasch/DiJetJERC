#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TProfile.h>

#include <vector>

using namespace std;

void CorrectionObject::Derive_Thresholds_DiJet(bool pt_check, bool useHF){
  cout << "--------------- Starting Derive_Thresholds_DiJet() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  
  //cout<<"debug 1\n";
  
  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/thresholds/");

  int trg_nr = n_triggerDi;
  int n_trigger = n_triggerDi;
  vector<int> triggerVal(triggerValDi, triggerValDi + sizeof(triggerValDi) / sizeof(triggerValDi[0]));

  int n_trigger_HF = n_triggerDi_HF;
  vector<int> triggerVal_HF(triggerValDi_HF, triggerValDi_HF + sizeof(triggerValDi_HF) / sizeof(triggerValDi_HF[0]));  
  
  TH1D *hdata_pt_ave[trg_nr-1];
  TH1D *hdata_pt_ave_wNext[trg_nr-1];

  TH1D *hdata_pt_1[trg_nr-1];
  TH1D *hdata_pt_2[trg_nr-1];

  TH1D *hdata_pt_1_wNext[trg_nr-1];
  TH1D *hdata_pt_2_wNext[trg_nr-1];

   
  TH1D *hdata_pt_ave_HF[n_trigger_HF-1];
  TH1D *hdata_pt_ave_wNext_HF[n_trigger_HF-1];

  TH1D *hdata_pt_1_HF[n_trigger_HF-1];
  TH1D *hdata_pt_2_HF[n_trigger_HF-1];

  TH1D *hdata_pt_1_wNext_HF[n_trigger_HF-1];
  TH1D *hdata_pt_2_wNext_HF[n_trigger_HF-1];
  
    //cout<<"debug 2\n";
    
  for(int j=0; j<trg_nr-1; j++){
    TString name = "pt_ave_trg"+to_string(triggerVal[j]);
    TString name2 = "pt_ave_wNext_trg"+to_string(triggerVal[j]);

    hdata_pt_ave[j]= new TH1D(name,"",nResponseBins*6,0,// j<7?600:
			      1200);
    hdata_pt_ave_wNext[j]= new TH1D(name2,"",nResponseBins*6,0,// j<7?600:
				    1200);


    TString name_1 = "pt_1_trg"+to_string(triggerVal[j]);
    TString name2_1 = "pt_1_wNext_trg"+to_string(triggerVal[j]);
    TString name_2 = "pt_2_trg"+to_string(triggerVal[j]);
    TString name2_2 = "pt_2_wNext_trg"+to_string(triggerVal[j]);
      
    hdata_pt_1[j]= new TH1D(name_1,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_1_wNext[j]= new TH1D(name2_1,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_2[j]= new TH1D(name_2,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_2_wNext[j]= new TH1D(name2_2,"",nResponseBins*6,0,j<7?600:1200);
    
  }

  for(int j=0; j<n_trigger_HF-1; j++){
    TString name_HF = "pt_ave_trg"+to_string(triggerVal_HF[j]);
    TString name2_HF = "pt_ave_wNext_trg"+to_string(triggerVal_HF[j]);

    hdata_pt_ave_HF[j]= new TH1D(name_HF,"",nResponseBins*6,0,// j<7?600:
				 1200);
    hdata_pt_ave_wNext_HF[j]= new TH1D(name2_HF,"",nResponseBins*6,0,// j<7?600:
				       1200);


    TString name_1_HF = "pt_1_trg"+to_string(triggerVal_HF[j]);
    TString name2_1_HF = "pt_1_wNext_trg"+to_string(triggerVal_HF[j]);
    TString name_2_HF = "pt_2_trg"+to_string(triggerVal[j]);
    TString name2_2_HF = "pt_2_wNext_trg"+to_string(triggerVal_HF[j]);
      
    hdata_pt_1_HF[j]= new TH1D(name_1_HF,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_1_wNext_HF[j]= new TH1D(name2_1_HF,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_2_HF[j]= new TH1D(name_2_HF,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_2_wNext_HF[j]= new TH1D(name2_2_HF,"",nResponseBins*6,0,j<7?600:1200);
    
  }  
  
  //cout<<"debug 3\n";
      
  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<int> trg40(myReader_DATA, "triggerDi40");
  TTreeReaderValue<int> trg60(myReader_DATA, "triggerDi60");
  TTreeReaderValue<int> trg80(myReader_DATA, "triggerDi80");
  TTreeReaderValue<int> trg140(myReader_DATA, "triggerDi140");
  TTreeReaderValue<int> trg200(myReader_DATA, "triggerDi200");
  TTreeReaderValue<int> trg260(myReader_DATA, "triggerDi260");
  TTreeReaderValue<int> trg320(myReader_DATA, "triggerDi320");
  TTreeReaderValue<int> trg400(myReader_DATA, "triggerDi400");
  TTreeReaderValue<int> trg500(myReader_DATA, "triggerDi500");

  TTreeReaderValue<int> trg60_HF(myReader_DATA, "trigger60_HF");
  TTreeReaderValue<int> trg80_HF(myReader_DATA, "trigger80_HF");
  TTreeReaderValue<int> trg100_HF(myReader_DATA, "trigger100_HF");  
  TTreeReaderValue<int> trg160_HF(myReader_DATA, "trigger160_HF");
  TTreeReaderValue<int> trg220_HF(myReader_DATA, "trigger220_HF");
  TTreeReaderValue<int> trg300_HF(myReader_DATA, "trigger300_HF");
  
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");

  TTreeReaderValue<Float_t> pt_1_data(myReader_DATA, "jet1_pt");
  TTreeReaderValue<Float_t> pt_2_data(myReader_DATA, "jet2_pt");

    //cout<<"debug 4\n";
  
  TTreeReaderValue<int> trg_arr[trg_nr] = {trg40,trg60,trg80,trg140,trg200,trg260,trg320,trg400,trg500};
  TTreeReaderValue<int> trg_arr_HF[n_triggerDi_HF] = {trg60_HF,trg80_HF,trg100_HF,trg160_HF,trg220_HF,trg300_HF};
  
  //cout<<"debug 5\n";
      
  int myCount = 0;
  int myCount_notX = 0;
  bool allExclusive = true;
  while (myReader_DATA.Next()) {
    bool exclusive = true;
     //cout<<"debug 5.1\n";
    
    exclusive = (*trg40)^(*trg60)^(*trg80)^(*trg140)^(*trg200)^(*trg260)^(*trg320)^(*trg400)^(*trg500)^(useHF ? (*trg60_HF)^(*trg80_HF)^(*trg100_HF)^(*trg160_HF)^(*trg220_HF)^(*trg300_HF) : 0);
    for(int j=0; j<trg_nr-1; j++){
       //cout<<"debug 5.2\n";
      if(*(trg_arr[j])){
	hdata_pt_ave[j]->Fill(*pt_ave_data);
	if(pt_check){
	  hdata_pt_1[j]->Fill(*pt_1_data);
	  hdata_pt_2[j]->Fill(*pt_2_data);
	}
      }
      if(*(trg_arr[j+1])){
	hdata_pt_ave_wNext[j]->Fill(*pt_ave_data);
	   if(pt_check){
	     hdata_pt_1_wNext[j]->Fill(*pt_1_data);
	     hdata_pt_2_wNext[j]->Fill(*pt_2_data);
	   }
      }
	
    }
    if(useHF){
      for(int j=0; j<n_triggerDi_HF-1; j++){
	 //cout<<"debug 5.3\n";
	if(*(trg_arr_HF[j])){
	 //cout<<"debug 5.31\n";	  
	  hdata_pt_ave_HF[j]->Fill(*pt_ave_data);
	  //cout<<"debug 5.32\n";
	  if(pt_check){
	    hdata_pt_1_HF[j]->Fill(*pt_1_data);
	    hdata_pt_2_HF[j]->Fill(*pt_2_data);
	  }
	}
	if(*(trg_arr_HF[j+1])){
	  	 //cout<<"debug 5.33\n";	
	  hdata_pt_ave_wNext_HF[j]->Fill(*pt_ave_data);
	  	 //cout<<"debug 5.34\n";	
	  if(pt_check){
	    hdata_pt_1_wNext_HF[j]->Fill(*pt_1_data);
	    hdata_pt_2_wNext_HF[j]->Fill(*pt_2_data);
	  }
	}
	
      }
    }

     //cout<<"debug 5.4\n";
    myCount++;
    if(!exclusive){
      myCount_notX++;
    }
  }

  std::cout<<"\ncount data "<<myCount<<"  count data trg not exclusive "<<myCount_notX<<std::endl;
  
  TH1D* ptave_data_eff[n_trigger-1];

  TH1D* pt1_data_eff[n_trigger-1];
  TH1D* pt2_data_eff[n_trigger-1];

  TH1D* ptave_data_eff_HF[n_trigger_HF-1];

  TH1D* pt1_data_eff_HF[n_trigger_HF-1];
  TH1D* pt2_data_eff_HF[n_trigger_HF-1];

  
  for(int j=0; j<trg_nr-1; j++){
    hdata_pt_ave_wNext[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[j])+"_pt_ave_wNext"+".root");
    hdata_pt_ave[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[j])+"_pt_ave"+".root");

    hdata_pt_ave_wNext[j]->Rebin(6);
    hdata_pt_ave[j]->Rebin(6);
    ptave_data_eff[j]= (TH1D*) hdata_pt_ave_wNext[j]->Clone();    
    ptave_data_eff[j]->Divide((TH1D*) hdata_pt_ave[j]->Clone());
    ptave_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[j+1])+".root");

    if(pt_check){
      hdata_pt_1_wNext[j]->Rebin(6);
      hdata_pt_1[j]->Rebin(6);
      pt1_data_eff[j]= (TH1D*) hdata_pt_1_wNext[j]->Clone();
      pt1_data_eff[j]->Divide((TH1D*) hdata_pt_1[j]->Clone());
      pt1_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[j+1])+"_pt1.root");

      hdata_pt_2_wNext[j]->Rebin(6);
      hdata_pt_2[j]->Rebin(6);
      pt2_data_eff[j]= (TH1D*) hdata_pt_2_wNext[j]->Clone();
      pt2_data_eff[j]->Divide((TH1D*) hdata_pt_2[j]->Clone());
      pt2_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[j+1])+"_pt2.root");
    }
    
  }

  if(useHF){
    for(int j=0; j<n_triggerDi_HF-1; j++){
      hdata_pt_ave_wNext_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[j])+"_HF_pt_ave_wNext"+".root");
      hdata_pt_ave_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[j])+"_HF_pt_ave"+".root");

      hdata_pt_ave_wNext_HF[j]->Rebin(6);
      hdata_pt_ave_HF[j]->Rebin(6);
      ptave_data_eff_HF[j]= (TH1D*) hdata_pt_ave_wNext_HF[j]->Clone();    
      ptave_data_eff_HF[j]->Divide((TH1D*) hdata_pt_ave_HF[j]->Clone());
      ptave_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[j+1])+"_HF.root");

      if(pt_check){
	hdata_pt_1_wNext_HF[j]->Rebin(6);
	hdata_pt_1_HF[j]->Rebin(6);
	pt1_data_eff_HF[j]= (TH1D*) hdata_pt_1_wNext_HF[j]->Clone();
	pt1_data_eff_HF[j]->Divide((TH1D*) hdata_pt_1_HF[j]->Clone());
	pt1_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[j+1])+"_HF_pt1.root");

	hdata_pt_2_wNext_HF[j]->Rebin(6);
	hdata_pt_2_HF[j]->Rebin(6);
	pt2_data_eff_HF[j]= (TH1D*) hdata_pt_2_wNext_HF[j]->Clone();
	pt2_data_eff_HF[j]->Divide((TH1D*) hdata_pt_2_HF[j]->Clone());
	pt2_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[j+1])+"_HF_pt2.root");
      }
    
    }    
  }
  


  cout << "fit the thresholds" << endl << endl;

  bool use_for_extrapol[n_trigger-1];
  int n_extrapol=0;
  
  TF1 *func[n_trigger-1];
  double thresholds[n_trigger-1];
  double thresholds09[n_trigger-1];
  double thresholds_err[n_trigger-1];
  double thresholds09_err[n_trigger-1];
  double thresholds_errUp[n_trigger-1];
  double thresholds09_errUp[n_trigger-1];
  double thresholds_errDown[n_trigger-1];
  double thresholds09_errDown[n_trigger-1];
  double fitrange_up[n_trigger-1] = {0, -100, 300, 0, 0, 0, 0, 0};
  double fitrange_down[n_trigger-1] = {0, 0, -50, 150, 100, 80, 0, 400}; 
  for(int i=0; i<n_trigger-1; i++){
    TString fitname = "fit";
    fitname +=  to_string(triggerVal[i+1]); 
    func[i] = new TF1(fitname,SmoothFit,triggerVal[i]-20 + fitrange_down[i],
		      triggerVal[i+1]+200 +fitrange_up[i] ,3);
    func[i]->SetParameters(triggerVal[i+1], 40., 1.);
    func[i]->SetParNames("p0", "p1", "N");
    ptave_data_eff[i]->Fit(func[i],"R");

    // use_for_extrapol[i] = func[i]->GetParError(0)<100. && func[i]->GetParError(1)<100.;
    use_for_extrapol[i] = true;//triggerVal[i+1]>190 || triggerVal[i+1]<70;
    cout<< (func[i]->GetParError(0)<100. && func[i]->GetParError(1)<100.)<<endl;
    if(use_for_extrapol[i]) n_extrapol++;
    
    thresholds[i] = func[i]->GetX(0.95*func[i]->GetParameter(2), triggerVal[i+1]-10, 800);
    thresholds09[i] = func[i]->GetX(0.9*func[i]->GetParameter(2), triggerVal[i+1]-10, 800);

  }


  cout<<"Debug: use_for_extrapolation array ";
  for(int i=0; i<n_trigger-1; i++){
    cout<<use_for_extrapol[i];
  }
  cout<<endl;
  
  const double triggerVal_noLow[n_trigger-1] = {60, 80, 140, 200, 260, 320, 400, 500};
  cout << "extrapolate the threshold of the lower trigger from "<< n_extrapol<<" fitted thresholds" << endl << endl;
    
  double extrapol_x[n_extrapol];
  double extrapol_y[n_extrapol];
  
  // cout<<"Debug1"<<endl;
    
  int gCount = 0;
  for(int i=0; i<n_trigger-1; i++){
    if(!use_for_extrapol[i]) continue;
    extrapol_x[gCount]=triggerVal_noLow[i];
    extrapol_y[gCount]=thresholds[i];
    gCount++;
  }
  // cout<<"Debug2"<<endl;;
  double extrapol_x09[n_extrapol];
  double extrapol_y09[n_extrapol];
  gCount = 0;
  for(int i=0; i<n_trigger-1; i++){
    if(!use_for_extrapol[i]) continue;
    extrapol_x09[gCount]=triggerVal_noLow[i];
    extrapol_y09[gCount]=thresholds09[i];
    gCount++;
  }
  
  cout<<"Debug: extrapol_x array ";
  for(int i=0; i<n_extrapol; i++){
    cout<<extrapol_x[i]<<" ";
  }
  cout<<endl;
   cout<<"Debug: extrapol_y array "; 
  for(int i=0; i<n_extrapol; i++){
    cout<<extrapol_y[i]<<" ";
  }
  cout<<endl;

  
  // cout<<"Debug3"<<endl;;
  double all_thresholds[n_trigger];
  double all_thresholds09[n_trigger]; 
    
  TGraphErrors* gr095 = new TGraphErrors(n_extrapol,  extrapol_x , extrapol_y);
  TF1* func095 = new TF1("func095", "pol1" , 0, 501);
  func095->SetLineColor(kBlue);
  TFitResultPtr r = gr095->Fit(func095);
  all_thresholds[0]=gr095->Eval(40);
  for(int i=0; i<n_trigger-1; i++){
    if(use_for_extrapol[i]) all_thresholds[i+1]=thresholds[i];
    else all_thresholds[i+1]=gr095->Eval(triggerVal_noLow[i]);
  }  

  // cout<<"Debug4"<<endl;;
  TGraphErrors* gr09 = new TGraphErrors(n_extrapol, extrapol_x09 , extrapol_y09);
  TF1* func09 = new TF1("func09", "pol1" , 0, 501);
  r = gr09->Fit(func09);  
  all_thresholds09[0]=gr09->Eval(40);
  for(int i=0; i<n_trigger-1; i++){
    if(use_for_extrapol[i]) all_thresholds09[i+1]=thresholds09[i];
    else all_thresholds09[i+1]=gr09->Eval(triggerVal_noLow[i]);
  }  

  cout<<endl;
  for(int i=0; i<n_trigger; i++){
    cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal[i] << "; 0.95 threshold: "<< all_thresholds[i]  << endl;;
    cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal[i] << "; 0.90 threshold: "<< all_thresholds09[i]  << endl;
    cout << "Old Threshold: " << pt_bins[i] <<endl <<endl;  
  }
    
  cout << "Draw plots to "<<CorrectionObject::_outpath+"plots/thresholds" << endl;

  ofstream myfile;
  myfile.open (CorrectionObject::_outpath+"plots/thresholds/"+"/thresholds.txt",ios::trunc);
  myfile << "trigger  old  \t  0.90 \t  0.95  \n";
  for(int i=0; i<n_trigger; i++){
    myfile<< std::fixed << std::setprecision(2) << triggerVal[i]<< " \t " << pt_bins[i]<< " \t "<< all_thresholds09[i]<< " \t "<< all_thresholds[i] << " \n"; 
  }
  myfile.close();


  TCanvas* c = new TCanvas("c");
  // TMultiGraph* mg = new TMultiGraph();
  // mg->Add(gr09);
  // mg->Add(gr095);
  // mg->Draw();
  gr095->SetMarkerStyle(3);
  gr095->Draw("ap");
  c->Print(CorrectionObject::_outpath+"plots/thresholds/"+"extrapolateLowestTrigger095"+".pdf","pdf");    
   TCanvas* c2 = new TCanvas("c2");
   gr09->SetMarkerStyle(3);
   gr095->SetMarkerStyle(5);
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gr09);
  mg->Add(gr095);
  mg->Draw("ap");
  mg->GetXaxis()->SetLimits(0.,600.);
  mg->GetYaxis()->SetLimits(0.,550.);
  mg->GetXaxis()->SetTitle("p_{t}^{ave}");
  mg->GetYaxis()->SetTitle("threshold");    
  mg->Draw("ap");  
  c2->Print(CorrectionObject::_outpath+"plots/thresholds/"+"extrapolateLowestTrigger"+".pdf","pdf");

    for(int i=0; i<n_trigger-1; i++){
      TCanvas* c1 = new TCanvas("c1");
      ptave_data_eff[i]->Draw();
      ptave_data_eff[i]->GetXaxis()->SetTitle("p_{t}^{ave}");
      ptave_data_eff[i]->GetYaxis()->SetTitle("#epsilon");      
      TLine *line = new TLine(thresholds[i],0,thresholds[i],
			      1.05*func[i]->GetParameter(2));
      line->SetLineColor(kBlack);
      line->Draw();
       TLine *line09 = new TLine(thresholds09[i],0,thresholds09[i],
				 1.05*func[i]->GetParameter(2));
      line09->SetLineColor(kBlack);
      line09->SetLineStyle(2);
      line09->Draw();     
      c1->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[i+1])+".pdf","pdf");

      if(pt_check){
	TCanvas* c4 = new TCanvas("c4");
	pt1_data_eff[i]->Draw();
	pt1_data_eff[i]->GetXaxis()->SetTitle("p_{t}^{1}");
	pt1_data_eff[i]->GetYaxis()->SetTitle("#epsilon");      
	c4->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[i+1])+"_pt1.pdf","pdf");
	
	TCanvas* c3 = new TCanvas("c3");
	pt2_data_eff[i]->Draw();
	pt2_data_eff[i]->GetXaxis()->SetTitle("p_{t}^{2}");
	pt2_data_eff[i]->GetYaxis()->SetTitle("#epsilon");      
	c3->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal[i+1])+"_pt2.pdf","pdf");	
      }

  }

    //all in one plot
    TLegend *leg3 = new TLegend(0.6,0.15,0.8,0.45);
    leg3 -> SetBorderSize(0);
    leg3 -> SetTextSize(0.035);
    leg3 -> SetFillColor(0);
    
    TCanvas* c5 = new TCanvas("trig_central","",0,0,800,600);
    gStyle->SetOptStat(0);
    c5->SetFrameFillColor(0);

    for(int i=0; i<n_trigger-1; i++){
           TString fitname = "fit";
      fitname +=  to_string(triggerVal[i+1]);
      double N_func = ptave_data_eff[i]->GetFunction(fitname)->GetParameter(2);
      ptave_data_eff[i]->Scale(1./N_func);
      ptave_data_eff[i]->RebinX(5);
      ptave_data_eff[i]->Scale(1./5.);
      ptave_data_eff[i]->GetFunction(fitname)->SetBit(TF1::kNotDraw);
      
      if(i==0){
	ptave_data_eff[i]->GetYaxis()->SetTitle("#epsilon");
	ptave_data_eff[i]->GetYaxis()->SetTitleOffset(1.1);
	ptave_data_eff[i]->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
	ptave_data_eff[i]->GetXaxis()->SetTitleOffset(1.1);
      	ptave_data_eff[i]->SetMaximum(1.5);
	ptave_data_eff[i]->SetMinimum(0);
      }
      
      ptave_data_eff[i]->SetMarkerStyle(20);
      ptave_data_eff[i]->SetMarkerColor(kRed-i*5);
      ptave_data_eff[i]->Draw("P X0 SAME");
      
      TString fitname_n= fitname+"_norm";
      TF1* func_1 =  new TF1(fitname_n,SmoothFit,0,
		      1200 ,3);
      func_1->SetParameters(ptave_data_eff[i]->GetFunction(fitname)->GetParameter(0),ptave_data_eff[i]->GetFunction(fitname)->GetParameter(1), 1.);

      func_1->SetLineColor(kRed-i*5);
      func_1->Draw("SAME");

      TLine * line = new TLine(all_thresholds[i+1], 0,all_thresholds[i+1], 1.5);
      line->SetLineStyle(2);
      line->Draw("SAME");
      
      leg3 -> AddEntry(ptave_data_eff[i],TString("HLT_DiPFJetAve"+to_string(triggerValDi[i+1])),"p");
      
    }
	leg3->Draw("SAME");
	c5->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve_all.pdf","pdf");



	
    if(useHF){



  cout << "fit the thresholds of HF" << endl << endl;
 
  bool use_for_extrapol_HF[n_trigger_HF-1];
  int n_extrapol=0;
  
  TF1 *func_HF[n_trigger_HF-1];
  double thresholds_HF[n_trigger_HF-1];
  double thresholds09_HF[n_trigger_HF-1];
  double thresholds_err_HF[n_trigger_HF-1];
  double thresholds09_err_HF[n_trigger_HF-1];
  double thresholds_errUp_HF[n_trigger_HF-1];
  double thresholds09_errUp_HF[n_trigger_HF-1];
  double thresholds_errDown_HF[n_trigger_HF-1];
  double thresholds09_errDown_HF[n_trigger_HF-1];
  double fitrange_up_HF[n_trigger_HF-1] = {-100, 0,-150, 0, -140};
  double fitrange_down_HF[n_trigger_HF-1] = {-20, 0, 0, 0, 0}; 
  for(int i=0; i<n_trigger_HF-1; i++){
    TString fitname = "fit";
    fitname +=  to_string(triggerVal_HF[i+1]); 
    func_HF[i] = new TF1(fitname,SmoothFit,triggerVal_HF[i]-20 +fitrange_down_HF[i] ,
		      triggerVal_HF[i+1]+200 +fitrange_up_HF[i],3);
    func_HF[i]->SetParameters(triggerVal_HF[i+1], 40., 1.);
    func_HF[i]->SetParNames("p0", "p1", "N");
    ptave_data_eff_HF[i]->Fit(func_HF[i],"R");

    // use_for_extrapol_HF[i] = func_HF[i]->GetParError(0)<100. && func_HF[i]->GetParError(1)<100.;
    use_for_extrapol_HF[i] = true;//triggerVal_HF[i+1]>190 || triggerVal_HF[i+1]<70;
    cout<< (func_HF[i]->GetParError(0)<100. && func_HF[i]->GetParError(1)<100.)<<endl;
    if(use_for_extrapol_HF[i]) n_extrapol++;
    
    thresholds_HF[i] = func_HF[i]->GetX(0.95*func_HF[i]->GetParameter(2), triggerVal_HF[i+1]-10, 800);
    thresholds09_HF[i] = func_HF[i]->GetX(0.9*func_HF[i]->GetParameter(2), triggerVal_HF[i+1]-10, 800);

  }


  cout<<"Debug: use_for_extrapolation array ";
  for(int i=0; i<n_trigger_HF-1; i++){
    cout<<use_for_extrapol_HF[i];
  }
  cout<<endl;
  
  const double triggerVal_noLow_HF[n_trigger_HF-1] = {80, 100, 160, 220, 300};
  cout << "extrapolate the threshold of the lower trigger from "<< n_extrapol<<" fitted thresholds" << endl << endl;
    
  double extrapol_x_HF[n_extrapol];
  double extrapol_y_HF[n_extrapol];
  
  // cout<<"Debug1"<<endl;
    
  int gCount = 0;
  for(int i=0; i<n_trigger_HF-1; i++){
    if(!use_for_extrapol_HF[i]) continue;
    extrapol_x_HF[gCount]=triggerVal_noLow_HF[i];
    extrapol_y_HF[gCount]=thresholds_HF[i];
    gCount++;
  }
  // cout<<"Debug2"<<endl;;
  double extrapol_x09_HF[n_extrapol];
  double extrapol_y09_HF[n_extrapol];
  gCount = 0;
  for(int i=0; i<n_trigger_HF-1; i++){
    if(!use_for_extrapol_HF[i]) continue;
    extrapol_x09_HF[gCount]=triggerVal_noLow_HF[i];
    extrapol_y09_HF[gCount]=thresholds09_HF[i];
    gCount++;
  }
  
  cout<<"Debug: extrapol_x array ";
  for(int i=0; i<n_extrapol; i++){
    cout<<extrapol_x_HF[i]<<" ";
  }
  cout<<endl;
   cout<<"Debug: extrapol_y array "; 
  for(int i=0; i<n_extrapol; i++){
    cout<<extrapol_y_HF[i]<<" ";
  }
  cout<<endl;

  
  // cout<<"Debug3"<<endl;;
  double all_thresholds_HF[n_trigger_HF];
  double all_thresholds09_HF[n_trigger_HF]; 
    
  gr095 = new TGraphErrors(n_extrapol,  extrapol_x , extrapol_y);
  func095 = new TF1("func095", "pol1" , 0, 501);
  func095->SetLineColor(kBlue);
  TFitResultPtr r = gr095->Fit(func095);
  all_thresholds_HF[0]=gr095->Eval(40);
  for(int i=0; i<n_trigger_HF-1; i++){
    if(use_for_extrapol_HF[i]) all_thresholds_HF[i+1]=thresholds_HF[i];
    else all_thresholds_HF[i+1]=gr095->Eval(triggerVal_noLow_HF[i]);
  }  

  // cout<<"Debug4"<<endl;;
  gr09 = new TGraphErrors(n_extrapol, extrapol_x09 , extrapol_y09);
  TF1* func09 = new TF1("func09", "pol1" , 0, 501);
  r = gr09->Fit(func09);  
  all_thresholds09_HF[0]=gr09->Eval(40);
  for(int i=0; i<n_trigger_HF-1; i++){
    if(use_for_extrapol_HF[i]) all_thresholds09_HF[i+1]=thresholds09_HF[i];
    else all_thresholds09_HF[i+1]=gr09->Eval(triggerVal_noLow_HF[i]);
  }  

  cout<<endl;
  for(int i=0; i<n_trigger_HF; i++){
    cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal_HF[i] << "; 0.95 threshold: "<< all_thresholds_HF[i]  << endl;;
    cout<< std::fixed << std::setprecision(2) << "Trigger value: "<< triggerVal_HF[i] << "; 0.90 threshold: "<< all_thresholds09_HF[i]  << endl;
    cout << "Old Threshold: " << pt_bins_HF[i] <<endl <<endl;  
  }
    
  cout << "Draw plots to "<<CorrectionObject::_outpath+"plots/thresholds" << endl;

  ofstream myfile;
  myfile.open (CorrectionObject::_outpath+"plots/thresholds/"+"/thresholds_HF.txt",ios::trunc);
  myfile << "trigger  old  \t  0.90 \t  0.95  \n";
  for(int i=0; i<n_trigger_HF; i++){
    myfile<< std::fixed << std::setprecision(2) << triggerVal_HF[i]<< " \t " << pt_bins_HF[i]<< " \t "<< all_thresholds09_HF[i]<< " \t "<< all_thresholds_HF[i] << " \n"; 
  }
  myfile.close();


  TCanvas* cHF = new TCanvas("cHF");
  // TMultiGraph* mg = new TMultiGraph();
  // mg->Add(gr09);
  // mg->Add(gr095);
  // mg->Draw();
  gr095->SetMarkerStyle(3);
  gr095->Draw("ap");
  cHF->Print(CorrectionObject::_outpath+"plots/thresholds/"+"extrapolateLowestTrigger095"+"_HF.pdf","pdf");    
   TCanvas* cHF2 = new TCanvas("cHF2");
   gr09->SetMarkerStyle(3);
   gr095->SetMarkerStyle(5);
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gr09);
  mg->Add(gr095);
  mg->Draw("ap");
  mg->GetXaxis()->SetLimits(0.,600.);
  mg->GetYaxis()->SetLimits(0.,550.);
  mg->GetXaxis()->SetTitle("p_{t}^{ave}");
  mg->GetYaxis()->SetTitle("threshold");    
  mg->Draw("ap");  
  cHF2->Print(CorrectionObject::_outpath+"plots/thresholds/"+"extrapolateLowestTrigger"+"_HF.pdf","pdf");

    for(int i=0; i<n_trigger_HF-1; i++){
      TCanvas* cHF1 = new TCanvas("cHF1");
      ptave_data_eff_HF[i]->Draw();
      ptave_data_eff_HF[i]->GetXaxis()->SetTitle("p_{t}^{ave}");
      ptave_data_eff_HF[i]->GetYaxis()->SetTitle("#epsilon");      
      TLine *line = new TLine(thresholds_HF[i],0,thresholds_HF[i],
			      1.05*func_HF[i]->GetParameter(2));
      line->SetLineColor(kBlack);
      line->Draw();
       TLine *line09 = new TLine(thresholds09_HF[i],0,thresholds09_HF[i],
				 1.05*func_HF[i]->GetParameter(2));
      line09->SetLineColor(kBlack);
      line09->SetLineStyle(2);
      line09->Draw();     
      cHF1->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[i+1])+"_HF.pdf","pdf");

      if(pt_check){
	TCanvas* cHF4 = new TCanvas("cHF4");
	pt1_data_eff_HF[i]->Draw();
	pt1_data_eff_HF[i]->GetXaxis()->SetTitle("p_{t}^{1}");
	pt1_data_eff_HF[i]->GetYaxis()->SetTitle("#epsilon");      
	cHF4->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[i+1])+"_pt1_HF.pdf","pdf");
	
	TCanvas* cHF3 = new TCanvas("cHF3");
	pt2_data_eff_HF[i]->Draw();
	pt2_data_eff_HF[i]->GetXaxis()->SetTitle("p_{t}^{2}");
	pt2_data_eff_HF[i]->GetYaxis()->SetTitle("#epsilon");      
	cHF3->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_DiPFJetAve"+to_string(triggerVal_HF[i+1])+"_pt2_HF.pdf","pdf");	
      }

  }

    }
    

}
    
    


