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

void CorrectionObject::Derive_Thresholds_SiJet(bool pt_check, bool useHF){
  cout << "--------------- Starting Derive_Thresholds_SiJet() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/thresholds/");

  int trg_nr=n_triggerSi;
  int n_trigger = n_triggerSi;
  vector<int> triggerVal(triggerValSi, triggerValSi + sizeof(triggerValSi) / sizeof(triggerValSi[0]));

  int n_trigger_HF = n_triggerSi_HF; 
  vector<int> triggerVal_HF(triggerValSi_HF, triggerValSi_HF + sizeof(triggerValSi_HF) / sizeof(triggerValSi_HF[0]));  
  
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
  
  
  for(int j=0; j<trg_nr-1; j++){
    TString name = "pt_ave_trg"+to_string(triggerVal[j]);
    TString name2 = "pt_ave_wNext_trg"+to_string(triggerVal[j]);

    hdata_pt_ave[j]= new TH1D(name,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_ave_wNext[j]= new TH1D(name2,"",nResponseBins*6,0,j<7?600:1200);


    TString name_1 = "pt_1_trg"+to_string(triggerVal[j]);
    TString name2_1 = "pt_1_wNext_trg"+to_string(triggerVal[j]);
    TString name_2 = "pt_2_trg"+to_string(triggerVal[j]);
    TString name2_2 = "pt_2_wNext_trg"+to_string(triggerVal[j]);
      
    hdata_pt_1[j]= new TH1D(name_1,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_1_wNext[j]= new TH1D(name2_1,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_2[j]= new TH1D(name_2,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_2_wNext[j]= new TH1D(name2_2,"",nResponseBins*3,0,j<7?600:1200);

    TString name_HF = "pt_ave_trg"+to_string(triggerVal_HF[j])+string("_HF");
    TString name2_HF = "pt_ave_wNext_trg"+to_string(triggerVal_HF[j])+string("_HF");

    hdata_pt_ave_HF[j]= new TH1D(name_HF,"",nResponseBins*6,0,j<7?600:1200);
    hdata_pt_ave_wNext_HF[j]= new TH1D(name2_HF,"",nResponseBins*6,0,j<7?600:1200);


    TString name_1_HF = "pt_1_trg"+to_string(triggerVal_HF[j])+"_HF";
    TString name2_1_HF = "pt_1_wNext_trg"+to_string(triggerVal_HF[j])+"_HF";
    TString name_2_HF = "pt_2_trg"+to_string(triggerVal_HF[j])+"_HF";
    TString name2_2_HF = "pt_2_wNext_trg"+to_string(triggerVal_HF[j])+"_HF";
      
    hdata_pt_1_HF[j]= new TH1D(name_1_HF,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_1_wNext_HF[j]= new TH1D(name2_1_HF,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_2_HF[j]= new TH1D(name_2_HF,"",nResponseBins*3,0,j<7?600:1200);
    hdata_pt_2_wNext_HF[j]= new TH1D(name2_2_HF,"",nResponseBins*3,0,j<7?600:1200);
    
  }
  
  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<int> trg40(myReader_DATA, "trigger40");
  TTreeReaderValue<int> trg60(myReader_DATA, "trigger60");
  TTreeReaderValue<int> trg80(myReader_DATA, "trigger80");
  TTreeReaderValue<int> trg140(myReader_DATA, "trigger140");
  TTreeReaderValue<int> trg200(myReader_DATA, "trigger200");
  TTreeReaderValue<int> trg260(myReader_DATA, "trigger260");
  TTreeReaderValue<int> trg320(myReader_DATA, "trigger320");
  TTreeReaderValue<int> trg400(myReader_DATA, "trigger400");
  TTreeReaderValue<int> trg450(myReader_DATA, "trigger450");
  TTreeReaderValue<int> trg500(myReader_DATA, "trigger500");

  cout<<"after central readers\n";

  TTreeReaderValue<int> trg60_HF(myReader_DATA, "trigger60_HF");
  TTreeReaderValue<int> trg80_HF(myReader_DATA, "trigger80_HF");
  TTreeReaderValue<int> trg140_HF(myReader_DATA, "trigger140_HF");  
  TTreeReaderValue<int> trg200_HF(myReader_DATA, "trigger200_HF");
  TTreeReaderValue<int> trg260_HF(myReader_DATA, "trigger260_HF");
  TTreeReaderValue<int> trg320_HF(myReader_DATA, "trigger320_HF");
  TTreeReaderValue<int> trg400_HF(myReader_DATA, "trigger400_HF");
  TTreeReaderValue<int> trg450_HF(myReader_DATA, "trigger450_HF");
  TTreeReaderValue<int> trg500_HF(myReader_DATA, "trigger500_HF");

  cout<<"after HF readers\n";  
  
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");

  TTreeReaderValue<Float_t> pt_1_data(myReader_DATA, "jet1_pt");
  TTreeReaderValue<Float_t> pt_2_data(myReader_DATA, "jet2_pt");
  

  TTreeReaderValue<int> trg_arr[trg_nr] = {trg40,trg60,trg80,trg140,trg200,trg260,trg320,trg400,trg450,trg500};

   TTreeReaderValue<int> trg_arr_HF[n_trigger_HF] = {trg60_HF,trg80_HF,trg140_HF,trg200_HF,trg260_HF,trg320_HF,trg400_HF,trg450_HF,trg500_HF}; 
  
  int myCount = 0;
  int myCount_notX = 0;
  bool allExclusive = true;
  while (myReader_DATA.Next()) {
    bool exclusive = true;
    exclusive = (*trg40)^(*trg60)^(*trg80)^(*trg140)^(*trg200)^(*trg260)^(*trg320)^(*trg400)^(*trg450)^(*trg500)^(useHF ? (*trg60_HF)^(*trg80_HF)^(*trg140_HF)^(*trg200_HF)^(*trg260_HF)^(*trg320_HF)^(*trg400_HF)^(*trg450_HF)^(*trg500_HF) : 0);
    for(int j=0; j<trg_nr-1; j++){
      if(*(trg_arr[j])){
	hdata_pt_ave[j]->Fill(*pt_ave_data);\
	// cout<<"debug after pt ave fill for central trg\n";
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
      // bool hftest = false;
      for(int j=0; j<n_trigger_HF-1; j++){
	 // cout<<"debug 5.3\n";
	// hftest = hftest||*(trg_arr_HF[j]);
       	if(*(trg_arr_HF[j])){
	 // cout<<"debug 5.31\n";
	 // cout<<(*(pt_ave_data))<<endl;
	 // cout<<"debug 5.311\n";	
	 
	  hdata_pt_ave_HF[j]->Fill(*pt_ave_data);
	  // cout<<"debug 5.32\n";
	  if(pt_check){
	    hdata_pt_1_HF[j]->Fill(*pt_1_data);
	    hdata_pt_2_HF[j]->Fill(*pt_2_data);
	  }
	}
	if(*(trg_arr_HF[j+1])){
	  	 // cout<<"debug 5.33\n";
	 // cout<<(*(pt_ave_data))<<endl;
	 // cout<<"debug 5.311\n";			 
	  hdata_pt_ave_wNext_HF[j]->Fill(*pt_ave_data);
	  	 // cout<<"debug 5.34\n";	
	  if(pt_check){
	    hdata_pt_1_wNext_HF[j]->Fill(*pt_1_data);
	    hdata_pt_2_wNext_HF[j]->Fill(*pt_2_data);
	  }
	}
	
      }
      // cout<<"is there hf stuff found in the readers?  "<<hftest<<endl;
    }    
    
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
    hdata_pt_ave_wNext[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[j])+"_pt_ave_wNext"+".root");
    hdata_pt_ave[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[j])+"_pt_ave"+".root");

    hdata_pt_ave_wNext[j]->Rebin(6);
    hdata_pt_ave[j]->Rebin(6);
    ptave_data_eff[j]= (TH1D*) hdata_pt_ave_wNext[j]->Clone();    
    ptave_data_eff[j]->Divide((TH1D*) hdata_pt_ave[j]->Clone());
    ptave_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[j+1])+".root");

    if(pt_check){
      hdata_pt_1_wNext[j]->Rebin(6);
      hdata_pt_1[j]->Rebin(6);
      pt1_data_eff[j]= (TH1D*) hdata_pt_1_wNext[j]->Clone();
      pt1_data_eff[j]->Divide((TH1D*) hdata_pt_1[j]->Clone());
      pt1_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[j+1])+"_pt1.root");

      hdata_pt_2_wNext[j]->Rebin(6);
      hdata_pt_2[j]->Rebin(6);
      pt2_data_eff[j]= (TH1D*) hdata_pt_2_wNext[j]->Clone();
      pt2_data_eff[j]->Divide((TH1D*) hdata_pt_2[j]->Clone());
      pt2_data_eff[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[j+1])+"_pt2.root");
    }
    
  }

   if(useHF){
     std::cout<<"\nstarted HF saving loop\n";
    for(int j=0; j<n_trigger_HF-1; j++){
      std::cout<<j<<std::endl;
      hdata_pt_ave_wNext_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[j])+"_HF_pt_ave_wNext"+".root");
      if(j==0)  std::cout<<"after first HF save"<<std::endl;
      hdata_pt_ave_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[j])+"_HF_pt_ave"+".root");

      hdata_pt_ave_wNext_HF[j]->Rebin(6);
      hdata_pt_ave_HF[j]->Rebin(6);
      ptave_data_eff_HF[j]= (TH1D*) hdata_pt_ave_wNext_HF[j]->Clone();    
      ptave_data_eff_HF[j]->Divide((TH1D*) hdata_pt_ave_HF[j]->Clone());
      ptave_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[j+1])+"_HF.root");

      if(pt_check){
	hdata_pt_1_wNext_HF[j]->Rebin(6);
	hdata_pt_1_HF[j]->Rebin(6);
	pt1_data_eff_HF[j]= (TH1D*) hdata_pt_1_wNext_HF[j]->Clone();
	pt1_data_eff_HF[j]->Divide((TH1D*) hdata_pt_1_HF[j]->Clone());
	pt1_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[j+1])+"_HF_pt1.root");

	hdata_pt_2_wNext_HF[j]->Rebin(6);
	hdata_pt_2_HF[j]->Rebin(6);
	pt2_data_eff_HF[j]= (TH1D*) hdata_pt_2_wNext_HF[j]->Clone();
	pt2_data_eff_HF[j]->Divide((TH1D*) hdata_pt_2_HF[j]->Clone());
	pt2_data_eff_HF[j]->SaveAs(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[j+1])+"_HF_pt2.root");
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
  double fitrange_up[n_trigger-1] = {0, -100, 0, 100, 100, 50, 0, 400, 0};
  double fitrange_down[n_trigger-1] = {0, 20, 0, 0, 50, 50, 0, -200, 0}; 
  for(int i=0; i<n_trigger-1; i++){
    TString fitname = "fit";
    fitname +=  to_string(triggerVal[i+1]); 
    func[i] = new TF1(fitname,SmoothFit,triggerVal[i]-(triggerVal[i+1]<490 ? 20 : 150) + fitrange_down[i],
		      triggerVal[i+1]+(triggerVal[i+1]<490 ? 200 : 500) + fitrange_up[i],3);
    func[i]->SetParameters(triggerVal[i+1], 40., 1.);
    func[i]->SetParNames("p0", "p1", "N");
    ptave_data_eff[i]->Fit(func[i],"R");

    use_for_extrapol[i] = ptave_data_eff[i]->GetEntries()>100 ;
    use_for_extrapol[i] = true;
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
  
  const double triggerVal_noLow[n_trigger-1] = {60, 80, 140, 200, 260, 320, 400, 450, 500};
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
  //  all_thresholds[0]=gr095->Eval(40);
  all_thresholds[0]=func095->Eval(40);
  for(int i=0; i<n_trigger-1; i++){
    if(use_for_extrapol[i]) all_thresholds[i+1]=thresholds[i];
    else all_thresholds[i+1]=gr095->Eval(triggerVal_noLow[i]);
  }  

  // cout<<"Debug4"<<endl;;
  TGraphErrors* gr09 = new TGraphErrors(n_extrapol, extrapol_x09 , extrapol_y09);
  TF1* func09 = new TF1("func09", "pol1" , 0, 501);
  r = gr09->Fit(func09);  
  //  all_thresholds09[0]=gr09->Eval(40);
  all_thresholds09[0]=func09->Eval(40);
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
   gr09->SetMarkerStyle(20);
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
      c1->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[i+1])+".pdf","pdf");

      if(pt_check){
	TCanvas* c4 = new TCanvas("c4");
	pt1_data_eff[i]->Draw();
	pt1_data_eff[i]->GetXaxis()->SetTitle("p_{t}^{1}");
	pt1_data_eff[i]->GetYaxis()->SetTitle("#epsilon");      
	c4->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[i+1])+"_pt1.pdf","pdf");
	
	TCanvas* c3 = new TCanvas("c3");
	pt2_data_eff[i]->Draw();
	pt2_data_eff[i]->GetXaxis()->SetTitle("p_{t}^{2}");
	pt2_data_eff[i]->GetYaxis()->SetTitle("#epsilon");      
	c3->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJet"+to_string(triggerVal[i+1])+"_pt2.pdf","pdf");	
      }

  }
    

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
  
  const double triggerVal_noLow_HF[n_trigger_HF-1] = {80,140,200,260,320,400,450,500};
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
      cHF1->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[i+1])+"_HF.pdf","pdf");

      if(pt_check){
	TCanvas* cHF4 = new TCanvas("cHF4");
	pt1_data_eff_HF[i]->Draw();
	pt1_data_eff_HF[i]->GetXaxis()->SetTitle("p_{t}^{1}");
	pt1_data_eff_HF[i]->GetYaxis()->SetTitle("#epsilon");      
	cHF4->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[i+1])+"_pt1_HF.pdf","pdf");
	
	TCanvas* cHF3 = new TCanvas("cHF3");
	pt2_data_eff_HF[i]->Draw();
	pt2_data_eff_HF[i]->GetXaxis()->SetTitle("p_{t}^{2}");
	pt2_data_eff_HF[i]->GetYaxis()->SetTitle("#epsilon");      
	cHF3->Print(CorrectionObject::_outpath+"plots/thresholds/"+"HLT_PFJetAve"+to_string(triggerVal_HF[i+1])+"_pt2_HF.pdf","pdf");	
      }

  }

    }
    



    
}
    
    


