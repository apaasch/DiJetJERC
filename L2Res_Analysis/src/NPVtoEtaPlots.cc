#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
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
#include <TH2D.h>

using namespace std;

void CorrectionObject::NPVtoEtaPlots(){
  cout << "--------------- StartingNPVtoEtaPlots() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  CorrectionObject::make_path(CorrectionObject::_outpath+"plots/control/");

  TH2D *hdata[n_pt-1];
  
  int count = 0;
  TString name1 = "hist_data_NPVtoEta_";

  for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
      TString name = name1  + pt_name; 
      hdata[k] = new TH2D(name,"",nResponseBins ,-3.4,-3.4,70, -0.1, 69.9);
      count++;
    }

  //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");  
  TTreeReaderValue<Float_t> probejet_pt_data(myReader_DATA, "probejet_pt");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<int> nvertices_data(myReader_DATA, "nvertices");  

  int myCount = 0;
  int myCount_cut = 0;
  while (myReader_DATA.Next()) {
    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;	  

      hdata[k]->Fill(*probejet_eta_data,*nvertices_data,*weight_data);
      myCount++;
	
    }
  }
  //DEBUG
  std::cout<<"\ncount data "<<myCount<<"  count cut data "<<myCount_cut<<std::endl<<std::endl;

  

  TFile* test_out_data = new TFile(CorrectionObject::_outpath+"plots/control/nvert_eta_2d_data.root","RECREATE");
    for(int k=0; k<n_pt-1; k++){
      hdata[k]->Write();
    }
    std::cout<<"Histograms were written to root files\n";
  test_out_data->Close();
  delete test_out_data;



   TFile* f_rel_data_nvert = new TFile(CorrectionObject::_outpath+"plots/control/nvert_eta_2d_data.root","READ");

    std::cout<<"Histograms were reread from root files\n";
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 

    
    for(int j=0; j<n_pt-1; j++){
      TCanvas* cFullA_nvert = new TCanvas("c1","c1",800,600);
      // tdrCanvas(cFullA_nvert,"cFullA_nvert",h,4,10,kSquare,CorrectionObject::_lumitag);
      TH2D* htemp_rel_data_nvert;
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      std::cout<<"Histograms to pdfs loop at "<<pt_name<<std::endl;
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_rel_data = "hist_data_NPVtoEta_"+pt_name;
      htemp_rel_data_nvert = (TH2D*)f_rel_data_nvert->Get(name_rel_data);
      htemp_rel_data_nvert->Draw("E");
      htemp_rel_data_nvert->GetXaxis()->SetTitle("#eta");
      htemp_rel_data_nvert->GetYaxis()->SetTitle("n primary vertices");
      htemp_rel_data_nvert->GetZaxis()->SetTitle("Entries per Bin");      
      htemp_rel_data_nvert->GetXaxis()->SetLimits(-3.4,3.4);
      htemp_rel_data_nvert->GetYaxis()->SetLimits(0,70);
      htemp_rel_data_nvert->Draw("COLZ");		
      tex->DrawLatex(0.47,0.85,"Data, ");
      tex->DrawLatex(0.54,0.8,legname);		
      cFullA_nvert->SaveAs(CorrectionObject::_outpath+"plots/control/nvert_eta_DATA_" + CorrectionObject::_generator_tag +"_pt_"+ pt_range[j] + "_" + pt_range[j+1] + ".pdf");
      delete cFullA_nvert;
      delete htemp_rel_data_nvert;
    }
 


  
  
}
