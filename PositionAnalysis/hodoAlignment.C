
// MACRO TO CALCULATE ALIGNMENT WITH HODOSCOPE - run positionAnalysis for all runs you want to analyze before

#include <iostream>
using namespace std;

float getpos(int run, int refrun, float refx, float runspacing){
  return refx+(run-refrun)*runspacing;
}

void hodoAlignment(bool dox, int start, int end, float runspacing, int refrun, float refpos, float left_cef3_corr_selection=2200, float right_cef3_corr_selection=3800){

  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(kFALSE);


  const int nfibers=8;
  const float shift = refpos - refrun*runspacing;

  TH1F *histo[nfibers];
  TH1F *histoall[nfibers];
  float hodo_corr[nfibers];
  TBranch *b_hodo_corr;
  float cef3_corr[4];
  TBranch *b_cef3_corr;
  float results[nfibers];
  float resultserr[nfibers];
  float sigmas[nfibers];
  float sigmaserr[nfibers];
  float integral[nfibers];
  TGraphErrors *position = new TGraphErrors(nfibers);
  TGraphErrors *beamRMS = new TGraphErrors(nfibers);
  TGraph *intercalibration = new TGraph(nfibers);

  float minpos = getpos(start,refrun,refpos,runspacing);
  float maxpos = getpos(end,refrun,refpos,runspacing);
  int runs = end-start+1;

for (int fib=0; fib<nfibers; fib++){
  histo[fib]=new TH1F(Form("histo%d",fib),Form("histo%d",fib),runs,minpos-runspacing/2,maxpos+runspacing/2);
  histoall[fib]=(TH1F*)(histo[fib]->Clone(Form("histoall%d",fib)));

    for (int i=start; i<=end; i++){
            system(Form("ln -sf PosAn_BTF_%d_*.root .templink",i));
            TFile *f = new TFile(".templink","read");
            TTree *t;
            f->GetObject("tree_passedEvents",t);
            t->SetBranchAddress(dox ? "hodox_corr" : "hodoy_corr", hodo_corr, &b_hodo_corr);
            t->SetBranchAddress("cef3_corr", cef3_corr, &b_cef3_corr);
            int entries = t->GetEntries();
	    float pos = getpos(i,refrun,refpos,runspacing);
            for (int j=0; j<entries; j++){
	      t->GetEntry(j);
	      float cef3=0;
	      for (int l=0; l<4; l++) cef3+=cef3_corr[l];
	      if (cef3>right_cef3_corr_selection) continue;
	      if (cef3<left_cef3_corr_selection) continue;
	      if (hodo_corr[fib]>0) histo[fib]->Fill(pos);
	      histoall[fib]->Fill(pos);
	    }
    }
    histo[fib]->Divide(histo[fib],histoall[fib],1,1,"B");
    
    TF1 *func = new TF1("model","[0]+[1]*TMath::Exp(-(x-[2])*(x-[2])/(2*[3]*[3]))",minpos-runspacing/2,maxpos+runspacing/2);
    func->SetParameters(1,1,0.5*(minpos+maxpos),1);
    func->SetParLimits(0,0,1000);
    func->SetParLimits(1,0,1000);
    func->SetParLimits(2,minpos,maxpos);
    func->SetParLimits(3,1.5,5);
    
    TCanvas *c = new TCanvas(Form("canv%d",fib),Form("canv%d",fib));
    c->cd();
    histo[fib]->Draw("E1");
    histo[fib]->GetYaxis()->SetRangeUser(0,1);
    func->SetParameter(2,histo[fib]->GetMean());
    histo[fib]->Fit(func);
    results[fib]=func->GetParameter(2);
    resultserr[fib]=func->GetParError(2);
    sigmas[fib]=func->GetParameter(3);
    sigmaserr[fib]=func->GetParError(3);
    integral[fib]=func->GetParameter(1)*func->GetParameter(3);
    position->SetPoint(fib,fib,results[fib]);
    position->SetPointError(fib,0,resultserr[fib]);
    beamRMS->SetPoint(fib,fib,sigmas[fib]);
    beamRMS->SetPointError(fib,0,sigmaserr[fib]);

 }

 TCanvas *c2 = new TCanvas("c2","c2");
 c2->cd();
 TH2F *dummy = new TH2F("align",Form("Beam alignment with hodoscope - %s direction",(dox ? "X" : "Y")),1,-0.5,nfibers-0.5,1,minpos,maxpos);
 dummy->GetXaxis()->SetTitle("Fiber");
 dummy->GetYaxis()->SetTitle("Position (mm)");
 dummy->Draw();
 position->SetLineWidth(2);
 position->Draw("P");

 TF1 *lin = new TF1("lin","[0]+x*[1]",-0.5,nfibers-0.5);
 position->Fit(lin);
 float cent = lin->Eval(0.5*(nfibers-1));
 cout << "Measured center position: " << cent << endl;

 TLine *line = new TLine(-0.5,cent,nfibers-0.5,cent);
 line->SetLineWidth(2);
 line->SetLineColor(kBlue);
 line->Draw();

 c2->SaveAs(Form("plot_hodoAlignment_%s.root",((dox)?"X":"Y")));
 c2->SaveAs(Form("plot_hodoAlignment_%s.pdf" ,((dox)?"X":"Y")));
 c2->SaveAs(Form("plot_hodoAlignment_%s.png" ,((dox)?"X":"Y")));




 TCanvas *c3 = new TCanvas("c3","c3");
 c3->cd();
 TH2F *dummy3 = new TH2F("myrms",Form("Beam RMS with hodoscope - %s direction",(dox ? "X" : "Y")),1,-0.5,nfibers-0.5,1,0,10);
 dummy3->GetXaxis()->SetTitle("Fiber");
 dummy3->GetYaxis()->SetTitle("Beam RMS (mm)");
 dummy3->Draw();
 beamRMS->SetLineWidth(2);
 beamRMS->Draw("P"); 

 TF1 *lin3 = new TF1("lin3","[0]",-0.5,nfibers-0.5);
 beamRMS->Fit(lin3);
 cout << "Measured beam RMS: " << lin3->GetParameter(0) << endl;

 TLine *line3 = new TLine(-0.5,lin3->GetParameter(0),nfibers-0.5,lin3->GetParameter(0));
 line3->SetLineWidth(2);
 line3->SetLineColor(kBlue);
 line3->Draw();


 c3->SaveAs(Form("plot_beamRMS_%s.root",((dox)?"X":"Y")));
 c3->SaveAs(Form("plot_beamRMS_%s.pdf" ,((dox)?"X":"Y")));
 c3->SaveAs(Form("plot_beamRMS_%s.png" ,((dox)?"X":"Y")));





 TCanvas *c4 = new TCanvas("c4","c4");
 c4->cd();
 TH2F *dummy4 = new TH2F("myintercalib",Form("Hodoscope fiber efficiency - %s direction",(dox ? "X" : "Y")),1,-0.5,nfibers-0.5,1,0,2);
 dummy4->GetXaxis()->SetTitle("Fiber");
 dummy4->GetYaxis()->SetTitle("Efficiency rel. to average");
 dummy4->Draw();
 float sum=0;
 for (int fib=0; fib<nfibers; fib++) sum+=integral[fib];
 sum/=nfibers;
 for (int fib=0; fib<nfibers; fib++) intercalibration->SetPoint(fib,fib,integral[fib]/sum);
 intercalibration->SetMarkerStyle(20);
 intercalibration->Draw("P");

 c4->SaveAs(Form("plot_hodoEfficiency_%s.root",((dox)?"X":"Y")));
 c4->SaveAs(Form("plot_hodoEfficiency_%s.pdf" ,((dox)?"X":"Y")));
 c4->SaveAs(Form("plot_hodoEfficiency_%s.png" ,((dox)?"X":"Y")));

 cout << endl << "Hodoscope efficiency rel. to average, along " << (dox ? "X" : "Y") << endl;
 for (int fib=0; fib<nfibers; fib++) cout << "Fiber " << fib << " = " << integral[fib]/sum << endl;


}
