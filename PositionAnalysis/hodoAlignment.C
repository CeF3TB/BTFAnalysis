
// MACRO TO CALCULATE ALIGNMENT WITH HODOSCOPE - run positionAnalysis for all runs you want to analyze before

#include <iostream>
using namespace std;

float getpos(int run, int refrun, float refx, float runspacing){
  return refx+(run-refrun)*runspacing;
}

void hodoAlignment(bool dox, int start, int end, float runspacing, int refrun, float refpos){

  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(kFALSE);


  const int nfibers=8;
  const float shift = refpos - refrun*runspacing;

  TH1F *histo[nfibers];
  TH1F *histoall[nfibers];
  float hodo_corr[nfibers];
  TBranch *b_hodo_corr;
  float results[nfibers];
  float resultserr[nfibers];
  TGraphErrors *position = new TGraphErrors(nfibers);

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
            int entries = t->GetEntries();
	    float pos = getpos(i,refrun,refpos,runspacing);
            for (int j=0; j<entries; j++){
	      t->GetEntry(j);
	      if (hodo_corr[fib]>0) histo[fib]->Fill(pos);
	      histoall[fib]->Fill(pos);
	    }
    }
    histo[fib]->Divide(histo[fib],histoall[fib],1,1,"B");
    
    TF1 *func = new TF1("model","[0]*x+[1]*TMath::Exp(-(x-[2])*(x-[2])/(2*[3]*[3]))",minpos-runspacing/2,maxpos+runspacing/2);
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
    position->SetPoint(fib,fib,results[fib]);
    position->SetPointError(fib,0,resultserr[fib]);
    
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

}
