#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../atlasstyle-00-04-02/AtlasUtils.h"
#include "../atlasstyle-00-04-02/AtlasStyle.h"
#include "../atlasstyle-00-04-02/AtlasLabels.h"
#include "../atlasstyle-00-04-02/AtlasStyle.C"

#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "../atlasstyle-00-04-02/AtlasLabels.C"
#include "../atlasstyle-00-04-02/AtlasUtils.C"
#endif

#ifdef __CINT__
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasLabels.C");
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasUtils.C");
#endif
const Float_t Weight[] = {6.7890E+07 * 6.1692E-06, 6.3996E+05 * 5.8420E-05, 4.7192E+03 * 1.1270E-04, 9.2038E-05 * 2.6602E+01};
Float_t FCal_range[] = {0, 0.063719, 0.14414, 0.289595, 0.525092, 0.87541, 1.36875, 2.04651, 2.98931, 5}; // fcal_cuts options

const int drbin = 20;
const int nFlav = 3;

const Float_t eta_selection = 2.1;
const bool unique_B = false;
const bool jetTruth = true;

char legend1[nFlav][50] = {"b-Jets", "c-jet", "light jet"};
char leg[nFlav][10] = {"B", "C", "U"};
const float min_distx = -0.5;
const float max_distx = -0.25;
const float min_disty = -0.95;
const float max_disty = -0.75;
const float min_distz = -200;
const float max_distz = 200;

const int myColor[] = {kBlue, kViolet, kMagenta, 418, kOrange, kYellow, kRed, kSpring,  kTeal, kCyan, kAzure, kGray, kGray + 1, kGray + 3};
const float min_distxs = -20;
const float max_distxs = 20;
const float min_distys = -20;
const float max_distys = 20;
const float min_distzs = -200;
const float max_distzs = 200;
const int cet[] = {0, 2, 2, 5, 5, 8}; //selected centrality sections
const int cet_N = (sizeof(cet) / sizeof(int)) / 2;

const int nData = 4;
const int dist_bin3 = 100;
const char *dataType = "WorkingDefault";
const char *dataTitle = "inclusive dijet samples";
char data[nData + 1][50] = {"Retrained pT > 0.5 GeV", "Retrained pT > 1.5 GeV", "Retrained pT > 2.0 GeV", "Retrained pT > 4.0 GeV", "Original"};
const bool PbPb = false;
char Type[2][10] = {"pp", "PbPb"};

const float cutMax = 40;
const float cutMin = -40;
const int numCuts = 200;
const float cutInt = (cutMax - cutMin) / numCuts;

const float ptLim[2] = {50., 100.};
const float trkLim[nData] = {0.5, 1.5, 2.0, 4.0};
const int msize = 1;

void Draw_llr()
{
    const char *uniqueness = unique_B ? "_Unique_B" : "";
    const char *jetTruthness = jetTruth ? "" : "_allJet";
    TFile *file[nData][2];
    TH1::AddDirectory(kFALSE);
    gStyle->SetLineWidth(1);

    for (int d = 0; d < nData; d++)
    {
        for (int pt = 0; pt < 2; pt++)
        {
            file[d][pt] = TFile::Open(Form("../IPxD/%s%s%srapidity%.1f_tuningIPEval_%.2f_%.1f.root", Type[PbPb], dataType, Type[PbPb], eta_selection, ptLim[pt], trkLim[d]));
        }
    }
    //TFile *out = TFile::Open(Form("%srapidity%s%.1f_LLR_IP_tot.root", dataType, Type[PbPb], eta_selection), "RECREATE");

    TH1F *llr_u[cet_N][nData + 1][2];
    TGraph* llr_u_g[cet_N][nData + 1][2];
    TH1F *llr_b[cet_N][nData + 1][2];
    TGraph* llr_b_g[cet_N][nData + 1][2];
    TH1F *llr_c[cet_N][nData + 1][2];
    TGraph* llr_c_g[cet_N][nData + 1][2];

    TH1F *total_u[cet_N][nData + 1][2];
    TH1F *total_b[cet_N][nData + 1][2];
    TH1F *total_c[cet_N][nData + 1][2];

    //TGraph *b[cet_N][nData + 1][2];
    //TMultiGraph *bllr[cet_N][nData + 1][2];
    //TMultiGraph *bpt[cet_N][nData + 1];
    //TMultiGraph *bdata[cet_N][2];

    int cent_N = PbPb ? cet_N : 1;
    //cout << "here" << endl;
    for (int d = 0; d < nData + 1; d++)
    {
        for (int pt = 0; pt < 2; pt++)
        {
            //int cent_N = PbPb ? cet_N : 1;
            //int cent_N = 1;
            for (int i = 0; i < cent_N; i++)
            {
                if (d == nData)
                {
                    llr_b[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("cent_%d_%s_bu", i, leg[1]))->Clone());
                    llr_b[i][d][pt]->SetMarkerColor(myColor[0]);
                    llr_b[i][d][pt]->SetLineColor(myColor[0]);
                    llr_b[i][d][pt]->SetLineStyle(1);
                    llr_b[i][d][pt]->SetLineWidth(2);
                    llr_b[i][d][pt]->SetMarkerStyle(msize);
                    llr_b[i][d][pt]->SetMarkerSize(1);

                    llr_c[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("cent_%d_%s_bu", i, leg[2]))->Clone());
                    llr_c[i][d][pt]->SetMarkerColor(myColor[3]);
                    llr_c[i][d][pt]->SetLineColor(myColor[3]);
                    llr_c[i][d][pt]->SetLineStyle(9);
                    llr_c[i][d][pt]->SetLineWidth(2);
                    llr_c[i][d][pt]->SetMarkerStyle(msize);
                    llr_c[i][d][pt]->SetMarkerSize(1);

                    llr_u[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("cent_%d_%s_bu", i, leg[0]))->Clone());
                    llr_u[i][d][pt]->SetMarkerColor(myColor[6]);
                    llr_u[i][d][pt]->SetLineColor(myColor[6]);
                    llr_u[i][d][pt]->SetLineStyle(3);
                    llr_u[i][d][pt]->SetLineWidth(2);
                    llr_u[i][d][pt]->SetMarkerStyle(msize);
                    llr_u[i][d][pt]->SetMarkerSize(1);

                    total_b[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("pT_cent_%d_%s", i, leg[1]))->Clone());
                    total_b[i][d][pt]->SetMarkerColor(myColor[0]);
                    total_b[i][d][pt]->SetLineColor(myColor[0]);
                    total_b[i][d][pt]->SetLineStyle(1);
                    total_b[i][d][pt]->SetLineWidth(2);
                    total_b[i][d][pt]->SetMarkerStyle(msize);
                    total_b[i][d][pt]->SetMarkerSize(1);

                    total_c[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("pT_cent_%d_%s", i, leg[2]))->Clone());
                    total_c[i][d][pt]->SetMarkerColor(myColor[3]);
                    total_c[i][d][pt]->SetLineColor(myColor[3]);
                    total_c[i][d][pt]->SetLineStyle(9);
                    total_c[i][d][pt]->SetLineWidth(2);
                    total_c[i][d][pt]->SetMarkerStyle(msize);
                    total_c[i][d][pt]->SetMarkerSize(1);

                    total_u[i][d][pt] = (TH1F *)(file[0][pt]->Get(Form("pT_cent_%d_%s", i, leg[0]))->Clone());
                    total_u[i][d][pt]->SetMarkerColor(myColor[6]);
                    total_u[i][d][pt]->SetLineColor(myColor[6]);
                    total_u[i][d][pt]->SetLineStyle(3);
                    total_u[i][d][pt]->SetLineWidth(2);
                    total_u[i][d][pt]->SetMarkerStyle(msize);
                    total_u[i][d][pt]->SetMarkerSize(1);
                }
                else
                {

                    llr_b[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("ret_cent_%d_%s_bu", i, leg[1]))->Clone());
                    llr_b[i][d][pt]->SetMarkerColor(myColor[0]);
                    llr_b[i][d][pt]->SetLineColor(myColor[0]);
                    llr_b[i][d][pt]->SetLineStyle(1);
                    llr_b[i][d][pt]->SetLineWidth(2);
                    llr_b[i][d][pt]->SetMarkerStyle(msize);
                    llr_b[i][d][pt]->SetMarkerSize(1);

                    llr_c[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("ret_cent_%d_%s_bu", i, leg[2]))->Clone());
                    llr_c[i][d][pt]->SetMarkerColor(myColor[3]);
                    llr_c[i][d][pt]->SetLineColor(myColor[3]);
                    llr_c[i][d][pt]->SetLineStyle(9);
                    llr_c[i][d][pt]->SetLineWidth(2);
                    llr_c[i][d][pt]->SetMarkerStyle(msize);
                    llr_c[i][d][pt]->SetMarkerSize(1);

                    llr_u[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("ret_cent_%d_%s_bu", i, leg[0]))->Clone());
                    llr_u[i][d][pt]->SetMarkerColor(myColor[6]);
                    llr_u[i][d][pt]->SetLineColor(myColor[6]);
                    llr_u[i][d][pt]->SetLineStyle(3);
                    llr_u[i][d][pt]->SetLineWidth(2);
                    llr_u[i][d][pt]->SetMarkerStyle(msize);
                    llr_u[i][d][pt]->SetMarkerSize(1);

                    total_b[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("pT_cent_%d_%s", i, leg[1]))->Clone());
                    total_b[i][d][pt]->SetMarkerColor(myColor[0]);
                    total_b[i][d][pt]->SetLineColor(myColor[0]);
                    total_b[i][d][pt]->SetLineStyle(1);
                    total_b[i][d][pt]->SetLineWidth(2);
                    total_b[i][d][pt]->SetMarkerStyle(msize);
                    total_b[i][d][pt]->SetMarkerSize(1);

                    total_c[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("pT_cent_%d_%s", i, leg[2]))->Clone());
                    total_c[i][d][pt]->SetMarkerColor(myColor[3]);
                    total_c[i][d][pt]->SetLineColor(myColor[3]);
                    total_c[i][d][pt]->SetLineStyle(9);
                    total_c[i][d][pt]->SetLineWidth(2);
                    total_c[i][d][pt]->SetMarkerStyle(msize);
                    total_c[i][d][pt]->SetMarkerSize(1);

                    total_u[i][d][pt] = (TH1F *)(file[d][pt]->Get(Form("pT_cent_%d_%s", i, leg[0]))->Clone());
                    total_u[i][d][pt]->SetMarkerColor(myColor[6]);
                    total_u[i][d][pt]->SetLineColor(myColor[6]);
                    total_u[i][d][pt]->SetLineStyle(3);
                    total_u[i][d][pt]->SetLineWidth(2);
                    total_u[i][d][pt]->SetMarkerStyle(msize);
                    total_u[i][d][pt]->SetMarkerSize(1);
                }
            }
        }
    }

    cout << "here" << endl;
    //out->Close();
    for (int d = 0; d < nData; d++)
    {
        for (int pt = 0; pt < 2; pt++)
        {
            file[d][pt]->Close();
        }
    }
    TCanvas *c0 = new TCanvas("c0", "c0",400,500);
    //bllr
    for (int c = 0; c < cent_N; c++)
    {
        for (int d = 0; d < nData + 1; d++)
        {
            for (int pt = 0; pt < 2; pt++)
            {
                //int cent_N = PbPb ? cet_N : 1;
                //int cent_N = 1;
                TPad *p0 = (TPad *)c0->cd();
                TH1F *h0 = (TH1F *)p0->DrawFrame(-20, 1.5*1e-4, 40, 1000);
                h0->GetXaxis()->SetTitle("IP2D log(P_{b}/P_{u})");
                h0->GetYaxis()->SetTitle("Unnormalized Weigthed Counts");
                h0->SetTitle(Form("LLR for IP2D Tagger %s", dataType));
                h0->Draw();

                gPad->SetTicks();
                gPad->SetGrid();
                c0->SetLogy();
                const char *Centrality = PbPb ? Form(" %d %% %d %%", 10 * cet[2 * c], 10 * cet[2 * c + 1]) : "";
                myText(0.4, 0.85, kBlack, Form("%s %s", Type[PbPb], Centrality), 0.04);
                myText(0.4, 0.8, kBlack, Form("pT_jet > %.0f GeV", ptLim[pt]), 0.04);
                myText(0.4, 0.75, kBlack, dataTitle, 0.04);
                llr_u[c][d][pt]->Draw("SAME hist");
                llr_b[c][d][pt]->Draw("SAME hist");
                llr_c[c][d][pt]->Draw("SAME hist");
                for (int f = 0; f < 3; f++)
                {
                    //myBoxText(0.7, 0.65 + 0.05 * f, 0.1, myColor[3 * f], 0, legend1[f], 0.3, myColor[3 * f], 1, false, 0.025);
                    int lsize = -1;
                    if (f == 0) lsize = 1;
                    if (f == 1) lsize = 10;
                    if (f == 2) lsize = 3;
                    myBoxText(0.6, 0.55 + f * 0.05, 0.05, 0, 0, legend1[f], lsize, myColor[3 * f], 1, true,0.05,1000,true,false);
                }
                c0->SaveAs(Form("Unnormalized_llr_%s_%s_%s_%s_%.2f.pdf", dataType, data[d], Type[PbPb], Centrality, ptLim[pt]));

h0 = (TH1F *)p0->DrawFrame(-20, 3*1e-7, 40, 10);
                h0->GetXaxis()->SetTitle("IP2D log(P_{b}/P_{u})");
                h0->GetYaxis()->SetTitle("Normalized Weigthed Fraction");
                h0->SetTitle(Form("LLR for IP2D Tagger %s", dataType));
                h0->Draw();

                gPad->SetTicks();
                gPad->SetGrid();
                c0->SetLogy();
                //const char *Centrality = PbPb ? Form(" %d %% %d %%", 10 * cet[2 * c], 10 * cet[2 * c + 1]) : "";
                myText(0.4, 0.85, kBlack, Form("%s %s", Type[PbPb], Centrality), 0.04);
                myText(0.4, 0.8, kBlack, Form("pT_jet > %.0f GeV", ptLim[pt]), 0.04);
                myText(0.4, 0.75, kBlack, dataTitle, 0.04);
                llr_u[c][d][pt]->Scale(1./llr_u[c][d][pt]->GetSumOfWeights());
                llr_u[c][d][pt]->Draw("SAME hist");
                llr_b[c][d][pt]->Scale(1./llr_b[c][d][pt]->GetSumOfWeights());
                llr_b[c][d][pt]->Draw("SAME hist");
                llr_c[c][d][pt]->Scale(1./llr_c[c][d][pt]->GetSumOfWeights());
                llr_c[c][d][pt]->Draw("SAME hist");
                for (int f = 0; f < 3; f++)
                {
                    //myBoxText(0.7, 0.65 + 0.05 * f, 0.1, myColor[3 * f], 0, legend1[f], 0.3, myColor[3 * f], 1, false, 0.025);
                    int lsize = -1;
                    if (f == 0) lsize = 1;
                    if (f == 1) lsize = 10;
                    if (f == 2) lsize = 3;
                    myBoxText(0.6, 0.7 - f * 0.05, 0.05, 0, 0, legend1[f], lsize, myColor[3 * f], 1, true,0.05,1000,true,false);
                }
                c0->SaveAs(Form("Normalized_llr_%s_%s_%s_%s_%.2f.pdf", dataType, data[d], Type[PbPb], Centrality, ptLim[pt]));

            }
        }
    }
}