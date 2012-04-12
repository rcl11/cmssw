
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

TH1* getHistogram(TFile* inputFile, const std::string& directory, const std::string& histogramName)
{  
  TString histogramName_full = directory.data();
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName.data());

  TH1* histogram = (TH1*)inputFile->Get(histogramName_full.Data());
  if ( !histogram) 
    std::cerr << "Failed to load histogram = " << histogramName_full 
	      << " from file = " << inputFile->GetName() << " !!" << std::endl;
  assert(histogram);

  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  histogram->Rebin(5);

  return histogram;
}

TH1* getHistogram(TFile* inputFile, const std::string& channel, double massPoint, 
		  const std::string& directory, const std::string& histogramName, double metResolution)
{
  std::string process_gg = "";
  if      ( massPoint >  95. ) process_gg = Form("ggHiggs%1.0f", massPoint);
  std::string process_qq = "";
  if      ( massPoint <  95. ) process_qq = "ZplusJets";
  else if ( massPoint < 135. ) process_qq = Form("vbfHiggs%1.0f", massPoint);
  std::string process_bb = "";
  if      ( massPoint > 155. ) process_bb = Form("bbHiggs%1.0f", massPoint);
 
  // CV: ONLY FOR TESTING
  process_gg = "";
  process_qq = "ZplusJets";
  process_bb = "";
  // CV: FOR TESTING ONLY
 
  std::string metResolution_label;
  if ( metResolution > 0. ) metResolution_label = Form("pfMEtRes%1.0f", metResolution);
  else metResolution_label = "pfMEtResMC";

  std::vector<TH1*> histograms;
  if ( process_gg != "" ) {
    std::string directory_gg = 
      Form("DQMData/%s/%s/%s/%s", process_gg.data(), channel.data(), metResolution_label.data(), directory.data());
    histograms.push_back(getHistogram(inputFile, directory_gg, histogramName));
  }
  if ( process_qq != "" ) {
    std::string directory_qq = 
      Form("DQMData/%s/%s/%s/%s/plotEntryType1", process_qq.data(), channel.data(), metResolution_label.data(), directory.data());
    histograms.push_back(getHistogram(inputFile, directory_qq, histogramName));
  }
  if ( process_bb != "" ) {
    std::string directory_bb = 
      Form("DQMData/%s/%s/%s/%s", process_bb.data(), channel.data(), metResolution_label.data(), directory.data());
    histograms.push_back(getHistogram(inputFile, directory_bb, histogramName));
  }
  TH1* histogramSum = NULL;
  for ( std::vector<TH1*>::const_iterator histogram = histograms.begin();
	histogram != histograms.end(); ++histogram ) {
    if ( !histogramSum ) {
      std::string histogramSumName = std::string((*histogram)->GetName()).append("_summed");
      histogramSum = (TH1*)(*histogram)->Clone(histogramSumName.data());
    } else {
      histogramSum->Add(*histogram);
    }
  }

  assert(histogramSum);

  if ( !histogramSum->GetSumw2N() ) histogramSum->Sumw2();
  if ( histogramSum->Integral() > 0. ) histogramSum->Scale(1./histogramSum->Integral());
  
  return histogramSum;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		    std::vector<std::string>& labelTextLines, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		    double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  histogram1->SetTitle("");
  histogram1->SetStats(false);
  histogram1->SetMinimum(yMin);
  histogram1->SetMaximum(yMax);

  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);

  int colors[4] = { 1, 2, 3, 4 };
  int lineStyles[4] = { 1, 7, 4, 3 };

  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->Draw("hist");

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->Draw("histsame");
  }

  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);
  legend->AddEntry(histogram1, legendEntry1.data(), "l");
  if ( histogram2 ) legend->AddEntry(histogram2, legendEntry2.data(), "l");
  legend->Draw();

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
  for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
	labelTextLine != labelTextLines.end(); ++labelTextLine ) {
    label->AddText(labelTextLine->data());
  }
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(labelTextSize);
  label->Draw();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete label;
  delete legend;
  delete canvas;  
}

struct histogram_vs_X_Type
{
  histogram_vs_X_Type(TH1* histogram, double x, double xErrUp, double xErrDown)
    : histogram_(histogram),
      x_(x),
      xErrUp_(xErrUp),
      xErrDown_(xErrDown)
  {}
   ~histogram_vs_X_Type() {}
  TH1* histogram_;
  double x_;
  double xErrUp_;
  double xErrDown_;
};

TGraph* makeGraph(std::vector<histogram_vs_X_Type>& histograms_vs_X, double y_true, const std::string& mode)
{
  enum { kResponse, kResolution };
  int mode_int = -1;
  if      ( mode == "response"   ) mode_int = kResponse;
  else if ( mode == "resolution" ) mode_int = kResolution;
  else assert(0);

  unsigned numPoints = histograms_vs_X.size();

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(numPoints);

  for ( unsigned iPoint = 0; iPoint < numPoints; ++iPoint ) {
    double x = histograms_vs_X[iPoint].x_;
    double xErrUp = histograms_vs_X[iPoint].xErrUp_;
    double xErrDown = histograms_vs_X[iPoint].xErrDown_;

    TH1* histogram = histograms_vs_X[iPoint].histogram_;

    double histogram_mean = histogram->GetMean();
    double histogram_meanErr = histogram->GetMeanError();
    double histogram_rms = histogram->GetRMS();
    double histogram_rmsErr = histogram->GetRMSError();
    
    double y, yErr;
    if ( mode_int == kResponse ) {
      y = histogram_mean/y_true;
      yErr = histogram_meanErr/y_true;
    } else if ( mode_int == kResolution ) {
      y = histogram_rms/histogram_mean;
      yErr = y*TMath::Sqrt(TMath::Power(histogram_rms/histogram_rmsErr, 2.) + TMath::Power(histogram_mean/histogram_meanErr, 2.));
    } else assert(0);

    graph->SetPoint(iPoint, x, y);
    graph->SetPointError(iPoint, xErrDown, xErrUp, yErr, yErr);
  }

  return graph;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
		std::vector<std::string>& labelTextLines, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  int colors[4] = { 1, 2, 3, 4 };
  int markerStyles[4] = { 20, 21, 22, 23 };

  graph1->SetLineColor(colors[0]);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->Draw("p");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->Draw("p");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->Draw("p");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->Draw("p");
  }
  
  TLegend* legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(legendTextSize);
  legend->AddEntry(graph1, legendEntry1.data(), "l");
  if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "l");
  if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "l");
  if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "l");
  legend->Draw();

  TPaveText* label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
  for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
	labelTextLine != labelTextLines.end(); ++labelTextLine ) {
    label->AddText(labelTextLine->data());
  }
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextColor(1);
  label->SetTextAlign(12);
  label->SetTextSize(labelTextSize);
  label->Draw();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete dummyHistogram;
  delete label;
  delete legend;
  delete canvas;  
}

void makeSVfitPerformancePlots()
{
  std::string inputFileName = "/data1/veelken/tmp/svFitStudies/AHtautau/2012Apr09/svFitPerformanceAnalysisPlots_all_2012Mar13.root";
  TFile* inputFile = new TFile(inputFileName.data());

  gROOT->SetBatch(true);

  std::string directory_PSkine_woLogM_Int      = "nSVfitAnalyzerOption1b";
  std::string directory_PSkine_woLogM_Fit      = "nSVfitAnalyzerOption1a";
  std::string directory_PSkine_wLogM_Int       = "nSVfitAnalyzerOption0b";
  std::string directory_PSkine_wLogM_Fit       = "nSVfitAnalyzerOption0a";
  std::string directory_MCkine_all_Int         = "nSVfitAnalyzerOption4b";
  std::string directory_MCkine_selected_Int    = "nSVfitAnalyzerOption5b";
  std::string directory_MEkine1_woPol_Int      = "nSVfitAnalyzerOption6b";
  std::string directory_MEkine1_woPol_Fit      = "nSVfitAnalyzerOption6a";
  std::string directory_MEkine12_wPolZorAH_Int = "nSVfitAnalyzerOption12b";
  std::string directory_MEkine12_wPolZorAH_Fit = "nSVfitAnalyzerOption12a";
  
  std::string histogramName_svFitMass = "svFitMass";
  std::string xAxisTitle_svFitMass    = "M_{#tau#tau} / GeV";
  std::string histogramName_dPhi12    = "dPhi12";

  typedef std::pair<double, double> pdouble;
  std::vector<pdouble> xRanges_dPhi;
  xRanges_dPhi.push_back(pdouble(  0.,  30.));
  xRanges_dPhi.push_back(pdouble( 30.,  60.));
  xRanges_dPhi.push_back(pdouble( 60.,  90.));
  xRanges_dPhi.push_back(pdouble( 90., 120.));
  xRanges_dPhi.push_back(pdouble(120., 140.));
  xRanges_dPhi.push_back(pdouble(140., 160.));
  xRanges_dPhi.push_back(pdouble(160., 170.));
  xRanges_dPhi.push_back(pdouble(170., 175.));
  xRanges_dPhi.push_back(pdouble(175., 180.));

  std::vector<double> massPoints;
  massPoints.push_back(90.);
  massPoints.push_back(120.);
  massPoints.push_back(130.);
  massPoints.push_back(160.);
  massPoints.push_back(200.);
  massPoints.push_back(300.);
  massPoints.push_back(450.);

  std::vector<double> massPointsToCompare;
  massPoints.push_back(90.);
  massPoints.push_back(130.);
  massPoints.push_back(300.);

  //double metResolution_nominal = -1.;
  double metResolution_nominal = 25.;

  std::string channel = "muTau";

  std::vector<std::string> label_ZplusJets;
  label_ZplusJets.push_back(std::string("Z #rightarrow #tau#tau MC"));
  label_ZplusJets.push_back(std::string("M = 90 GeV"));

  TH1* histogram_ZplusJets = 
    getHistogram(inputFile, channel, 90., directory_PSkine_woLogM_Int, 
		 histogramName_svFitMass, metResolution_nominal);
  showHistograms(800, 800, 
		 histogram_ZplusJets, "PS model", NULL, "",
		 0.04, 0.61, 0.74, 0.28, 0.15,
		 label_ZplusJets, 0.04, 0.175, 0.78, 0.24, 0.11, 
		 0., 250., xAxisTitle_svFitMass, 1.2,
		 0., 0.165, "a.u.", 1.6,
		 "svFitPerformance_ZplusJets_PSkine_woLogM_Int.eps");

  for ( std::vector<pdouble>::const_iterator xRange_dPhi = xRanges_dPhi.begin();
	xRange_dPhi != xRanges_dPhi.end(); ++xRange_dPhi ) {
    double dPhi_min = xRange_dPhi->first;
    double dPhi_max = xRange_dPhi->second;
    std::string dPhi_label;
    if      ( dPhi_min >  0. && dPhi_max <  180. ) dPhi_label = Form("dPhi%1.0fto%1.0f", dPhi_min, dPhi_max);
    else if ( dPhi_min == 0. && dPhi_max <  180. ) dPhi_label = Form("dPhiLt%1.0f", dPhi_max);
    else if ( dPhi_min >  0. && dPhi_max == 180. ) dPhi_label = Form("dPhiGt%1.0f", dPhi_min);
    else assert(0);
    std::string histogramName_svFitMass_full = Form("%s/%s", dPhi_label.data(), histogramName_svFitMass.data());
    TH1* histogram_ZplusJets = 
      getHistogram(inputFile, channel, 90., directory_PSkine_woLogM_Int, 
		   histogramName_svFitMass_full, metResolution_nominal);
    std::string dPhi_label2;
    if      ( dPhi_min >  0. && dPhi_max <  180. ) dPhi_label2 = Form("%1.0f < #Delta#phi < %1.0f", dPhi_min, dPhi_max);
    else if ( dPhi_min == 0. && dPhi_max <  180. ) dPhi_label2 = Form("#Delta#phi < %1.0f", dPhi_max);
    else if ( dPhi_min >  0. && dPhi_max == 180. ) dPhi_label2 = Form("#Delta#phi > %1.0f", dPhi_min);
    else assert(0);
    std::vector<std::string> label_ZplusJets_full = label_ZplusJets;
    label_ZplusJets_full.push_back(dPhi_label2);
    std::string outputFileName = Form("svFitPerformance_ZplusJets_PSkine_woLogM_Int_%s.eps", dPhi_label.data());
    showHistograms(800, 800, 
		   histogram_ZplusJets, "PS model", NULL, "",
		   0.04, 0.61, 0.74, 0.28, 0.15,
		   label_ZplusJets_full, 0.04, 0.175, 0.725, 0.24, 0.165, 
		   0., 250., xAxisTitle_svFitMass, 1.2,
		   0., 0.165, "a.u.", 1.6,
		   outputFileName.data());
  }

  std::vector<histogram_vs_X_Type> histograms_ZplusJets_vs_dPhi;
  for ( std::vector<pdouble>::const_iterator xRange_dPhi = xRanges_dPhi.begin();
	xRange_dPhi != xRanges_dPhi.end(); ++xRange_dPhi ) {
    double dPhi_min = xRange_dPhi->first;
    double dPhi_max = xRange_dPhi->second;
    std::string dPhi_label;
    if      ( dPhi_min >  0. && dPhi_max <  180. ) dPhi_label = Form("dPhi%1.0fto%1.0f", dPhi_min, dPhi_max);
    else if ( dPhi_min == 0. && dPhi_max <  180. ) dPhi_label = Form("dPhiLt%1.0f", dPhi_max);
    else if ( dPhi_min >  0. && dPhi_max == 180. ) dPhi_label = Form("dPhiGt%1.0f", dPhi_min);
    else assert(0);
    std::string histogramName_dPhi12_full = Form("%s/%s", dPhi_label.data(), histogramName_dPhi12.data());
    TH1* histogram_ZplusJets_dPhi12 = 
      getHistogram(inputFile, channel, 90., directory_PSkine_woLogM_Int, 
		   histogramName_dPhi12_full, metResolution_nominal);
    double dPhi_mean = histogram_ZplusJets_dPhi12->GetMean();
    std::string histogramName_svFitMass_full = Form("%s/%s", dPhi_label.data(), histogramName_svFitMass.data());
    TH1* histogram_ZplusJets_svFitMass = 
      getHistogram(inputFile, channel, 90., directory_PSkine_woLogM_Int, 
		   histogramName_svFitMass_full, metResolution_nominal);
    histograms_ZplusJets_vs_dPhi.push_back(
      histogram_vs_X_Type(histogram_ZplusJets_svFitMass, dPhi_mean, dPhi_max - dPhi_mean, dPhi_mean - dPhi_min));
  }    
  std::vector<std::string> label_dPhi;
  TGraph* graph_ZplusJets_response = makeGraph(histograms_ZplusJets_vs_dPhi, 90., "response");
  showGraphs(800, 600,
	     graph_ZplusJets_response, "M = 90 GeV", NULL, "", NULL, "", NULL, "",
	     0.04, 0.61, 0.74, 0.28, 0.15,
	     label_dPhi, 0.04, 0.175, 0.725, 0.24, 0.165, 
	     0., 180., "#Delta#phi / #circ", 1.2,
	     0., 0.50, "<M_{#tau#tau}>/M", 1.6,
	     "svFitPerformance_ZplusJets_PSkine_woLogM_Int_response_vs_dPhi.eps");
  TGraph* graph_ZplusJets_resolution = makeGraph(histograms_ZplusJets_vs_dPhi, 90., "resolution");
  showGraphs(800, 600,
	     graph_ZplusJets_resolution, "M = 90 GeV", NULL, "", NULL, "", NULL, "",
	     0.04, 0.61, 0.74, 0.28, 0.15,
	     label_dPhi, 0.04, 0.175, 0.725, 0.24, 0.165, 
	     0., 180., "#Delta#phi / #circ", 1.2,
	     0., 0.50, "#sigmaM_{#tau#tau}/<M_{#tau#tau}>", 1.6,
	     "svFitPerformance_ZplusJets_PSkine_woLogM_Int_resolution_vs_dPhi.eps");

  delete inputFile;
}
