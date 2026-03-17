//! ////////////////////////////////////////////////////////////////////////////
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 9th, 2025                                           //   
//!                                                                           // 
//! @details: Plots the spectra produced by generateSpectra_icarus.C          //
//! Setting effPur to true will plot a split canvas, with the spectra at the  //
//! top and efficiency and purity at the bottom.                              //
//! ////////////////////////////////////////////////////////////////////////////
#pragma once

#include "TROOT.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"

//! To do common ROOT operations (e.g. center a histogram title/axis label):
#include "rootlogon.C"

#include <vector>      //! to use std::vector objects
#include <type_traits> //! to use std::is_same
#include <string>      //! to use C++ string objects

//! CAFAna objects used in this file:
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

//! Imports user-defined general tools developed over a research career:
#include "auxTools.h"

//! Specifics of the electron neutrino selection are located in this file:
#include "nuESelectionHelper_icarus.h" 

// -----------------------------------------------------------------------------
// Function declarations
//!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
TString getCellColorCode(const double val);

void printTableHeader(
    const std::string& outFile = "./outputs_tables/table_output.tex", 
    const bool includeOther = true,
    const int precision = 3
);

void printTableFooter(
    const std::string& outFile = "./outputs_tables/table_output.tex"
);

void printEventsLine(
    const std::string& cutName, 
    const float nNuE, 
    const float nNuMu, 
    const float nNC, 
    const float nCosmics, 
    const float nOther, 
    const float eff, 
    const float pur,
    const std::string& outFile = "./outputs_tables/table_output.tex",
    const int precision = 3, 
    const int fieldWidth = 6
);

void drawComponentsLegend(
    TH1* hNuE, 
    TH1* hNuMu, 
    TH1* hNC, 
    TH1* hCosmics, 
    TH1* hOther
);

inline Spectrum* LoadSpectrum(
    const std::string& input_file,
    const std::string& label,
    const std::string& type,
    bool CRTVeto
);

template <typename PlotT, typename SelT>
void plotSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const bool spillMode,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const bool CRTVeto = false,
    const int rank = 0,
    const int nRanks = 1,
    const bool logScale = false,
    const bool effPur = true
);

template <typename PlotT, typename SelT>
void plotSequentialSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const bool spillMode,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const bool CRTVeto = false,
    const int rank = 0,
    const int nRanks = 1,
    const bool logScale = false,
    const bool effPur = true
);

template <typename PlotT, typename SelT>
void plotCRTSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const int rank = 0,
    const int nRanks = 1,
    const bool logScale = false
);



/** @fn plotSpectra_icarus()
 * @param input [optional] Type of input sample used (e.g. nus = "just neutrino 
 *        events" and nucosmics = "neutrino events and cosmic events")
 * @param combo [optional] Indicates if building plots from multiple input files
 * @param rank [optional] Rank of process (used in batch processing)
 * @param nRanks [optional] Total number of processes (used in batch processing)
*/
void plotSpectra_icarus( const std::string inFile, 
                                const std::string outDir = "./outputs_plots",
                                const std::string input = "nus", 
                                const bool spillMode = false,
                                const int rank = 0,
                                const int nRanks = 1)
{
  //! @note These are the original names of the input files used in the
  //! generateSpectra_icarus.C script. Ensure that you are using the correct
  //! finsample so this plotting script output the plots you expect.
 
  //! Testing sample file:
  // const std::string inFile = "refactored_g4_icarus_step2_g4step2_detsim_"
  //                               "stage0_stage1_13629324_0.flat.caf"
  //                               "-e8a8c7bb-7ad3-46b2-a5a6-2496ffc7199f.root"; 

  //! @note Nominal samweb definition for this electron neutrino selection:
  //! icaruspro_production_v09_89_01_02p01_2024A_ICARUS_MC_Sys_NuCos_2024A_MC_Sys_NuCos_respunCV_2ndV_flatcaf

  plotSpectra( inFile, outDir, input, spillMode, vars_slice, sels_slice, false, 
    rank, nRanks); //! non-CRT-veto plots
  plotSequentialSpectra( inFile, outDir, input, spillMode, vars_slice, individual_cuts_slice, false,
    rank, nRanks);

  plotSpectra( inFile, outDir, input, spillMode, vars_slice, sels_slice, true,
    rank, nRanks); //! CRT veto plots
  plotSequentialSpectra( inFile, outDir, input, spillMode, vars_slice, individual_cuts_slice, true,
    rank, nRanks);
 
  plotCRTSpectra( inFile, outDir, input, crtvars_spill, crtsels_spill, 
    rank, nRanks);
}

//!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

/** @fn getCellColorCode()
 * @brief Returns a LaTeX color code for a table cell based on the input value.
 *
 * This function encodes a numeric value as a LaTeX `\cellcolor{}` directive.
 * The color represents the magnitude of the value:
 * - |val| < 0.05 → no color (empty string)
 * - 0.05 ≤ |val| < 0.10 → light green
 * - 0.10 ≤ |val| < 0.15 → light yellow
 * - |val| ≥ 0.15 → light red
 *
 * @param val: Numeric value used to determine the color code.
 * @return TString containing the corresponding LaTeX color command, or an 
 * empty string if |val| < 0.05.
 */
TString getCellColorCode(const double val)
{
  TString str = ""; //! initialize as empty string for absVal < 0.05 (see below)
  double absVal = abs( val);

  //! @note: color-coding starts at 0.05
  if( absVal >= 0.05 && absVal < 0.10) str = "\\cellcolor{green!25}";
  else if( absVal >= 0.10 && absVal < 0.15) str = "\\cellcolor{yellow!25}";
  else if( absVal>=0.15) str = "\\cellcolor{red!25}";
  return str;
}


/** @fn printTableHeader()
 * @brief Prints the LaTeX table header for an event summary table.
 *
 * This function initializes a LaTeX table structure in the specified output file.
 * The header format depends on whether the "Other" background category is included.
 *
 * @param outFile [optional]: Path to the output `.tex` file.
 * @param includeOther [optional]: Whether to include the "Other" background column in the table.
 * @param precision [optional]: Number of digits to print after the decimal point.
 */
void printTableHeader( 
  const std::string& outFile, 
  const bool includeOther,
  const int precision)
{
  ofstream outStream( outFile, std::ios::trunc); //! overwrite or start new file
  outStream << std::fixed << std::setprecision( precision);

  outStream << "\\begin{table}[H]\n";
  outStream << "\\centering\n";
  outStream << "\\resizebox{\\textwidth}{!}{\n";

  if( includeOther){ //! include "other backgrounds" entry in table
    outStream << "\\begin{tabular}{|l||c|c|c|c|c||c|c|}\n";
    outStream << "\\hline \n";
  
    outStream << "\\multicolumn{1}{|c||}{} & \\multicolumn{5}{c||}{Number of Entries (\\% of Total)} & \\multicolumn{2}{c|}{Integrated} \n";
    outStream << "\\hline \n";

    outStream << "\\multicolumn{1}{|c||}{Cut} & $\\nu_{e}$ CC & $\\nu_{\\mu}$ CC & NC & Cosmic & Other & Efficiency & Purity \n";
    outStream << "\\hline \n";
  }

  else {
    outStream << "\\begin{tabular}{|l||c|c|c|c||c|c|}\n";
    outStream << "\\hline \n";
  
    outStream << "\\multicolumn{1}{|c||}{} & \\multicolumn{4}{c||}{Number of interactions (\\% total)} & \\multicolumn{2}{c|}{Integrated} \n";
    outStream << "\\hline \n";

    outStream << "\\multicolumn{1}{|c||}{Cut} & $\\nu_{e}$ CC & $\\nu_{\\mu}$ CC & NC & Cosmic & Efficiency & Purity \n";
    outStream << "\\hline \n";
  }
}


/** @fn printTableFooter()
 * @brief Prints the closing lines of a LaTeX event summary table.
 *
 * This function finalizes the LaTeX table by writing the closing
 * `\end{tabular}` and `\end{table}` commands to the specified output file.
 *
 * @param outFile [optional]: Path to the output `.tex` file.
 */
void printTableFooter( const std::string& outFile)
{
  ofstream outStream( outFile, std::ios::app); //! append output to end of file

  outStream << "\\end{tabular}} \n";
  outStream << "\\end{table}";
  outStream << "\n\n\n";
}


/**
 * @brief Writes a single row of event counts and percentages to a LaTeX table.
 *
 * Each entry is written in the format "N (P%)", where `N` is the raw count
 * and `P` is the percentage of the total. The line also includes efficiency
 * and purity values at the end.
 *
 * @param cutName: Label for the current selection or cut.
 * @param nNuE: Number of ν_e CC events.
 * @param nNuMu: Number of ν_μ CC events.
 * @param nNC: Number of neutral current events.
 * @param nCosmics: Number of cosmic events.
 * @param nOther: Number of "Other" background events.
 * @param eff: Selection efficiency value.
 * @param pur: Selection purity value.
 * @param outFile [optional]: Path to the output `.tex` file.
 * @param precision [optional]: Number of digits to print after the decimal point.
 * @param fieldWidth [optional]: Field width for numerical output alignment.
 */
void printEventsLine( const std::string& cutName, 
    const float nNuE, 
    const float nNuMu, 
    const float nNC, 
    const float nCosmics, 
    const float nOther, 
    const float eff, 
    const float pur,
    const std::string& outFile,
    const int precision, 
    const int fieldWidth)
{
  ofstream outStream( outFile, std::ios::app); //! append to end of file
  outStream << std::fixed << std::setprecision(precision);

  float total = nNuE + nNuMu + nNC + nCosmics + nOther;

  float pNuE     = 100 * nNuE     / total;
  float pNuMu    = 100 * nNuMu    / total;
  float pNC      = 100 * nNC      / total;
  float pCosmics = 100 * nCosmics / total;
  float pOther   = 100 * nOther   / total;

  float pTotalBackground = 100 * (total - nNuE) / total;

  outStream << std::fixed << std::setw( fieldWidth) << std::setprecision( precision);
  outStream << cutName  << "&";
  outStream << nNuE     << " ("<<pNuE<<")"     << "&";
  outStream << nNuMu    << " ("<<pNuMu<<")"    << "&";
  outStream << nNC      << " ("<<pNC<<")"      << "&";
  outStream << nCosmics << " ("<<pCosmics<<")" << "&";
  outStream << nOther   << " ("<<pOther<<")"   << "&";
  outStream << eff << "&" << pur << "\n";
  outStream << "\\hline \n";
}


/** @fn drawComponentsLegend()
 * @brief Create a TLegend to identify components specified by parameters.
 *
 * Creates and renders a ROOT TLegend object showing labels for
 * various event categories (ν_e CC, ν_μ CC, NC, cosmics, other).
 * 
 * @param hNuE: Histogram for electron neutrino component
 * @param hNuMu: Histogram for muon neutrino component
 * @param hNC: Histogram for neutral current interaction component
 * @param hCosmics: Histogram for cosmic component
 * @param hOther: Histogram for components not list above
*/
void drawComponentsLegend(
    TH1* hNuE, 
    TH1* hNuMu, 
    TH1* hNC, 
    TH1* hCosmics, 
    TH1* hOther)
{
  TLegend *leg = new TLegend( 0.75, 0.675, 0.90, 0.85);

  leg->AddEntry( hNuE,     "#nu_{e} CC",   "l");
  leg->AddEntry( hNuMu,    "#nu_{#mu} CC", "l");
  leg->AddEntry( hNC,      "NC",           "l");
  leg->AddEntry( hCosmics, "Cosmics",      "l");
  leg->AddEntry( hOther,   "Other Bkg",    "l");

//  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->Draw();
}

/**
 * @brief Load a Spectrum from file with a standardized label format.
 *
 * Constructs the full object label in the form:
 *     <label>_type=<type>_<CRTVetoState>
 *
 * where <CRTVetoState> is either "kCRTHitVetoFD" or "noVetoApplied".
 *
 * @param input_file: input ROOT filename containing the Spectrum
 * @param label: base label used to identify the Spectrum object
 * @param type: string inserted into "_type=<type>_" in the label
 * @param CRTVeto: if true, adds "kCRTHitVetoFD"; otherwise "noVetoApplied"
 *
 * @return Spectrum: raw pointer to the loaded Spectrum object
 */
inline Spectrum* LoadSpectrum(const std::string& input_file,
                              const std::string& label,
                              const std::string& type,
                              bool CRTVeto) {
  const std::string fullLabel =
      label + "_type=" + type + "_" + (CRTVeto ? "kCRTHitVetoFD" : "noVetoApplied");

  return LoadFromFile<Spectrum>(input_file, fullLabel).release();
}


/** @fn plotSpectra()
  @brief Plot Spectra on TCanvas and save the plots as PDFs in outDir.

  @param inFile: Name of individual input file
  @param outDir: Path where the output plot PDFs will go
  @param spillMode: flag to indicate if you're plotting slice/spill Spectra
  @param vars: Variables you want to plot
  @param sels: Cuts on the variables to plot
  @param CRTVeto [optional]: Indicates if plotting CRT veto variables; default is false
  @param rank [optional]: Rank of process (used in batch processing)
  @param nRanks [optional]: Total number of processes (used in batch processing)
  @param logScale [optional]: Set y-axis log scale; default is true
  @param effPur [optional]: Flag to generate efficiency and purity graphs; default is true
*/
template <typename PlotT, typename SelT>
void plotSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const bool spillMode,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const bool CRTVeto,
    const int rank,
    const int nRanks,
    const bool logScale,
    const bool effPur) 
{
  //! Place LaTeX table header creation outside any loops. For more information
  //! see the @SummaryTable tag below.
  if( !CRTVeto) {
    printTableHeader("./outputs_tables/single-cut_summaryTable.tex");
    printTableHeader("./outputs_tables/N-1-cut_summaryTable.tex");
  }

  const size_t kNVar = vars.size();
  const size_t kNSel = sels.size();
  int total = kNVar * kNSel;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  //! Make a plot for each variable and cut.
  //! 
  //! @note there are more TH1Ds than Spectra (see details 
  //! below) so we do not detail at the Spectrum level.
  //!
  //! Definitions analyzers might find useful throughout this file:
  //! 1. hAllNuE: ALL of the electron neutrino events; used for POT accounting,
  //!       efficiency, and purity
  //! 2. hTrueSigNuE: effective truth signal definition (restrictions on 
  //!       containment, fiducial volume, and reconstructed shower energy) 
  //!       applied to the truth-level electron neutrino designation.
  //! 
  //! @note any effective truth signal definition TH1Ds are created to model how
  //! reconstruction efforts inherantly have some bias associated with them 
  //! given we don't know if we reconstructed all showers that "went through" 
  //! our detector.
  for(size_t iVar = 0; iVar < kNVar; ++iVar){
    Spectrum *sAllNuE = LoadFromFile<Spectrum>(inFile, 
      "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=kNoCut_type=NuECC_" \
        + (CRTVeto ? "kCRTHitVetoFD" : "noVetoApplied") ).release();

    TH1D *hAllNuE = sAllNuE->ToTH1( sAllNuE->POT()); 
    float iAllNuE = hAllNuE->Integral();

    //! ////////////////////////////////////////////////////////////////////////
    //!  ^
    //!  |  POT-accounting Spectra/TH1D only created once per variable, but
    //!  \- POT shouldn't change according to variable plotted anyway.
    //! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    //!  /- Spectra/TH1D that are created for every variable-selection 
    //!  |  combination.
    //!  V
    //! ////////////////////////////////////////////////////////////////////////

    for(size_t iSel = 0; iSel < kNSel; ++iSel, ++idx){
      if( idx < start || idx >= end) continue; //! for batch processing on grid 
      
      std::string cornerTag = sels[iSel].label + (CRTVeto ? " CRT Hit Veto Applied" : "");
;
      std::string label = "var=" + vars[iVar].label \
                            + "_unit=" + vars[iVar].unit \
                            + "_sel=" + sels[iSel].label;

      //! @note: TH1D hMyName is a copy of Spectrum sMyName if such a Spectrum exists.
      TH1D *hNuE,              //! electron neutrinos
           *hTrueSigNuE,       //! 'true signal' electron neutrinos [see above]
           *hNuMu,             //! muon neutrinos
           *hNC,               //! neutral current neutrinos
           *hCosmics,          //! in-time cosmics
           *hTotal,            //! anything with a reconstructed neutrino
           *hOther,            //! anything that's not NC, NuECC, NuMuCC, or in-time cosmic
           *hBackground;       //! anything that's not NuECC

      Spectrum *sNuE        = LoadSpectrum(inFile, label, "NuECC", CRTVeto);
      Spectrum *sTrueSigNuE = LoadSpectrum(inFile, label, "NuECCTrueSignal", CRTVeto);
      Spectrum *sNuMu       = LoadSpectrum(inFile, label, "NuMuCC", CRTVeto);        
      Spectrum *sNC         = LoadSpectrum(inFile, label, "NuNC", CRTVeto);        
      Spectrum *sCosmics    = LoadSpectrum(inFile, label, "NuCosmic", CRTVeto);        
      Spectrum *sTotal      = LoadSpectrum(inFile, label, "NuTotal", CRTVeto);        
          
      hNuE        = sNuE->ToTH1( sNuE->POT());    
      hTrueSigNuE = sTrueSigNuE->ToTH1( sTrueSigNuE->POT());
      hNuMu       = sNuMu->ToTH1( sNuMu->POT());   
      hNC         = sNC->ToTH1( sNC->POT());     
      hCosmics    = sCosmics->ToTH1( sCosmics->POT());
      hTotal      = sTotal->ToTH1( sTotal->POT());  

      hOther = (TH1D*)hTotal->Clone();
      hOther->Add(hNuE, -1);
      hOther->Add(hNuMu, -1);
      hOther->Add(hNC, -1);
      hOther->Add(hCosmics, -1);

      hBackground = (TH1D*)hTotal->Clone();
      hBackground->Add(hNuE, -1);

      //! Extract plotting info (e.g. axes names, titles, etc.) about plot:
      TitleInfo info = parseTitle( hNuE->GetXaxis()->GetTitle());
      std::string xTitle = info.unit.empty() ? info.variable : 
          Form("%s [%s]", info.variable.c_str(), info.unit.c_str());

      //! Make signal to background ratio distribution:
      TH1D *hSigBkgRatio = (TH1D*)hNuE->Clone();
      hSigBkgRatio->Divide( hBackground);

      //! Make efficiency and purity distributions:
      TH1D *hEff = getEffHistogram( hNuE, hBackground, hAllNuE);
      TH1D *hPur = getPurHistogram( hNuE, hBackground);

      //! Format histograms:
      hEff->SetTitle(""); hPur->SetTitle(""); //! strip titles

      formatHist( hNuE, colorNuE, lineNuE, 2);
      formatHist( hTrueSigNuE, colorNuE, lineTrueSigNuE, 2);
      formatHist( hNuMu, colorNuMu, lineNuMu, 2);
      formatHist( hNC, colorNC, lineNC, 2);
      formatHist( hCosmics, colorCosmics, lineCosmics, 2);
      formatHist( hOther, colorOther, lineOther, 3);
      formatHist( hEff, colorEff, kSolid, 2);
      formatHist( hPur, colorPur, kSolid, 2);

      CenterTitles( hNuE);
      CenterTitles( hTrueSigNuE);
      CenterTitles( hNuMu);
      CenterTitles( hNC);
      CenterTitles( hCosmics);
      CenterTitles( hOther);
      CenterTitles( hEff);
      CenterTitles( hPur);

      //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      //! @SummaryTable: calculations only needed for one variable!
      //!                it shouldn't matter if you test iVar == j for any j
      if( !CRTVeto && iVar == 0) { //! can effectively cut on any variable
        float iNuE        = hNuE->Integral();
        float iNuMu       = hNuMu->Integral();
        float iNC         = hNC->Integral();
        float iCosmics    = hCosmics->Integral();
        float iOther      = hOther->Integral();
        float iBackground = hBackground->Integral();
        float iTotal      = hTotal->Integral();

        float iEff = iNuE / iAllNuE;
        float iPur = iNuE / (iNuE + iBackground);
        
        if( iSel == 0) { //! kNoCut/kNoSpillCut cases added to both Summary Tables
          printEventsLine( sels[iSel].label, 
            iNuE, iNuMu, iNC, iCosmics, iOther, 0.0, 0.0, 
            "./outputs_tables/single-cut_summaryTable.tex");
          printEventsLine( sels[iSel].label, 
            iNuE, iNuMu, iNC, iCosmics, iOther, 0.0, 0.0, 
            "./outputs_tables/N-1-cut_summaryTable.tex");
        }
        else if( iSel < nCuts){ //! Single-cut Summary Table:
          printEventsLine( sels[iSel].label, 
            iNuE, iNuMu, iNC, iCosmics, iOther, iEff, iPur, 
            "./outputs_tables/single-cut_summaryTable.tex");
        }
        else if( iSel < kNSel - 1) { //! N-1-cut Summary Table:
          printEventsLine( sels[iSel].label,
            iNuE, iNuMu, iNC, iCosmics, iOther, iEff, iPur,
            "./outputs_tables/N-1-cut_summaryTable.tex");
        }
        else { //! full selection added to both Summary Tables
          printEventsLine( sels[iSel].label, 
            iNuE, iNuMu, iNC, iCosmics, iOther, iEff, iPur, 
            "./outputs_tables/single-cut_summaryTable.tex");
          printEventsLine( sels[iSel].label, 
            iNuE, iNuMu, iNC, iCosmics, iOther, iEff, iPur, 
            "./outputs_tables/N-1-cut_summaryTable.tex");
        }
      } //! end @SummaryTable
      //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

      //! Make single-pad TCanvas for plotting different event samples:
      TCanvas *cEvents = new TCanvas( vars[iVar].label.c_str(), 
                                      vars[iVar].label.c_str(),
                                      cXRes, cYRes);
      if( iVar == 0) { //! fake stacking for Number of Slices/Spills variable
        hNuMu->Add( hNuE);
        hNC->Add( hNuMu);
        hCosmics->Add( hNC);
        hOther->Add( hCosmics);

        fillWithDimmedColor( hNuE, 0);
        fillWithDimmedColor( hNuMu, 0);
        fillWithDimmedColor( hNC, 0);
        fillWithDimmedColor( hCosmics, 0);
        fillWithDimmedColor( hOther, 0);
      }

      //! Set y-axis scale:
      if( logScale) {
        gPad->SetLogy();
        hOther->GetYaxis()->SetRangeUser( 1.0, 
                100*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther }));
      }
      else {
        hOther->GetYaxis()->SetRangeUser( 0.0, 
                1.3*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther }));
      }

      //! We work with hOther as the "master" histogram. Since it is (hopefully)
      //! near-zero in every bin, we plot this histogram first and use it to
      //! set the axes titles.
      hOther->GetXaxis()->SetTitle( xTitle.c_str());
      hOther->GetYaxis()->SetTitle( (spillMode ? 
        "Number of Spills" : "Number of Slices") );

      //! Draw hNuE, hNuMu, hNC, hCosmics, and hOther on the same axis.
      hOther->Draw("hist same");
      hCosmics->Draw("hist same");
      hNC->Draw("hist same");
      hNuMu->Draw("hist same");
      hNuE->Draw("hist same");

      gPad->SetGrid(1, 1); //! Set grid on both X and Y axes

      //! Extract the POT to print onto output plots.
      //! @note: POT is an amalgamation of spills. This quantity is not changed
      //! when making selection cuts, separating by interaction type, or
      //! plotting different variables. The POTTag is thus the total amount of
      //! POT from the dataset processed.
      std::string POTTag = makePOTTag( sAllNuE->POT());

      //! Draw plaintext information on events TCanvas.
      DrawSigBkgIntText( hNuE, hBackground, POTTag, 0.04);
      drawComponentsLegend( hNuE, hNuMu, hNC, hCosmics, hOther);
      isSimulation( true);
      myCornerLabel( cornerTag.c_str());

      //! Save the events TCanvas
      cEvents->Update();
      std::string eventsPlotName = outDir + "/" + input + "_" + label \
          + (CRTVeto ? "_CRTHitVetoApplied" : "") \
          + (logScale ? "__logScale" : "") \
          + ".pdf";
      cEvents->Print( eventsPlotName.c_str());
      cEvents->Close();
     
      //! Make split-pad TCanvas for plotting different event samples and
      //! their efficiency and/or purity:
      //! @note: Top pad is variable plotted, bottom is efficiency/purity.
      if( effPur) {
        TCanvas *cEventsEffPur;
        TPad *padEvents, *padEffPur;
        splitCanvasInTwo( cEventsEffPur, padEvents, padEffPur);

        padEvents->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1); 

        //! Set y-axis log scale:
        if( logScale) {
          gPad->SetLogy(1); //! toggle y-axis log scale on
          hOther->GetYaxis()->SetRangeUser( 1.0, 
                100*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther } ));
        }
        else {
          gPad->SetLogy(0); //! toggle y-axis log scale off
          hOther->GetYaxis()->SetRangeUser( 0.0,
                1.3*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther } ));
        }

        //! ////////////////////////////////////////////////////////////////////
        //! Again we treat hOther as the "master" histogram. Since it's
        //! near-zero in every bin, we plot this histogram first and use it to
        //! set the axes titles.
        //! ////////////////////////////////////////////////////////////////////

        //! @note: Pads share x-axis, so there's no need to put an x-axis
        //!        title on the top pad.
        hOther->GetXaxis()->SetLabelSize(0); //! remove hOther x-axis labels
        hOther->GetXaxis()->SetTitle("");
        
        //! Draw hNuE, hNuMu, hNC, hCosmics, and hOther on the same axis.
        hOther->Draw("hist same");
        hCosmics->Draw("hist same");
        hNC->Draw("hist same");
        hNuMu->Draw("hist same");
        hNuE->Draw("hist same");

        //! Draw plaintext information on events TCanvas.
        DrawSigBkgIntText( hNuE, hBackground, POTTag, 0.04);
        drawComponentsLegend( hNuE, hNuMu, hNC, hCosmics, hOther);
        isSimulation( true);

        myCornerLabel( cornerTag.c_str());
        cEventsEffPur->Update();
        padEffPur->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogy(0); //! never have log scale for efficiency/purity

        //! @note: The TH2 constructor used takes the following args, in order:
        //!        name, title, nXBins, xLow, xUp, nYBins, yLow, yUp
        //! @note: TH1::GetSize() includes one underflow and one overflow bin
        TH2 *axEffPur = new TH2F("", "",
          hNuE->GetSize()-2, 
          hNuE->GetXaxis()->GetXmin(), 
          hNuE->GetXaxis()->GetXmax(),
          30, 0.0, 1.1);
        TH2 *axPur = new TH2F("", "",
          hNuE->GetSize()-2, 
          hNuE->GetXaxis()->GetXmin(), 
          hNuE->GetXaxis()->GetXmax(),
          30, 0.0, 1.1*hSigBkgRatio->GetMaximum());

        //! Set efficiency/purity (bottom) pad's x-axis titles:
        axEffPur->GetXaxis()->SetTitle( xTitle.c_str());
        axPur->GetXaxis()->SetTitle( xTitle.c_str());

        //! Strip y-axis titles for the efficiency/purity pad since it's inferred
        //! that efficiency/purity is a fraction on the interval [0, 1]
        axEffPur->GetYaxis()->SetTitle("");
        axPur->GetYaxis()->SetTitle("");

        if( iSel == 0) { //! no cut? then we don't care about efficiency
          axPur->Draw();
          hPur->Draw("hist");
          DrawPurLegend( hPur, "purity");
        }
        else {
          axEffPur->Draw();
          hEff->Draw("hist same");
          hPur->Draw("hist same");
          DrawIntEffPurLegend( hEff, "efficiency", hPur, "purity");
        }

        gPad->SetGrid(1, 1); //! Set grid on both X and Y axes
        cEventsEffPur->Update();
        gPad->RedrawAxis();

        std::string eventsEffPurPlotName = outDir + "/" + input + "_" + label \
            + (CRTVeto ? "_CRTHitVetoApplied" : "") \
            + (logScale ? "__logScale" : "") \
            + "_effPur.pdf";
        cEventsEffPur->Print( eventsEffPurPlotName.c_str());
        cEventsEffPur->Close();
      } //! end if effPur
    } //! end iVar
  } //! end iSel

  //! Place LaTeX table footer creation outside any loops. For more information
  //! see the @SummaryTable tag above.
  if( !CRTVeto) {
    printTableFooter("./outputs_tables/single-cut_summaryTable.tex");
    printTableFooter("./outputs_tables/N-1-cut_summaryTable.tex");
  }
}


/** @fn plotSequentialSpectra()
  @brief Plot Spectra on TCanvas and save the plots as PDFs in outDir.

  @param inFile: Name of individual input file
  @param outDir: Path where the output plot PDFs will go
  @param spillMode: flag to indicate if you're plotting slice/spill Spectra
  @param vars: Variables you want to plot
  @param sels: Cuts on the variables to plot
  @param CRTVeto [optional]: Indicates if plotting CRT veto variables; default is false
  @param rank [optional]: Rank of process (used in batch processing)
  @param nRanks [optional]: Total number of processes (used in batch processing)
  @param logScale [optional]: Set y-axis log scale; default is true
  @param effPur [optional]: Flag to generate efficiency and purity graphs; default is true
*/
template <typename PlotT, typename SelT>
void plotSequentialSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const bool spillMode,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const bool CRTVeto,
    const int rank,
    const int nRanks,
    const bool logScale,
    const bool effPur) 
{
  //! Place LaTeX table header creation outside any loops. For more information
  //! see the @SummaryTable tag below.
  if( !CRTVeto) {
    printTableHeader("./outputs_tables/sequential-cut_summaryTable.tex");
  }

  const size_t kNVar = vars.size();
  const size_t kNSel = sels.size();
  int total = kNVar * kNSel;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  //! Make a plot for each variable and cut.
  //! 
  //! @note there are more TH1Ds than Spectra (see details 
  //! below) so we do not detail at the Spectrum level.
  //!
  //! Definitions analyzers might find useful throughout this file:
  //! 1. hAllNuE: ALL of the electron neutrino events; used for POT accounting,
  //!       efficiency, and purity
  //! 2. hTrueSigNuE: effective truth signal definition (restrictions on 
  //!       containment, fiducial volume, and reconstructed shower energy) 
  //!       applied to the truth-level electron neutrino designation.
  //! 
  //! @note any effective truth signal definition TH1Ds are created to model how
  //! reconstruction efforts inherantly have some bias associated with them 
  //! given we don't know if we reconstructed all showers that "went through" 
  //! our detector.
  for(size_t iVar = 0; iVar < kNVar; ++iVar){
    Spectrum *sAllNuE = LoadFromFile<Spectrum>(inFile, 
      "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=kNoCut_type=NuECC_" \
        + (CRTVeto ? "kCRTHitVetoFD" : "noVetoApplied") ).release();

    TH1D *hAllNuE = sAllNuE->ToTH1( sAllNuE->POT()); 
    float iAllNuE = hAllNuE->Integral();

    //! ////////////////////////////////////////////////////////////////////////
    //!  ^
    //!  |  POT-accounting Spectra/TH1D only created once per variable, but
    //!  \- POT shouldn't change according to variable plotted anyway.
    //! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    //!  /- Spectra/TH1D that are created for every variable-selection 
    //!  |  combination.
    //!  V
    //! ////////////////////////////////////////////////////////////////////////

    //! @attention: we start at index 2 in the selection loop below
    //! since the first index of the individual_cuts_slice vector will
    //! be a singular cut, e.g. kNuEContainedFDCut, and the zero-th index will
    //! always be the kNoCut case. This avoides the following error:
    //!
    //! Error in <TFile::mkdir>: An object with name <name> exists already
    //!
    //! Moreover, we initialize the first kSeqCutName to be the first
    //! index of sels to avoid this duplication logic.
    //!
    //! Additionally, we track the cornerTag and label strings by incrementing
    //! a counter variable so as to not generate long file names and cornerTags.

    unsigned int seqCounter = 1;
    std::string kSeqCutName = sels[ seqCounter].label;

    for(size_t iSel = 2; iSel < kNSel - 2; ++iSel, ++idx){
      if( idx < start || idx >= end) continue; //! for batch processing on grid 

      if constexpr (std::is_same<SelT, SelDefSlice>::value) { //! only do sequential
                                                              //! cuts in slice mode
        kSeqCutName += "_" + sels[iSel].label;
        seqCounter++;
      }
      
      std::string cornerTag = "kSEQ: " \
        + std::to_string( seqCounter) + " many cuts applied" \
        + (CRTVeto ? " CRT Hit Veto Applied" : "");
;
      std::string label = "var=" + vars[iVar].label \
                            + "_unit=" + vars[iVar].unit \
                            + "_sel=SEQ_" + kSeqCutName;
      std::string outLabel = "var=" + vars[iVar].label \
                            + "_unit=" + vars[iVar].unit \
                            + "_sel=SEQ" + std::to_string( seqCounter);

      //! @note: TH1D hMyName is a copy of Spectrum sMyName if such a Spectrum exists.
      TH1D *hNuE,              //! electron neutrinos
           *hTrueSigNuE,       //! 'true signal' electron neutrinos [see above]
           *hNuMu,             //! muon neutrinos
           *hNC,               //! neutral current neutrinos
           *hCosmics,          //! in-time cosmics
           *hTotal,            //! anything with a reconstructed neutrino
           *hOther,            //! anything that's not NC, NuECC, NuMuCC, or in-time cosmic
           *hBackground;       //! anything that's not NuECC

      Spectrum *sNuE        = LoadSpectrum(inFile, label, "NuECC", CRTVeto);
      Spectrum *sTrueSigNuE = LoadSpectrum(inFile, label, "NuECCTrueSignal", CRTVeto);
      Spectrum *sNuMu       = LoadSpectrum(inFile, label, "NuMuCC", CRTVeto);        
      Spectrum *sNC         = LoadSpectrum(inFile, label, "NuNC", CRTVeto);        
      Spectrum *sCosmics    = LoadSpectrum(inFile, label, "NuCosmic", CRTVeto);        
      Spectrum *sTotal      = LoadSpectrum(inFile, label, "NuTotal", CRTVeto);        
          
      hNuE        = sNuE->ToTH1( sNuE->POT());    
      hTrueSigNuE = sTrueSigNuE->ToTH1( sTrueSigNuE->POT());
      hNuMu       = sNuMu->ToTH1( sNuMu->POT());   
      hNC         = sNC->ToTH1( sNC->POT());     
      hCosmics    = sCosmics->ToTH1( sCosmics->POT());
      hTotal      = sTotal->ToTH1( sTotal->POT());  

      hOther = (TH1D*)hTotal->Clone();
      hOther->Add(hNuE, -1);
      hOther->Add(hNuMu, -1);
      hOther->Add(hNC, -1);
      hOther->Add(hCosmics, -1);

      hBackground = (TH1D*)hTotal->Clone();
      hBackground->Add(hNuE, -1);

      //! Extract plotting info (e.g. axes names, titles, etc.) about plot:
      TitleInfo info = parseTitle( hNuE->GetXaxis()->GetTitle());
      std::string xTitle = info.unit.empty() ? info.variable : 
          Form("%s [%s]", info.variable.c_str(), info.unit.c_str());

      //! Make signal to background ratio distribution:
      TH1D *hSigBkgRatio = (TH1D*)hNuE->Clone();
      hSigBkgRatio->Divide( hBackground);

      //! Make efficiency and purity distributions:
      TH1D *hEff = getEffHistogram( hNuE, hBackground, hAllNuE);
      TH1D *hPur = getPurHistogram( hNuE, hBackground);

      //! Format histograms:
      hEff->SetTitle(""); hPur->SetTitle(""); //! strip titles

      formatHist( hNuE, colorNuE, lineNuE, 2);
      formatHist( hTrueSigNuE, colorNuE, lineTrueSigNuE, 2);
      formatHist( hNuMu, colorNuMu, lineNuMu, 2);
      formatHist( hNC, colorNC, lineNC, 2);
      formatHist( hCosmics, colorCosmics, lineCosmics, 2);
      formatHist( hOther, colorOther, lineOther, 3);
      formatHist( hEff, colorEff, kSolid, 2);
      formatHist( hPur, colorPur, kSolid, 2);

      CenterTitles( hNuE);
      CenterTitles( hTrueSigNuE);
      CenterTitles( hNuMu);
      CenterTitles( hNC);
      CenterTitles( hCosmics);
      CenterTitles( hOther);
      CenterTitles( hEff);
      CenterTitles( hPur);

      //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      //! @SummaryTable: calculations only needed for one variable!
      //!                it shouldn't matter if you test iVar == j for any j
      if( !CRTVeto && iVar == 0) { //! can effectively cut on any variable
        float iNuE        = hNuE->Integral();
        float iNuMu       = hNuMu->Integral();
        float iNC         = hNC->Integral();
        float iCosmics    = hCosmics->Integral();
        float iOther      = hOther->Integral();
        float iBackground = hBackground->Integral();
        float iTotal      = hTotal->Integral();

        float iEff = iNuE / iAllNuE;
        float iPur = iNuE / (iNuE + iBackground);
        
        printEventsLine( kSeqCutName, 
          iNuE, iNuMu, iNC, iCosmics, iOther, iEff, iPur, 
          "./outputs_tables/sequential-cut_summaryTable.tex");
      } //! end @SummaryTable
      //!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

      //! Make single-pad TCanvas for plotting different event samples:
      TCanvas *cEvents = new TCanvas( vars[iVar].label.c_str(), 
                                      vars[iVar].label.c_str(),
                                      cXRes, cYRes);
      if( iVar == 0) { //! fake stacking for Number of Slices/Spills variable
        hNuMu->Add( hNuE);
        hNC->Add( hNuMu);
        hCosmics->Add( hNC);
        hOther->Add( hCosmics);

        fillWithDimmedColor( hNuE, 0);
        fillWithDimmedColor( hNuMu, 0);
        fillWithDimmedColor( hNC, 0);
        fillWithDimmedColor( hCosmics, 0);
        fillWithDimmedColor( hOther, 0);
      }

      //! Set y-axis scale:
      if( logScale) {
        gPad->SetLogy();
        hOther->GetYaxis()->SetRangeUser( 1.0, 
                100*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther }));
      }
      else {
        hOther->GetYaxis()->SetRangeUser( 0.0, 
                1.3*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther }));
      }

      //! We work with hOther as the "master" histogram. Since it is (hopefully)
      //! near-zero in every bin, we plot this histogram first and use it to
      //! set the axes titles.
      hOther->GetXaxis()->SetTitle( xTitle.c_str());
      hOther->GetYaxis()->SetTitle( (spillMode ? 
        "Number of Spills" : "Number of Slices") );

      //! Draw hNuE, hNuMu, hNC, hCosmics, and hOther on the same axis.
      hOther->Draw("hist same");
      hCosmics->Draw("hist same");
      hNC->Draw("hist same");
      hNuMu->Draw("hist same");
      hNuE->Draw("hist same");

      gPad->SetGrid(1, 1); //! Set grid on both X and Y axes

      //! Extract the POT to print onto output plots.
      //! @note: POT is an amalgamation of spills. This quantity is not changed
      //! when making selection cuts, separating by interaction type, or
      //! plotting different variables. The POTTag is thus the total amount of
      //! POT from the dataset processed.
      std::string POTTag = makePOTTag( sAllNuE->POT());

      //! Draw plaintext information on events TCanvas.
      DrawSigBkgIntText( hNuE, hBackground, POTTag, 0.04);
      drawComponentsLegend( hNuE, hNuMu, hNC, hCosmics, hOther);
      isSimulation( true);
      myCornerLabel( cornerTag.c_str());

      //! Save the events TCanvas
      cEvents->Update();
      std::string eventsPlotName = outDir + "/" + input + "_" + outLabel \
          + (CRTVeto ? "_CRTHitVetoApplied" : "") \
          + (logScale ? "__logScale" : "") \
          + ".pdf";
      cEvents->Print( eventsPlotName.c_str());
      cEvents->Close();
     
      //! Make split-pad TCanvas for plotting different event samples and
      //! their efficiency and/or purity:
      //! @note: Top pad is variable plotted, bottom is efficiency/purity.
      if( effPur) {
        TCanvas *cEventsEffPur;
        TPad *padEvents, *padEffPur;
        splitCanvasInTwo( cEventsEffPur, padEvents, padEffPur);

        padEvents->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1); 

        //! Set y-axis log scale:
        if( logScale) {
          gPad->SetLogy(1); //! toggle y-axis log scale on
          hOther->GetYaxis()->SetRangeUser( 1.0, 
                100*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther } ));
        }
        else {
          gPad->SetLogy(0); //! toggle y-axis log scale off
          hOther->GetYaxis()->SetRangeUser( 0.0,
                1.3*getHistMaxY( { hNuE, hNuMu, hNC, hCosmics, hOther } ));
        }

        //! ////////////////////////////////////////////////////////////////////
        //! Again we treat hOther as the "master" histogram. Since it's
        //! near-zero in every bin, we plot this histogram first and use it to
        //! set the axes titles.
        //! ////////////////////////////////////////////////////////////////////

        //! @note: Pads share x-axis, so there's no need to put an x-axis
        //!        title on the top pad.
        hOther->GetXaxis()->SetLabelSize(0); //! remove hOther x-axis labels
        hOther->GetXaxis()->SetTitle("");
        
        //! Draw hNuE, hNuMu, hNC, hCosmics, and hOther on the same axis.
        hOther->Draw("hist same");
        hCosmics->Draw("hist same");
        hNC->Draw("hist same");
        hNuMu->Draw("hist same");
        hNuE->Draw("hist same");

        //! Draw plaintext information on events TCanvas.
        DrawSigBkgIntText( hNuE, hBackground, POTTag, 0.04);
        drawComponentsLegend( hNuE, hNuMu, hNC, hCosmics, hOther);
        isSimulation( true);

        myCornerLabel( cornerTag.c_str());
        cEventsEffPur->Update();
        padEffPur->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogy(0); //! never have log scale for efficiency/purity

        //! @note: The TH2 constructor used takes the following args, in order:
        //!        name, title, nXBins, xLow, xUp, nYBins, yLow, yUp
        //! @note: TH1::GetSize() includes one underflow and one overflow bin
        TH2 *axEffPur = new TH2F("", "",
          hNuE->GetSize()-2, 
          hNuE->GetXaxis()->GetXmin(), 
          hNuE->GetXaxis()->GetXmax(),
          30, 0.0, 1.1);
        TH2 *axPur = new TH2F("", "",
          hNuE->GetSize()-2, 
          hNuE->GetXaxis()->GetXmin(), 
          hNuE->GetXaxis()->GetXmax(),
          30, 0.0, 1.1*hSigBkgRatio->GetMaximum());

        //! Set efficiency/purity (bottom) pad's x-axis titles:
        axEffPur->GetXaxis()->SetTitle( xTitle.c_str());
        axPur->GetXaxis()->SetTitle( xTitle.c_str());

        //! Strip y-axis titles for the efficiency/purity pad since it's inferred
        //! that efficiency/purity is a fraction on the interval [0, 1]
        axEffPur->GetYaxis()->SetTitle("");
        axPur->GetYaxis()->SetTitle("");

        if( iSel == 0) { //! no cut? then we don't care about efficiency
          axPur->Draw();
          hPur->Draw("hist");
          DrawPurLegend( hPur, "purity");
        }
        else {
          axEffPur->Draw();
          hEff->Draw("hist same");
          hPur->Draw("hist same");
          DrawIntEffPurLegend( hEff, "efficiency", hPur, "purity");
        }

        gPad->SetGrid(1, 1); //! Set grid on both X and Y axes
        cEventsEffPur->Update();
        gPad->RedrawAxis();

        std::string eventsEffPurPlotName = outDir + "/" + input + "_" + outLabel \
            + (CRTVeto ? "_CRTHitVetoApplied" : "") \
            + (logScale ? "__logScale" : "") \
            + "_effPur.pdf";
        cEventsEffPur->Print( eventsEffPurPlotName.c_str());
        cEventsEffPur->Close();
      } //! end if effPur
    } //! end iVar
  } //! end iSel

  //! Place LaTeX table footer creation outside any loops. For more information
  //! see the @SummaryTable tag above.
  if( !CRTVeto) {
    printTableFooter("./outputs_tables/sequential-cut_summaryTable.tex");
  }
}


/** @fn plotCRTSpectra()
  @brief: Plot CRT Spectra on TCanvas and save the plots as PDFs in outDir.
  @note: This function should only be executed on one input file. This function
         should NOT be called if you are working with multiple input files that
         contain different types of input (e.g. a file for cosmics and a file
         for electron neutrino events).

  @param inFile: Name of the input file
  @param outDir: Path where the output plot PDFs should go
  @param input: Type of input contained within inFile (e.g. cosmics)
  @param vars: Variables you want to plot
  @param sels: Cuts on the variables to plot
  @param rank [optional]: Rank of process (used in batch processing)
  @param nRanks [optional]: Total number of processes (used in batch processing)
  @param logScale [optional]: Set y-axis log scale; default is true
*/
template <typename PlotT, typename SelT>
void plotCRTSpectra(
    const std::string& inFile,
    const std::string& outDir,
    const std::string& input,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const int rank,
    const int nRanks,
    const bool logScale) 
{
  const size_t kNVarCRTSpill = vars.size();
  const size_t kNSelCRTSpill = sels.size();

  int total = kNVarCRTSpill * kNSelCRTSpill;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  //! Make a plot for each variable and cut.
  for(size_t iVar = 0; iVar < kNVarCRTSpill; ++iVar){
    for(size_t iSel = 0; iSel < kNSelCRTSpill; ++iSel, ++idx){
      if( idx < start || idx >= end) continue; //! for batch processing on grid 

      std::string cornerTag = sels[iSel].label;
      std::string label = "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=" + sels[iSel].label \
        + "_type=CRT";

      Spectrum *sCRT = LoadFromFile<Spectrum>( inFile, label).release();
      TH1D *hCRT = sCRT->ToTH1( sCRT->POT());

      formatHist( hCRT, colorCRTHitVetoFDSpillCut, kSolid, 2);
      CenterTitles( hCRT);

      float iCRT = hCRT->Integral();

      TCanvas *cEventsCRT = new TCanvas( vars[iVar].label.c_str(), 
                                         vars[iVar].label.c_str(),
                                         cXRes, cYRes);

      if( logScale) {
        gPad->SetLogy();
        hCRT->GetYaxis()->SetRangeUser( 1.0, 100*hCRT->GetMaximum());
      }
      else {
        hCRT->GetYaxis()->SetRangeUser( 0.0, 1.3*hCRT->GetMaximum());
      }

      TitleInfo info = parseTitle( hCRT->GetTitle());
      std::string xTitle = info.unit.empty() ?
          info.variable : Form("%s [%s]", info.variable.c_str(), info.unit.c_str());

      hCRT->GetXaxis()->SetTitle( xTitle.c_str());
      hCRT->GetYaxis()->SetTitle("Number of Spills");

      hCRT->Draw("hist");

      gPad->SetGrid(1, 1); //! Set grid on both X and Y axes

      myCornerLabel( cornerTag.c_str());
      isSimulation( true);

      //! Save the CRT TCanvas
      cEventsCRT->Update();
      std::string eventsCRTPlotName = outDir + "/" + input + "_" + label \
          + (logScale ? "__logScale" : "") \
          + ".pdf";
      cEventsCRT->Print( eventsCRTPlotName.c_str());
      cEventsCRT->Close();
    } //! end iVar
  } //! end iSel
}