/*
 * THypGeSpectrumAnalyser.h
 * 
 * Copyright 2013 Marcell Steinen <steinen@kph.uni-mainz.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

// "IsA()" Error is given by ClassDef line. This line must be used for Go4 and PANDAroot. 

#ifndef THYPGESPECTRUMANALYSER_H
#define THYPGESPECTRUMANALYSER_H

#include "TH1D.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

using namespace std;

class THypGeSpectrumAnalyser
{
	public:
		THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,Int_t nPeaks =2,Int_t FuncWidthExt = 100);
		THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext, vector<double> *Energies_ext,Int_t FuncWidthExt = 100);			// constructor with energies from external array
		THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,TString Nuclei,Int_t FuncWidthExt = 100);						// constructor with energies from internal database
		virtual ~THypGeSpectrumAnalyser();
	
		TH1D*											GetEnergySpectrum() {return fhEnergySpectrum;}
		void 											SetEnergySpectrum(TH1D* hEnergySpec_ext);
	
		Int_t 										AnalyseSpectrum();
		Int_t											FindPeaks();
		Int_t											FitPeaks();

		Int_t											DrawCalibratedSpectrum();
		Int_t											DoEnergyCalibration();
		Int_t											CalculateFWHM();
		Int_t											CalculateFWTM();
		
		Int_t											CompareNuclei(TString NucleiName);						// database comparison
		
		void											SetTxtFileOutputName(TString TxtFilename_ext);
		void											SetRootFileOutputName(TString RootFilename_ext);
		void											SetOutputPath(TString OutputPath_ext);
		
		Int_t 										ExportToTextFile(TString TxtFilename_ext = "test.txt");
		Int_t											ExportToRootFile(TString RootFilename_ext = "test.root");
		
		void											SetSearchRange(Double_t RangeMin,Double_t RangeMax);
		
		

	private:

		Double_t 									FindUpperSigmaFitLimit(Double_t threshold, Float_t StartingPoint);
		void 											PeakSort();
		Double_t 									Calibrate(Double_t Channel);
		
		/* add your private declarations */
		TH1D*											fhEnergySpectrum;
		
		TCanvas										*fCalSpecCanvas;
		TH1D*											fhCalibratedSpectrum;
		TSpectrum*								fSpectrum;
		Float_t										*PeaksPosX;
		Float_t										*PeaksPosY;
		
		Int_t 										nPeaks;
		Int_t 										nPeaksFound;		//number of Peaks found by fSpectrum->Search("...")

		TF1 											*FitFunc[50];
		TF1 											*FitFuncWithoutBg[50];
		TF1												*CalibratedFunction[50];

		Int_t 										FuncWidth;
		vector<double> 					PeakFitX;
		vector<double> 					PeakFitXError;
		vector<double>						PeakCounts;
		
		TGraphErrors							*fgCalibration;
		TCanvas										*fCalCanvas;
		TF1 											*CalFunc;
		
		Double_t 									Amplitude[50];
		Double_t									FWHMHeight[50];
		Double_t 									FWHMlow[50];
		Double_t 									FWHMhigh[50];
		Double_t 									FWHMlowEnergy[50];
		Double_t 									FWHMhighEnergy[50];
		
		Double_t									FWTMHeight[50];
		Double_t 									FWTMlow[50];
		Double_t 									FWTMhigh[50];
		Double_t 									FWTMlowEnergy[50];
		Double_t									FWTMhighEnergy[50];
		

		map<double,double>			FWHM_ch;
		map<double,double>			FWHM_en;
		map<double,double>			FWTM_ch;
		map<double,double>			FWTM_en;
		
		vector<double>						Energies;				//vector of peak energies
		vector<double>						EnergyErrors;	//vector of peak energy errors
		
		TString 									OutputPath;
		
		ofstream 									TxtFile;
		TString										TxtFilename;
		
		TFile 										*RootFile;
		TString										RootFilename;
		
		Bool_t										BreakProgram;
		
		//ClassDef(THypGeSpectrumAnalyser,1)					// This line must be used for Go4 and PANDAroot !!!!!!!!!!!   gives "IsA()" Error if compiled with g++
};

#endif /* THYPGESPECTRUMANALYSER_H */ 
