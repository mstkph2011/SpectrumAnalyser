/*
 * THypGeSpectrumAnalyser.cxx
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


#include "THypGeSpectrumAnalyser.h"
#include "THypGePeakFitFunction.h"

#include "TMath.h"
#include "TROOT.h"
#include "TArrow.h"
#include "TText.h"
#include "TNtuple.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"


THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergyspec_ext, Int_t nPeaks_ext, Int_t FuncWidthExt )
{
	fhEnergySpectrum = hEnergyspec_ext;
	nPeaks = nPeaks_ext;
	FuncWidth = FuncWidthExt;
	fSpectrum = new TSpectrum(nPeaks);
	if (CompareNuclei( "eu152")==-1)
		BreakProgram =-1;
	fFittingMethod=1;
		
}

THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,vector<double> *Energies_ext, Int_t FuncWidthExt)			// constructor with energies from external array
{
	fhEnergySpectrum = hEnergySpec_ext;
	nPeaks = Energies.size();
	FuncWidth = FuncWidthExt;
	fSpectrum = new TSpectrum(nPeaks);
	fFittingMethod=1;
}

THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,TString Nuclei, Int_t FuncWidthExt )						// constructor with energies from internal database
{
	fhEnergySpectrum = hEnergySpec_ext;
	FuncWidth = FuncWidthExt;
	nPeaks =CompareNuclei(Nuclei);
	
	fSpectrum = new TSpectrum(nPeaks);
	fFittingMethod=1;
}

// deconstrucor
THypGeSpectrumAnalyser::~THypGeSpectrumAnalyser()
{
	delete fSpectrum;
}

void THypGeSpectrumAnalyser::SetEnergySpectrum(TH1D* hEnergySpec_ext)
{
	fhEnergySpectrum = hEnergySpec_ext;
}

Int_t THypGeSpectrumAnalyser::AnalyseSpectrum()
{
	if (nPeaks == -1 || BreakProgram == -1)
	{
		std::cerr << "*******Nuclei not found!" << endl;
		return -1;
	}
	FindPeaks();
	FitPeaks();
	DoEnergyCalibration();

	CalculateFWHM();
	CalculateFWTM();
	DrawCalibratedSpectrum();		// after FWHM and FWTM for TArrow s
	ExportToTextFile();
	ExportToRootFile();
	return 0;
}

Int_t THypGeSpectrumAnalyser::FindPeaks()
{
	nPeaksFound = fSpectrum->Search(fhEnergySpectrum,5,"",0.001); //must (maybe) be adapted to spectrum; peaks found are sorted by their height (y value)
	PeaksPosX = (Float_t*) fSpectrum->GetPositionX();
	PeaksPosY = (Float_t*) fSpectrum->GetPositionY();
	std::cout << "Found " << nPeaksFound << " Peaks:" <<endl;
	for(Int_t i =0;i< nPeaksFound;i++)
	{
		std::cout << "Peak "<< i+1 <<": X " << PeaksPosX[i]<<"\t\tY " << PeaksPosY[i] << endl;
	}
	PeakSort();				// sorts peaks by their x position to make energy calibration easier
	return 0;
}

Int_t THypGeSpectrumAnalyser::FitPeaks()
{
	char buf[40];
	//TVirtualFitter::SetDefaultFitter("Minuit2");
  //TVirtualFitter::SetDefaultFitter("Fumili2");
  TVirtualFitter::SetPrecision(1e-5);
	THypGePeakFitFunction *func = new THypGePeakFitFunction();
	//FitFunc = new *TF1[nPeaksFound];
	//TF1 *FitFunc[nPeaksFound];
	Int_t npar =0;
	for(Int_t i = 0; i < nPeaksFound;i++)
	{
		
		if (fFittingMethod == 0)		// simple gaussian + step + lin bg
		// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude step, par[4]: gradient of lin bg
		{
			npar = 6;
			sprintf(buf,"FitFunc_%d",i+1);
			FitFunc[i] = new TF1(buf,func,&THypGePeakFitFunction::PeakFuncGaussian,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","PeakFuncGaussian");	// width has to be adapted
		}
		if (fFittingMethod == 1)
		{
			npar = 8;
			sprintf(buf,"FitFunc_%d",i+1);
			FitFunc[i] = new TF1(buf,func,&THypGePeakFitFunction::PeakFunc,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","PeakFunc");	// width has to be adapted
		}
		
		if (fFittingMethod == 2)		// more parameters!!!
		{
			npar = 10;
			sprintf(buf,"FitFunc_%d",i+1);
			FitFunc[i] = new TF1(buf,func,&THypGePeakFitFunction::PeakFuncFreeSkewedPosition,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","PeakFuncFreeSkewedPosition");	// width has to be adapted
		}
		
		if (fFittingMethod == 3)		// second gausian as tail!!!
		{
			npar = 9;
			sprintf(buf,"FitFunc_%d",i+1);
			FitFunc[i] = new TF1(buf,func,&THypGePeakFitFunction::PeakFuncDoubleGausian,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","PeakFuncDoubleGausian");	// width has to be adapted
		}
		
		FitFunc[i]->SetNpx(100000);
		
		//setting parameter names, start values and limits
		sprintf(buf,"AmplGaus_%d",i+1);
		FitFunc[i]->SetParName(0,buf);
		FitFunc[i]->SetParameter(0,PeaksPosY[i]*1.001-fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)));
		
		cout << "ASdasdasdasdasdasd" << fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)) << endl;	
		
		//if (PeaksPosY[i]>1000)			// take small peaks (bigger portion of background) into account (lower lower limit needed)
			FitFunc[i]->SetParLimits(0,PeaksPosY[i]*0.95-fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)),PeaksPosY[i]*1.01-fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)));
		//else
		//	FitFunc[i]->SetParLimits(0,PeaksPosY[i]*0.8,PeaksPosY[i]*1.01);
		
		sprintf(buf,"x0_%d",i+1);
		FitFunc[i]->SetParName(1,buf);
		FitFunc[i]->SetParameter(1,PeaksPosX[i]);
		FitFunc[i]->SetParLimits(1,PeaksPosX[i]*0.98,PeaksPosX[i]*1.02);
		
		sprintf(buf,"sigma_%d",i+1);
		FitFunc[i]->SetParName(2,buf);
		FitFunc[i]->SetParameter(2,1);
		Double_t UpperSigmaFitLimit = FindUpperSigmaFitLimit(PeaksPosY[i]*0.3, PeaksPosX[i]);
		FitFunc[i]->SetParLimits(2,0,UpperSigmaFitLimit);
		//std::cout << FindUpperSigmaFitLimit(PeaksPosY[i]*0.5, PeaksPosX[i])<< endl;
		
		sprintf(buf,"AmplSmoStep_%d",i+1);
		FitFunc[i]->SetParName(3,buf);
		FitFunc[i]->SetParameter(3,10);
		//FitFunc[i]->FixParameter(3,0);
		FitFunc[i]->SetParLimits(3,0.000001,10000);
		
		sprintf(buf,"linBg_%d",i+1);
		FitFunc[i]->SetParName(4,buf);
		FitFunc[i]->SetParameter(4,-0.3);
		FitFunc[i]->SetParLimits(4,-10,0);
		
		sprintf(buf,"constBg_%d",i+1);
		FitFunc[i]->SetParName(5,buf);
		FitFunc[i]->SetParameter(5,100);
		//FitFunc[i]->SetParLimits(5,-1,0);
		if (fFittingMethod == 1 || fFittingMethod == 2)
		{
			// additional parameters for this fitting function
			// par[5]: amplitude skewed gaus, par[6]: beta
			
			sprintf(buf,"AmplSkGaus_%d",i+1);
			FitFunc[i]->SetParName(6,buf);
			FitFunc[i]->SetParameter(6,1000);
			//FitFunc[i]->FixParameter(6,0);
			FitFunc[i]->SetParLimits(6,0.0000001,10000);
			
			sprintf(buf,"beta_%d",i+1);
			FitFunc[i]->SetParName(7,buf);
			FitFunc[i]->SetParameter(7,1);
			//FitFunc[i]->FixParameter(7,1);
			FitFunc[i]->SetParLimits(7,0.000001,1);
			
		}
		if (fFittingMethod == 2)
		{
			// additional parameters for this fitting function
			// par[7]: x0_s, par[8]: sigma_s 
			
			sprintf(buf,"x0SkewedGaus_%d",i+1);
			FitFunc[i]->SetParName(8,buf);
			FitFunc[i]->SetParameter(8,PeaksPosX[i]*0.99);
			//FitFunc[i]->FixParameter(8,0);
			FitFunc[i]->SetParLimits(8,0.98*PeaksPosX[i],PeaksPosX[i]*1.01);
			
			sprintf(buf,"sigmaSkewedGausn_%d",i+1);
			FitFunc[i]->SetParName(9,buf);
			FitFunc[i]->SetParameter(9,1);
			FitFunc[i]->SetParLimits(9,0,UpperSigmaFitLimit);
			//FitFunc[i]->FixParameter(9,1);
			
		}
		
		if (fFittingMethod == 3)
		{
			// additional parameters for this fitting function
			// par[5]: amplitude second gaus, par[6]: x0_ second gaus, par[7]: sigma second gaus
			
			sprintf(buf,"AmplSecondGaus_%d",i+1);
			FitFunc[i]->SetParName(6,buf);
			FitFunc[i]->SetParameter(6,PeaksPosY[i]*0.01);
			//FitFunc[i]->FixParameter(6,0);
			
			FitFunc[i]->SetParLimits(6,PeaksPosY[i]*0.02-fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)),PeaksPosY[i]*0.2-fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(PeaksPosX[i]-FuncWidth)));
			
			sprintf(buf,"x0SecondGaus_%d",i+1);
			FitFunc[i]->SetParName(7,buf);
			FitFunc[i]->SetParameter(7,PeaksPosX[i]-2);
			//FitFunc[i]->FixParameter(7,0);
			FitFunc[i]->SetParLimits(7,PeaksPosX[i]-4,PeaksPosX[i]-1.5);
			
			sprintf(buf,"sigmaSecondGaus_%d",i+1);
			FitFunc[i]->SetParName(8,buf);
			FitFunc[i]->SetParameter(8,1);
			FitFunc[i]->SetParLimits(8,0.1,UpperSigmaFitLimit);
			//FitFunc[i]->FixParameter(8,1);
		}
		std::cout << "Start Fitting Peak " << i << endl;
		if(i==0)
			fhEnergySpectrum->Fit(FitFunc[i],"R");
		else
			fhEnergySpectrum->Fit(FitFunc[i],"R+");
		
		cout << npar << endl;
		//TVirtualFitter *fitter = TVirtualFitter::GetFitter();
		//TMatrixD CovMatrix(npar,npar,fitter->GetCovarianceMatrix());
		//TMatrixD Matrix2(npar,npar,fitter->GetCovarianceMatrix());
		
		
		//CovMatrix.Print();
		//cout << "Determinant: " << CovMatrix.Determinant() << endl;
		
		cout << "Chi²" << FitFunc[i]->GetChisquare() << endl;
		cout << "NDF" << FitFunc[i]-> GetNDF() << endl;
		 cout << "Chi²_red" << FitFunc[i]->GetChisquare()/FitFunc[i]-> GetNDF() << endl;
		PeakFitX.push_back(FitFunc[i]->GetParameter(1));
		PeakFitXError.push_back(FitFunc[i]->GetParError(1));
		
		FitFuncWithoutBg[i] = new TF1(*FitFunc[i]);
		FitFuncWithoutBg[i]->SetParameter(3,0);
		FitFuncWithoutBg[i]->SetParameter(4,0);
		cout << "test" << endl;
		//FitFuncWithoutBg[i]->Draw("same");
		PeakCounts.push_back(FitFuncWithoutBg[i]->Integral(FitFuncWithoutBg[i]->GetXmin(),FitFuncWithoutBg[i]->GetXmax()));
		//std::cout << "Integral:\t" <<FitFuncWithoutBg[i]->Integral(FitFuncWithoutBg[i]->GetXmin(),FitFuncWithoutBg[i]->GetXmax())<< endl;
		
		
		
		
	}
	// draw single components of fit -  new loop is needed to make all of them visible (drawing must be done after all the fiting above, otherwise only components of last peak are drawn (possibly due to "+" option of multiple fits???))
	for(Int_t i = 0; i < nPeaksFound;i++)
	{
		// basic components
		
		sprintf(buf,"FuncGaus_%d",i+1);
		FuncGaus[i] = new TF1(buf,func,&THypGePeakFitFunction::GausOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","GausOnly");
		FuncGaus[i]->SetParameters(FitFunc[i]->GetParameters());
		FuncGaus[i]->SetNpx(100000);
		FuncGaus[i]->SetLineColor( kGreen );
		FuncGaus[i]->Draw( "SAME" );
		
		sprintf(buf,"FuncSmoothedStep_%d",i+1);
		FuncSmoothedStep[i] = new TF1(buf,func,&THypGePeakFitFunction::SmoothedStepOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","SmoothedStepOnly");
		FuncSmoothedStep[i]->SetParameters(FitFunc[i]->GetParameters());
		FuncSmoothedStep[i]->SetNpx(100000);
		FuncSmoothedStep[i]->SetLineColor( kMagenta );
		FuncSmoothedStep[i]->Draw( "SAME" );
		
		sprintf(buf,"FuncLinear_%d",i+1);
		FuncLinear[i] = new TF1(buf,func,&THypGePeakFitFunction::LinearOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","LinearOnly");
		FuncLinear[i]->SetParameters(FitFunc[i]->GetParameters());
		FuncLinear[i]->SetNpx(100000);
		FuncLinear[i]->SetLineColor( kYellow );
		FuncLinear[i]->Draw( "SAME" );
		
		
		//components depending on the fitting model
		
		if (fFittingMethod == 1)
		{
			sprintf(buf,"SkewedGausian_%d",i+1);
			FuncTail[i] = new TF1(buf,func,&THypGePeakFitFunction::SkewedOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","SkewedOnly");
			FuncTail[i]->SetParameters(FitFunc[i]->GetParameters());
			FuncTail[i]->SetNpx(100000);
			FuncTail[i]->SetLineColor( kBlack );
			FuncTail[i]->Draw( "SAME" );
		}
		
		if (fFittingMethod == 2)
		{
			sprintf(buf,"SkewedGausian_%d",i+1);
			FuncTail[i] = new TF1(buf,func,&THypGePeakFitFunction::FreeSkewedOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","FreeSkewedOnly");
			FuncTail[i]->SetParameters(FitFunc[i]->GetParameters());
			FuncTail[i]->SetNpx(100000);
			FuncTail[i]->SetLineColor( kBlack );
			FuncTail[i]->Draw( "SAME" );
		}
		
		if (fFittingMethod == 3)
		{
			sprintf(buf,"SecondGausian_%d",i+1);
			FuncTail[i] = new TF1(buf,func,&THypGePeakFitFunction::SecondGausOnly,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,npar,"THypGePeakFitFunction","SecondGausOnly");
			FuncTail[i]->SetParameters(FitFunc[i]->GetParameters());
			FuncTail[i]->SetNpx(100000);
			FuncTail[i]->SetLineColor( kBlack );
			FuncTail[i]->Draw( "SAME" );
		}
	}
	
	
	fhSubstractSpectrum = new TH1D("fhSubstractSpectrum","Substracted spectrum; Channel; Counts", fhEnergySpectrum->GetNbinsX(),fhEnergySpectrum->GetXaxis()->GetXmin(),fhEnergySpectrum->GetXaxis()->GetXmax());
	fSubstractCanvas= new TCanvas("fSubstractCanvas","Substracted spectrum",800,600);
	fSubstractCanvas->cd();
	for (Int_t i=0;i < nPeaksFound; i++)
	{
		for (Double_t iValue = (Double_t) PeaksPosX[i]-FuncWidth; iValue <= (Double_t) PeaksPosX[i]+FuncWidth; iValue += fhEnergySpectrum->GetBinWidth(100))
		{
			fhSubstractSpectrum->Fill(fhEnergySpectrum->GetBinCenter(fhEnergySpectrum->FindBin(iValue)), fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(iValue))-FitFunc[i]->Eval(fhEnergySpectrum->GetBinCenter(fhEnergySpectrum->FindBin(iValue))));
			fhSubstractSpectrum->SetBinError(fhEnergySpectrum->FindBin(iValue), sqrt(fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(iValue))));
			//cout <<fhEnergySpectrum->GetBinContent(fhEnergySpectrum->FindBin(iValue))-FitFunc[i]->Eval(fhEnergySpectrum->GetBinCenter(fhEnergySpectrum->FindBin(iValue))) <<" " <<FitFunc[i]->Eval(fhEnergySpectrum->GetBinCenter(fhEnergySpectrum->FindBin(iValue))) <<endl;
		}
	}
	fhSubstractSpectrum->Draw("E1");
	return 0;
}
Int_t THypGeSpectrumAnalyser::DoEnergyCalibration()
{
	std::cout << "Starting energy calibration" << endl;
	fCalCanvas = new TCanvas("fCalCanvas","Energy Calibration",800,600);				//all "canvas stuff" may need changes for go4
	fCalCanvas->cd();

	fgCalibration = new TGraphErrors(nPeaksFound,&(PeakFitX[0]),&(Energies[0]),&(PeakFitXError[0]),&(EnergyErrors[0]));
	CalFunc = new TF1("CalFunc","pol2",0,8000);
	fgCalibration->Draw("AP");
	fgCalibration->Fit(CalFunc);
	
	return 0;
}

Int_t THypGeSpectrumAnalyser::CalculateFWHM()
{
	for(Int_t i = 0; i <nPeaksFound; i++)
	{
		Amplitude[i] = FitFunc[i]->GetMaximum()-FitFunc[i]->GetMinimum();				// Easy way to substract the background, only the lower part on the high energy side is substracted
		FWHMHeight[i] = FitFunc[i]->GetMaximum()-Amplitude[i]/2;
		FWHMlow[i] = FitFunc[i]->GetX(FWHMHeight[i],FitFunc[i]->GetXmin(),FitFunc[i]->GetMaximumX());
		FWHMhigh[i] = FitFunc[i]->GetX(FWHMHeight[i],FitFunc[i]->GetMaximumX(),FitFunc[i]->GetXmax());
		FWHMlowEnergy[i] = Calibrate(FWHMlow[i]);
		FWHMhighEnergy[i] = Calibrate(FWHMhigh[i]);
		FWHM_ch.insert(std::pair<double,double>(Energies[i],FWHMhigh[i]-FWHMlow[i]));
		FWHM_en.insert(std::pair<double,double>(Energies[i],FWHMhighEnergy[i]-FWHMlowEnergy[i]));
		//std::cout << "FWHMlow" << FWHMlow << "\t\tFWHMhigh" << FWHMhigh << endl;
		//std::cout << "FWHMlowEn" << FWHMlowEnergy << "\t\tFWHMhighEn" << FWHMhighEnergy << endl;
		//std::cout << FWHMhighEnergy-FWHMlowEnergy<< endl;
	}
	//std::map<double,double>::const_iterator it=FWHM_en.begin();
	//for (it=FWHM_en.begin(); it!=FWHM_en.end(); ++it)
		//std::cout << it->first << " => " << it->second << '\n';
	return 0;
}

Int_t THypGeSpectrumAnalyser::CalculateFWTM()
{
	for(Int_t i = 0; i <nPeaksFound; i++)
	{
		Amplitude[i] = FitFunc[i]->GetMaximum()-FitFunc[i]->GetMinimum();				// Easy way to substract the background, only the lower part on the high energy side is substracted
		FWTMHeight[i] = FitFunc[i]->GetMaximum()-Amplitude[i]*9/10;
		FWTMlow[i] = FitFunc[i]->GetX(FWTMHeight[i],FitFunc[i]->GetXmin(),FitFunc[i]->GetMaximumX());
		FWTMhigh[i] = FitFunc[i]->GetX(FWTMHeight[i],FitFunc[i]->GetMaximumX(),FitFunc[i]->GetXmax());
		FWTMlowEnergy[i] = Calibrate(FWTMlow[i]);
		FWTMhighEnergy[i] = Calibrate(FWTMhigh[i]);
		FWTM_ch.insert(std::pair<double,double>(Energies[i],FWTMhigh[i]-FWTMlow[i]));
		FWTM_en.insert(std::pair<double,double>(Energies[i],FWTMhighEnergy[i]-FWTMlowEnergy[i]));
		//std::cout << "FWTMlow" << FWTMlow << "\t\tFWTMhigh" << FWTMhigh << endl;
		//std::cout << "FWTMlowEn" << FWTMlowEnergy << "\t\tFWTMhighEn" << FWTMhighEnergy << endl;
		//std::cout << FWTMhighEnergy-FWTMlowEnergy<< endl;
	}
	//std::map<double,double>::const_iterator it=FWTM_en.begin();
	//for (it=FWTM_en.begin(); it!=FWTM_en.end(); ++it)
		//std::cout << it->first << " => " << it->second << '\n';
	
	return 0;
}


Int_t THypGeSpectrumAnalyser::CompareNuclei(TString NucleiName)
{
	NucleiName.ToLower();
	if(!NucleiName.CompareTo("co60") )
	{
		//npeaks = 2;
		std::cout << "Found Co60"<< endl;
		Energies.push_back(1172.23);
		EnergyErrors.push_back(0.030);
		Energies.push_back(1332.51);
		EnergyErrors.push_back(0.015);
		return Energies.size();
	}
	
		if(!NucleiName.CompareTo("co60bg") )
	{
		//npeaks = 2;
		std::cout << "Found Co60 and K40"<< endl;			
		Energies.push_back(1172.23);				//keV
		EnergyErrors.push_back(0.030);			//keV
		Energies.push_back(1332.51);				//keV
		EnergyErrors.push_back(0.015);			//keV
		Energies.push_back(1460.83);				//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		return Energies.size();					// 3 peaks
	}
	
	if( !NucleiName.CompareTo( "eu152") )
	{
		std::cout << "Found Eu152"<< endl;
		Energies.push_back(121.78);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(244.70);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(344.28); 				//keV	 
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(411.1);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(444.0); 					//keV	 
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(778.90);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(867.4);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(964.08); 				//keV	 
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(1085.9);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(1112.1); 				//keV	 
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		Energies.push_back(1408.0);					//keV
		EnergyErrors.push_back(0.1);				// placeholder, real value unknown
		return Energies.size();				// 11 peaks
	}
	if( !NucleiName.CompareTo( "eu152bg") )
	{
		std::cout << "Found Eu152 and K40"<< endl;
		return Energies.size();
	}
	
	if(!NucleiName.CompareTo("cs137") )
	{
		std::cout << "Found Cs137"<< endl;
		return Energies.size();
	}
	
	if(!NucleiName.CompareTo("cs137bg") )
	{
		std::cout << "Found Cs137 and K40"<< endl;
		return Energies.size();
	}
	
	if(!NucleiName.CompareTo("na22") )
	{
		std::cout << "Found Na22"<< endl;
		return Energies.size();
	}
	
	if(!NucleiName.CompareTo("na22bg") )
	{
		std::cout << "Found Na22 and K40"<< endl;
		return Energies.size();
	}
	//std::cerr << "****************Nuclei not found!***************" <<endl;
	return -1;
}

Double_t THypGeSpectrumAnalyser::FindUpperSigmaFitLimit(Double_t threshold, Float_t StartingPoint) 
{
	//find fitting range of sigma. this is done by starting at "starting point" and than finding the position when the peak goes bellow a certain threshold.
	Int_t StartBin = fhEnergySpectrum->FindBin((Double_t)StartingPoint);
   //std::cout << "StartBin" << StartBin;
   
   for (Int_t bin=StartBin;bin<=fhEnergySpectrum->GetNbinsX();bin++) 
   {
		 //std::cout << bin<<"\t\t\tSpecValue\t\t"<<fhEnergySpectrum->GetBinContent(bin) << endl;
      if (fhEnergySpectrum->GetBinContent(bin) < threshold)
      {
				//std::cout <<  fhEnergySpectrum->GetBinCenter(bin)<< "\t\t"<<StartingPoint << endl;
				return fhEnergySpectrum->GetBinCenter(bin)-StartingPoint;
			}
   }
   return -1;
}

void THypGeSpectrumAnalyser::PeakSort()
{
	std::cout << "Sorting Peaks" << endl;
	Float_t *arrayOrig = new Float_t[nPeaksFound];
	for(int i = 0;i<nPeaksFound; ++i)
		arrayOrig[i] = PeaksPosX[i];
  std::sort(PeaksPosX, PeaksPosX + nPeaksFound); //sorting Peak Positions (lower to greater x)
        
	Float_t buffy[nPeaksFound];
   for (Int_t i = 0;i < nPeaksFound; ++i)		//placing Amplitude to the coresponding Peak Position
   {
			for (Int_t j = 0; j < nPeaksFound; ++j)
			{
				if (arrayOrig[i]==PeaksPosX[j])
				{
					buffy[j] = PeaksPosY[i];
				}
			}
		}
		for(int i = 0;i<nPeaksFound; ++i)
		{
			PeaksPosY[i] = buffy[i];
		}
		for (int i = 0; i < nPeaksFound; ++i) 
     std::cout <<"("<< PeaksPosX[i]<<"," << PeaksPosY[i] << ") ";
    std::cout << endl;	
}


Int_t THypGeSpectrumAnalyser::DrawCalibratedSpectrum()
{
	std::cout << "Drawing of CalibratedSpectrum" <<endl;
	fCalSpecCanvas = new TCanvas("fCalSpecCanvas", "Calibrated spectrum",800,600);				//all "canvas stuff" may need changes for go4
	fCalSpecCanvas->cd();
	fhEnergySpectrum->GetXaxis()->UnZoom();
	//fhCalibratedSpectrum = new TH1D("fhCalibratedSpectrum","Calibrated spectrum;Energy [keV];Counts [a.u.]",fhEnergySpectrum->GetNbinsX(),Calibrate(fhEnergySpectrum->GetXaxis()->GetBinCenter(fhEnergySpectrum->GetXaxis()->GetFirst())),Calibrate(fhEnergySpectrum->GetXaxis()->GetBinCenter(fhEnergySpectrum->GetXaxis()->GetLast())));
	
	//create new bins for calibrated spectrum. These bins have variable bin size due to the non-linear energy calibration
	const int nbins = fhEnergySpectrum->GetNbinsX(); 
	double new_bins[nbins+1]; 
	for(int i=0; i <= nbins; i++)
	{ 
		new_bins[i] = Calibrate( fhEnergySpectrum-> GetBinLowEdge(i+1) ); 
		//std::cout << new_bins[i]<<endl;
	} 
	fhCalibratedSpectrum = new TH1D("fhCalibratedSpectrum","Calibrated spectrum;Energy [keV];Counts [a.u.]",fhEnergySpectrum->GetNbinsX(),new_bins);
	for(Int_t i = 0; i<=fhEnergySpectrum->GetNbinsX(); i++)
	{
		fhCalibratedSpectrum->Fill(Calibrate(fhEnergySpectrum->GetXaxis()->GetBinCenter(i)),fhEnergySpectrum->GetBinContent(i));
	}
	fhCalibratedSpectrum->SetStats(0);
	fhCalibratedSpectrum->Draw();
	
	std::map<double,double>::const_iterator it=FWHM_en.begin();
	std::map<double,double>::const_iterator it2=FWTM_en.begin();
	for(Int_t i = 0;i< nPeaksFound; i++)		
	{
		
		CalibratedFunction[i] = new TF1(*FitFunc[i]);
		//adapt parameters to energy spectrum
		Double_t OldSigma = CalibratedFunction[i]->GetParameter(2);
		CalibratedFunction[i]->SetParameter(2,Calibrate(CalibratedFunction[i]->GetParameter(1))-Calibrate(CalibratedFunction[i]->GetParameter(1)-CalibratedFunction[i]->GetParameter(2)));
		CalibratedFunction[i]->SetParameter(4,Calibrate(CalibratedFunction[i]->GetParameter(1))-Calibrate(CalibratedFunction[i]->GetParameter(1)-CalibratedFunction[i]->GetParameter(4)));
		CalibratedFunction[i]->SetParameter(1,Calibrate(CalibratedFunction[i]->GetParameter(1)));
		//adapt amplitudes, this has to be done, because of the par[2] (sigma) in the denominator in front of the gausian part. Since the gausian ampl. (par[0]) is in the denominator of any other part their ampl. must be adapted too. The lin. bg gradient must not be adapted, because the factors from the change of x and par[0] cancel each other
		//CalibratedFunction[i]->SetParameter(0,CalibratedFunction[i]->GetParameter(0)*CalibratedFunction[i]->GetParameter(2)/OldSigma);
		//std::cout << "new ampl" << CalibratedFunction[i]->GetParameter(0) << endl;
		//CalibratedFunction[i]->SetParameter(3,CalibratedFunction[i]->GetParameter(3)*CalibratedFunction[i]->GetParameter(2)/OldSigma);
		//CalibratedFunction[i]->SetParameter(5,CalibratedFunction[i]->GetParameter(5)*CalibratedFunction[i]->GetParameter(2)/OldSigma);
		//function range must be adapted
		CalibratedFunction[i]->SetRange(Calibrate(CalibratedFunction[i]->GetXmin()),Calibrate(CalibratedFunction[i]->GetXmax()));
		
		CalibratedFunction[i]->Draw("same");
		TArrow *FWHMArrow = new TArrow(FWHMlowEnergy[i], FWHMHeight[i],FWHMhighEnergy[i],FWHMHeight[i],0.02,"<|>");
		FWHMArrow->Draw("");
		TArrow *FWTMArrow = new TArrow(FWTMlowEnergy[i], FWTMHeight[i],FWTMhighEnergy[i],FWTMHeight[i],0.02,"<|>");
		FWTMArrow->Draw("");
		char TextBuffer[50];
		sprintf(TextBuffer,"%.2f keV",it->second);
		//TextBuffer = "test";
		TText *FWHMText = new TText(CalibratedFunction[i]->GetParameter(1),FWHMHeight[i]*1.05, TextBuffer);
		FWHMText->SetTextAlign(22);
		FWHMText->SetTextFont(63);
		FWHMText->SetTextSizePixels(22);
		FWHMText->Draw();
		sprintf(TextBuffer,"%.2f keV",it2->second);
		TText *FWTMText = new TText(CalibratedFunction[i]->GetParameter(1),FWTMHeight[i]+FWHMHeight[i]*0.05, TextBuffer);
		FWTMText->SetTextAlign(22);
		FWTMText->SetTextFont(63);
		FWTMText->SetTextSizePixels(22);
		FWTMText->Draw();
		it++;
		it2++;
	}
	return 0;
}

Double_t THypGeSpectrumAnalyser::Calibrate(Double_t Channel)
{
	Double_t a = CalFunc->GetParameter(0);
	Double_t b = CalFunc->GetParameter(1);
	Double_t c = CalFunc->GetParameter(2);
	return a + b *Channel + c * Channel * Channel;
}


Int_t THypGeSpectrumAnalyser::ExportToTextFile(TString TxtFilename_ext )
{
	TString TxtFilenameWithPath;
	if (OutputPath)
		TxtFilenameWithPath = OutputPath;
	if (!TxtFilename)
		TxtFilename = TxtFilename_ext;
	TxtFilenameWithPath += TxtFilename;
	TxtFile.open(TxtFilenameWithPath.Data());
	TxtFile << "File:\t" << fhEnergySpectrum->GetTitle()<< endl;
	TxtFile << "Nuclei:\t" << endl;
	TxtFile << "Energy [keV]\t\tFWHM [keV]\t\tFWTM [keV]\t\tFWTM/FWHM\t\tPeakCounts" << endl;
	
	std::map<double,double>::const_iterator it=FWHM_en.begin();
	std::map<double,double>::const_iterator it2=FWTM_en.begin();
	Int_t j=0;
	for (it=FWHM_en.begin(); it!=FWHM_en.end(); ++it)
	{
		TxtFile << it->first << "\t\t" << it->second <<"\t\t"<< it2->second <<"\t\t" << it2->second/it->second << "\t\t" << PeakCounts[j] << endl;
		it2++;
		j++;
	}
	TxtFile.close();
	return 0;
}
Int_t THypGeSpectrumAnalyser::ExportToRootFile(TString RootFilename_ext )
{
	TString RootFilenameWithPath;
	if (OutputPath)
		RootFilenameWithPath = OutputPath;
	if (!RootFilename)
		RootFilename += RootFilename_ext;
	RootFilenameWithPath += RootFilename;
	
	RootFile = (TFile*)gROOT->FindObject(RootFilenameWithPath); 
	if (RootFile) 
		RootFile->Close();
	RootFile = new TFile(RootFilenameWithPath,"RECREATE",RootFilenameWithPath);
		
	TNtuple* DataNTuple = new TNtuple("DataNTuple","Data ntuple","PeakNumber:Energy:FWHM:FWTM:FWTMtoFWHMratio:PeakCounts");
	Double_t wEnergy,wFWHM,wFWTM,wRatio,wPeakCounts;				// w means "write"
	Int_t wPeakNumber;
	std::map<double,double>::const_iterator it=FWHM_en.begin();
	std::map<double,double>::const_iterator it2=FWTM_en.begin();
	for(Int_t i = 0; i < nPeaksFound; i++)
	{
		wPeakNumber = i;
		wEnergy= it->first;
		wFWHM = it->second;
		wFWTM = it2->second;
		wRatio = it2->second/it->second;
		wPeakCounts = PeakCounts[i];
		DataNTuple->Fill(wPeakNumber,wEnergy,wFWHM,wFWTM,wRatio,wPeakCounts);
		it++;
		it2++;
	}
	DataNTuple->Write();
	fCalCanvas->Write();
	fCalSpecCanvas->Write();
	RootFile->Close();
	return 0;
}

void THypGeSpectrumAnalyser::SetSearchRange(Double_t RangeMin,Double_t RangeMax)
{
	fhEnergySpectrum->GetXaxis()->SetRange(RangeMin,RangeMax);
}

void THypGeSpectrumAnalyser::SetTxtFileOutputName(TString TxtFilename_ext)
{
	TxtFilename = TxtFilename_ext; 
}
void THypGeSpectrumAnalyser::SetRootFileOutputName(TString RootFilename_ext)
{
	 RootFilename = RootFilename_ext;
}

void THypGeSpectrumAnalyser::SetOutputPath(TString OutputPath_ext)
{
	OutputPath = OutputPath_ext;
}


void THypGeSpectrumAnalyser::SetGaussianFitting()
{
	fFittingMethod =0;
}

Bool_t THypGeSpectrumAnalyser::IsGaussianFitting()
{
	if (fFittingMethod == 0)
		return 1;
	else
		return 0;
}


void THypGeSpectrumAnalyser::SetFreeSkewedFitting()
{
	fFittingMethod =2;
}

Bool_t THypGeSpectrumAnalyser::IsFreeSkewedFitting()
{
	if (fFittingMethod == 2)
		return 1;
	else
		return 0;
}

void THypGeSpectrumAnalyser::SetSecondGausianFitting()
{
	fFittingMethod =3;
}

Bool_t THypGeSpectrumAnalyser::IsSecondGausianFitting()
{
	if (fFittingMethod == 3)
		return 1;
	else
		return 0;
}
