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


THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergyspec_ext, Int_t nPeaks_ext, Int_t FuncWidthExt )
{
	fhEnergySpectrum = hEnergyspec_ext;
	nPeaks = nPeaks_ext;
	FuncWidth = FuncWidthExt;
	fSpectrum = new TSpectrum(nPeaks);
	if (CompareNuclei( "eu152")==-1)
		BreakProgram =-1;
		
}

THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,vector<double> *Energies_ext, Int_t FuncWidthExt)			// constructor with energies from external array
{
	fhEnergySpectrum = hEnergySpec_ext;
	nPeaks = Energies.size();
	FuncWidth = FuncWidthExt;
	fSpectrum = new TSpectrum(nPeaks);
}

THypGeSpectrumAnalyser::THypGeSpectrumAnalyser(TH1D* hEnergySpec_ext,TString Nuclei, Int_t FuncWidthExt )						// constructor with energies from internal database
{
	fhEnergySpectrum = hEnergySpec_ext;
	FuncWidth = FuncWidthExt;
	nPeaks =CompareNuclei(Nuclei);
	
	fSpectrum = new TSpectrum(nPeaks);
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
	PeaksPosX = fSpectrum->GetPositionX();
	PeaksPosY = fSpectrum->GetPositionY();
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
	char buf[20];
	//TVirtualFitter::SetDefaultFitter("Minuit2");
  //TVirtualFitter::SetDefaultFitter("Fumili2");
  TVirtualFitter::SetPrecision(1e-5);
	THypGePeakFitFunction *func = new THypGePeakFitFunction();
	//FitFunc = new *TF1[nPeaksFound];
	//TF1 *FitFunc[nPeaksFound];
	for(Int_t i = 0; i < nPeaksFound;i++)
	{
		sprintf(buf,"FitFunc_%d",i+1);
		FitFunc[i] = new TF1(buf,func,&THypGePeakFitFunction::PeakFunc,(Double_t) PeaksPosX[i]-FuncWidth,(Double_t) PeaksPosX[i]+FuncWidth,7,"THypGePeakFitFunction","PeakFunc");	// width has to be adapted
		FitFunc[i]->SetNpx(100000);
		// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude skewed gaus, par[4]: beta,par[5]: amplitude smooted step, par[6]: constant
		sprintf(buf,"AmplGaus_%d",i+1);
		FitFunc[i]->SetParName(0,buf);
		FitFunc[i]->SetParameter(0,PeaksPosY[i]);
		FitFunc[i]->SetParLimits(0,PeaksPosY[i]*0.1,PeaksPosY[i]*10);
		
		sprintf(buf,"x0_%d",i+1);
		FitFunc[i]->SetParName(1,buf);
		FitFunc[i]->SetParameter(1,PeaksPosX[i]);
		FitFunc[i]->SetParLimits(1,PeaksPosX[i]*0.9,PeaksPosX[i]*1.1);
		
		sprintf(buf,"sigma_%d",i+1);
		FitFunc[i]->SetParName(2,buf);
		FitFunc[i]->SetParameter(2,1);
		Double_t UpperSigmaFitLimit = FindUpperSigmaFitLimit(PeaksPosY[i]*0.3, PeaksPosX[i]);
		FitFunc[i]->SetParLimits(2,0,UpperSigmaFitLimit);
		//std::cout << FindUpperSigmaFitLimit(PeaksPosY[i]*0.5, PeaksPosX[i])<< endl;
		
		sprintf(buf,"AmplSkGaus_%d",i+1);
		FitFunc[i]->SetParName(3,buf);
		FitFunc[i]->SetParameter(3,0e-5);
		//FitFunc[i]->FixParameter(3,0);
		FitFunc[i]->SetParLimits(3,0,100);
		
		sprintf(buf,"beta_%d",i+1);
		FitFunc[i]->SetParName(4,buf);
		FitFunc[i]->SetParameter(4,1);
		//FitFunc[i]->FixParameter(4,1);
		FitFunc[i]->SetParLimits(4,0,10);
		
		sprintf(buf,"AmplSmoStep_%d",i+1);
		FitFunc[i]->SetParName(5,buf);
		FitFunc[i]->SetParameter(5,10);
		//FitFunc[i]->FixParameter(5,0);
		//FitFunc[i]->SetParLimits(5,0,PeaksPosY[i]*0.1*1e6);
		
		sprintf(buf,"cont_%d",i+1);
		FitFunc[i]->SetParName(6,buf);
		//FitFunc[i]->SetParameter(6,100);
		//FitFunc[i]->SetParLimits(6,0,PeaksPosY[i]*0.1*1e6);
		std::cout << "Start Fitting Peak " << i << endl;
		if(i==0)
			fhEnergySpectrum->Fit(FitFunc[i],"R");
		else
			fhEnergySpectrum->Fit(FitFunc[i],"R+");
		PeakFitX.push_back(FitFunc[i]->GetParameter(1));
		PeakFitXError.push_back(FitFunc[i]->GetParError(1));
		
		FitFuncWithoutBg[i] = new TF1(*FitFunc[i]);
		FitFuncWithoutBg[i]->SetParameter(5,0);
		FitFuncWithoutBg[i]->SetParameter(6,0);
		//FitFuncWithoutBg[i]->Draw("same");
		PeakCounts.push_back(FitFuncWithoutBg[i]->Integral(FitFuncWithoutBg[i]->GetXmin(),FitFuncWithoutBg[i]->GetXmax()));
		//std::cout << "Integral:\t" <<FitFuncWithoutBg[i]->Integral(FitFuncWithoutBg[i]->GetXmin(),FitFuncWithoutBg[i]->GetXmax())<< endl;
	}
	
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

//Double_t THypGeSpectrumAnalyser::PeakFunc (Double_t *x,Double_t *par)
//{
			//// Peak shape is taken from gf2 http://www.phy.anl.gov/gammasphere/doc/gf2.hlp
			//// composed by the sum of:
			////	(1) a Gaussian,
			////	(2) a skewed Gaussian,
			////	(3) a smoothed step function to increase the background on the low-energy
			////			side of the peak. Components (2) and/or (3) can easily be set to zero if not
			////			required, and
			////	(4)	a constant				// not in gf2
			////
			////	(2)	y = constant
			////			* EXP( (x-c)/beta )
			////			* ERFC( (x-c)/(SQRT(2)*sigma) + sigma/(SQRT(2)*beta) )
			////	(3)	y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
			
			//// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude skewed gaus, par[4]: beta,par[5]: amplitude smooted step, par[6]: constant

			//Double_t Gausian = par[0]/TMath::Sqrt(2 * TMath::Pi())/par[2] * TMath::Exp(- 0.5*pow((x[0] - par[1])/par[2],2)); // par[0]: amplitude, par[1]: x[0]0; par[2]: sigma
			//Double_t SkewedGausian = par[3] 
							//* TMath::Exp((x[0]-par[1])/par[4]) 
							//* TMath::Erfc( ((x[0] - par[1]) /(TMath::Sqrt(2)*par[2]) ) + par[2]/TMath::Sqrt(2)/par[4]); //par[3]: amplitude, par[4]: beta
			//Double_t SmoothedStep = par[5] * TMath::Erfc( (x[0] - par[1]) /(TMath::Sqrt(2)*par[2])); // par[5]: amplitude
			//Double_t Constant = par[6]; // par[6]: constant
	//return Gausian + SkewedGausian + SmoothedStep +  Constant;
//}


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
