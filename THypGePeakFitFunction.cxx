/*
 * THypGePeakFitFunction.cxx
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

	//support class for peak fitting, ROOT needs this external construction (?)
#include "THypGePeakFitFunction.h"


THypGePeakFitFunction::THypGePeakFitFunction()
{
	
}

Double_t THypGePeakFitFunction::PeakFuncGaussian (Double_t *x,Double_t *par)
{
			// Peak shape is taken from gf2 http://www.phy.anl.gov/gammasphere/doc/gf2.hlp
			// composed by the sum of:
			//	(1) a Gaussian,
			//	(2) a smoothed step function to increase the background on the low-energy
			//			side of the peak. Components (2) and/or (3) can easily be set to zero if not
			//			required, and
			//	(3)	a constant				// not in gf2
			//

			//	(2)	y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
			
			// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude step, par[4]: gradient of lin bg
			
			// with sigma in amplitude
			//Double_t Gausian = par[0]*TMath::Sqrt(2 * TMath::Pi())/par[2] * TMath::Exp(- 0.5*pow((x[0] - par[1])/par[2],2)); // par[0]: amplitude, par[1]: x[0]0; par[2]: sigma
			// or without sigma in amplitude
			
			//Double_t Gausian = par[0]/*TMath::Sqrt(2 * TMath::Pi())/par[2]*/ * TMath::Exp(- 0.5*pow((x[0] - par[1])/par[2],2)); // par[0]: amplitude, par[1]: x[0]0; par[2]: sigma
			//Double_t SmoothedStep = par[3] * TMath::Erfc( (x[0] - par[1]) /(TMath::Sqrt(2)*par[2])); // par[3]: amplitude
			//Double_t Constant = par[4]*x[0]; // par[4]: gradient of lin b
	//return Gausian + SmoothedStep +  Constant;
	return GausOnly(x,par) + SmoothedStepOnly(x,par) + LinearOnly(x,par);
}

Double_t THypGePeakFitFunction::PeakFunc (Double_t *x,Double_t *par)
{
			// Peak shape is taken from gf2 http://www.phy.anl.gov/gammasphere/doc/gf2.hlp
			// composed by the sum of:
			//	(1) a Gaussian,
			//	(2) a skewed Gaussian,
			//	(3) a smoothed step function to increase the background on the low-energy
			//			side of the peak. Components (2) and/or (3) can easily be set to zero if not
			//			required, and
			//	(4)	a constant				// not in gf2
			//
			//	(2)	y = constant
			//			* EXP( (x-c)/beta )
			//			* ERFC( (x-c)/(SQRT(2)*sigma) + sigma/(SQRT(2)*beta) )
			//	(3)	y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
			
			// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude step, par[4]: gradient of lin bg
			// par[5]: amplitude skewed gaus, par[6]: beta
			
	return PeakFuncGaussian (x,par) + SkewedOnly(x,par);
}

Double_t THypGePeakFitFunction::PeakFuncFreeSkewedPosition (Double_t *x,Double_t *par)
{
			// Peak shape is taken from gf2 http://www.phy.anl.gov/gammasphere/doc/gf2.hlp
			// composed by the sum of:
			//	(1) a Gaussian,
			//	(2) a skewed Gaussian,
			//	(3) a smoothed step function to increase the background on the low-energy
			//			side of the peak. Components (2) and/or (3) can easily be set to zero if not
			//			required, and
			//	(4)	a constant				// not in gf2
			//
			//	(2)	y = constant
			//			* EXP( (x-c)/beta )
			//			* ERFC( (x-c)/(SQRT(2)*sigma) + sigma/(SQRT(2)*beta) )
			//	(3)	y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
			
			// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude step, par[4]: gradient of lin bg
			// par[5]: amplitude skewed gaus, par[6]: beta,  par[7]: x0_s, par[8]: sigma_s  
			
			//Double_t SkewedGausian = par[5] 
			//				* TMath::Exp((x[0]-par[7])/par[6]) 
			//				* TMath::Erfc( ((x[0] - par[7]) /(TMath::Sqrt(2)*par[8]) ) + par[8]/TMath::Sqrt(2)/par[6]); //par[5]: amplitude skewed gaus, par[6]: beta,  par[7]: x0_s, par[8]: sigma_s 
	return PeakFuncGaussian (x,par) + FreeSkewedOnly(x,par);//SkewedGausian;
}

Double_t THypGePeakFitFunction::PeakFuncDoubleGausian (Double_t *x,Double_t *par)
{
			// Peak shape is taken from gf2 http://www.phy.anl.gov/gammasphere/doc/gf2.hlp
			// composed by the sum of:
			//	(1) a Gaussian,
			//	(2) a second Gaussian,
			//	(3) a smoothed step function to increase the background on the low-energy
			//			side of the peak. Components (2) and/or (3) can easily be set to zero if not
			//			required, and
			//	(4)	a constant				// not in gf2
			//
			//	(3)	y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
			
			// par[0]: amplitude gaus, par[1]: x0; par[2]: sigma, par[3]: amplitude step, par[4]: gradient of lin bg
			// par[5]: amplitude second gaus, par[6]: x0_2,  par[7]: sigma_2  
			
	return PeakFuncGaussian (x,par) + SecondGausOnly(x,par);
}

// components
Double_t THypGePeakFitFunction::GausOnly (Double_t *x,Double_t *par)
{
			
			Double_t Gausian = par[0] * TMath::Exp(- 0.5*pow((x[0] - par[1])/par[2],2)); // par[0]: amplitude, par[1]: x[0]0; par[2]: sigma
			
	return Gausian ;
}

Double_t THypGePeakFitFunction::SmoothedStepOnly (Double_t *x,Double_t *par)
{
			Double_t SmoothedStep = par[3] * TMath::Erfc( (x[0] - par[1]) /(TMath::Sqrt(2)*par[2])); // par[3]: amplitude
	return  SmoothedStep;
}

Double_t THypGePeakFitFunction::LinearOnly (Double_t *x,Double_t *par)
{
			Double_t Constant = par[4]*x[0]+par[5]; // par[4]: gradient of lin bg, par[5]: const bg
	return  Constant;
}
Double_t THypGePeakFitFunction::SkewedOnly (Double_t *x,Double_t *par)
{
	Double_t SkewedGausian = par[6] 
					* TMath::Exp((x[0]-par[1])/par[7]) 
					* TMath::Erfc( ((x[0] - par[1]) /(TMath::Sqrt(2)*par[2]) ) + par[2]/TMath::Sqrt(2)/par[7]); // par[6]: amplitude skewed gaus, par[7]: beta
	return SkewedGausian;
}
Double_t THypGePeakFitFunction::FreeSkewedOnly (Double_t *x,Double_t *par)
{
			Double_t SkewedGausian = par[6] 
							* TMath::Exp((x[0]-par[8])/par[7]) 
							* TMath::Erfc( ((x[0] - par[8]) /(TMath::Sqrt(2)*par[9]) ) + par[9]/TMath::Sqrt(2)/par[7]); //par[6]: amplitude skewed gaus, par[7]: beta,  par[8]: x0_s, par[9]: sigma_s 
	return  SkewedGausian;
}

Double_t THypGePeakFitFunction::SecondGausOnly (Double_t *x,Double_t *par)
{
			Double_t SecondGausian = par[6] * TMath::Exp(- 0.5*pow((x[0] - par[7])/par[8],2)); // par[6]: amplitude, par[7]: x0_2; par[8]: sigma_2
	return  SecondGausian;
}
