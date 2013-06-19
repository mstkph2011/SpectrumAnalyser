/*
 * THypGePeakFitFunction.h
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
#ifndef THYPGEPEAKFITFUNCTION_H
#define THYPGEPEAKFITFUNCTION_H

#include "TMath.h"

class THypGePeakFitFunction
{
	public:
		THypGePeakFitFunction();
		Double_t 			PeakFunc (Double_t *x,Double_t *par);
		Double_t 			PeakFuncGaussian (Double_t *x,Double_t *par);
		Double_t			PeakFuncFreeSkewedPosition (Double_t *x,Double_t *par);
		Double_t			PeakFuncDoubleGausian (Double_t *x,Double_t *par);
		
		Double_t			GausOnly (Double_t *x,Double_t *par);
		Double_t			SmoothedStepOnly (Double_t *x,Double_t *par);
		Double_t			LinearOnly (Double_t *x,Double_t *par);
		Double_t			SkewedOnly(Double_t *x,Double_t *par);
		Double_t			FreeSkewedOnly(Double_t *x,Double_t *par);
		Double_t			SecondGausOnly (Double_t *x,Double_t *par);
		
	private:
		/* add your private declarations */
};

#endif /* THYPGEPEAKFITFUNCTION_H */ 
