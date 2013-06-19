#include "TH1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TRint.h"							// starts the interactive mode of ROOT (needed for canvas etc.)

#include "THypGeSpectrumAnalyser.h"

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main(int argc, char* argv[] )
{
	Int_t DataInput =3; //choose way of creating data spectrum: 1 = generate with rand ; 2 = read file
	TRint *App = new TRint("ROOT",0,0); //&argc,argv);
	TCanvas *c1 = new TCanvas("c1","Canvas",600,400);
	c1-> cd();
	TH1D *hSpec;
	TString InputFile;
	
	char buffer[100];
	
	if(DataInput==1)
	{
		hSpec = new TH1D("hSpec","hSpec",1000,1,101);
		
		TRandom *rand = new TRandom(123);
		for (Int_t i = 1; i <100001; i++)
		{
			hSpec->Fill(rand->Gaus(50,2));
			hSpec->Fill(rand->Gaus(10,1));
			if (i < 1001)
			{
				hSpec->Fill(i/10.-0.05,200);
				//cout << hSpec->GetBinWidth(i) <<endl;
			}
		}
		hSpec->Draw();
	}
	Int_t Channel, Data;
	if(DataInput ==2)
	{
		 InputFile= "Co60_290811.csv";
		
		if (argc > 1)
		{
			cout << "Argc: " <<argc << "\t" << argv[1] <<endl;
			InputFile = argv[1];
		}
		else
			cout << "Argc: 1" << endl;
		hSpec = new TH1D("hSpec",InputFile.Data(),8192,1,8192);		// The name of the spectrum is send to the Txt output file
		ifstream DataFile/*("Co60_290811.csv");//*/(InputFile.Data());
		if(!DataFile)
		{
			cout << "ERROR - file does NOT exist" << endl;
			return -1;
		}
		else
		{
			cout << "Start reading" << endl;
			Int_t i =0;
			while(DataFile.good())
			{
				if (i >= 9)
				{
					DataFile.getline(buffer,100);
					sscanf(buffer,"%d;%d",&Channel,&Data);
					//cout << i << "\t"<< Channel << "\t" << Data<<endl;;
					hSpec->Fill(Channel,Data);
				}
				else
				{
					DataFile.getline(buffer,100);
					//cout << i <<"\t" << buffer << endl;
				}
				i++;
			}
		}
	}
	if(DataInput==3)
	{
		InputFile= "Co60_spec.root";
		
		if (argc > 1)
		{
			cout << "Argc: " <<argc << "\t" << argv[1] <<endl;
			InputFile = argv[1];
		}
		
		TFile *InFile= new TFile(InputFile);
		InFile->GetObject("Energy",hSpec);
		//hSpec = new TH1D("hSpec","hSpec",1000,1,101);
		
		hSpec->Draw();
	}
	
	TString Path, InfileName,RootFilename, TxtFilename;
	Int_t LengthOfPath;
	if (DataInput > 1)
	{
		// get path and filename from input
		LengthOfPath =InputFile.Last('/') +1;
		Path = InputFile;
		cout << Path.Remove(LengthOfPath,Path.Length()-LengthOfPath) << endl;
		InfileName = InputFile;
		cout <<InfileName.Remove(0,LengthOfPath) << endl;
		//cout << InputFile.Remove(0,InputFile.Last('/')+1) << " "<< InputFile << endl;
		InfileName.Remove(InfileName.Length()-5, 5);
		
		RootFilename = "Ana_";
		RootFilename += InfileName;
		RootFilename += ".root";
		
		TxtFilename = "Ana_";
		TxtFilename += InfileName;
		TxtFilename += ".txt";
		
		cout << TxtFilename << endl;
	}
	
	
	//hSpec->GetXaxis()->SetRange(2000,8192);
	THypGeSpectrumAnalyser *Ana = new THypGeSpectrumAnalyser(hSpec,"co60bg", 10 );
	Ana->SetSearchRange(2000,8192);
	Ana->SetOutputPath(Path);
	Ana->SetTxtFileOutputName(TxtFilename);
	Ana->SetRootFileOutputName(RootFilename);
	Ana->SetGaussianFitting();
	Ana->SetFreeSkewedFitting();
	Ana->SetSecondGausianFitting();
	Ana->AnalyseSpectrum();
	//hSpec->Draw("");
	hSpec->Print();
	
	
	App->Run();
	return 0;
}
