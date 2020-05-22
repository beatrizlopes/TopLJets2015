// ROOT include
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

// cpp include
#include <iostream>

using namespace std; 
class PixEff
{
	public: 
	
	//Constructor 
	PixEff(){
		std::cout << "Reads eff histogram from "<<eff_file.Data()<<std::endl;
		_file0 = TFile::Open(eff_file);
		
		_hef201745[0] = (TH1D *)_file0->Get("Pixel/2017/2017B/h45_220_2017B_all_1D");
		_hef201756[0] = (TH1D *)_file0->Get("Pixel/2017/2017B/h56_220_2017B_all_1D");
		_hef201745[1] = (TH1D *)_file0->Get("Pixel/2017/2017C1/h45_220_2017C1_all_1D");
		_hef201756[1] = (TH1D *)_file0->Get("Pixel/2017/2017C1/h56_220_2017C1_all_1D");
		_hef201745[2] = (TH1D *)_file0->Get("Pixel/2017/2017E/h45_220_2017E_all_1D");
		_hef201756[2] = (TH1D *)_file0->Get("Pixel/2017/2017E/h56_220_2017E_all_1D");
		_hef201745[3] = (TH1D *)_file0->Get("Pixel/2017/2017F1/h45_220_2017F1_all_1D");
		_hef201756[3] = (TH1D *)_file0->Get("Pixel/2017/2017F1/h56_220_2017F1_all_1D");
		
		_nbins = _hef201756[0]->GetNbinsX();
		_bw    = _hef201756[0]->GetBinWidth(1);
		_xmin  = _hef201756[0]->GetBinCenter(1) - 0.5 * _bw;
		_xmax  = _hef201756[0]->GetBinCenter(_nbins) + 0.5 * _bw;
		
	}
	
	// Destructor 
	~PixEff() {}
	
	float getEff(float xi, int rpid, int runNumber){
		if( xi < _xmin || xi > _xmax) return 0;
		int ibin = int( (xi-_xmin)/_bw );
		if(rpid==23){
			return _hef201745[0]->GetBinContent(ibin);
		}
		else if(rpid==123){
			return _hef201756[0]->GetBinContent(ibin);
		}
		else{ cout << "wrong arm number (expect 23 or 123)"<<endl; return 0;}
	}
	
	
	private: 
	int _nbins; float _bw, _xmin, _xmax;
	TFile * _file0;
	TString eff_file = "/eos/project/c/ctpps/subsystems/Pixel/RPixTracking/pixelEfficiencies.root";
	int _year = 2017;
	TH1D *_hef201745[4], *_hef201756[4];
	
	
};

