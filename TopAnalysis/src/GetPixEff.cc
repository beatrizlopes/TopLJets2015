// ROOT include
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"

// cpp include
#include <iostream>

using namespace std; 

enum class Periods {
	era2017B  = 0,
	era2017C1 = 1,
	era2017C2 = 2,
	era2017D  = 3,
	era2017E  = 4,
	era2017F1 = 5,
	era2017F2 = 6,
	era2017F3 = 7,
};

class PixEff
{

	public: 
	
	//Constructor 
	PixEff(){
		std::cout << "Reads eff histogram from "<<eff_file.Data()<<std::endl;
		_file0 = TFile::Open(eff_file);
		
		_hef201745[static_cast<int>(Periods::era2017B)]  = (TH1D *)_file0->Get("Pixel/2017/2017B/h45_220_2017B_all_1D");
		_hef201756[static_cast<int>(Periods::era2017B)]  = (TH1D *)_file0->Get("Pixel/2017/2017B/h56_220_2017B_all_1D");
		_hef201745[static_cast<int>(Periods::era2017C1)] = (TH1D *)_file0->Get("Pixel/2017/2017C1/h45_220_2017C1_all_1D");
		_hef201756[static_cast<int>(Periods::era2017C1)] = (TH1D *)_file0->Get("Pixel/2017/2017C1/h56_220_2017C1_all_1D");
		_hef201745[static_cast<int>(Periods::era2017C2)] = (TH1D *)_file0->Get("Pixel/2017/2017C2/h45_220_2017C2_all_1D");
		_hef201756[static_cast<int>(Periods::era2017C2)] = (TH1D *)_file0->Get("Pixel/2017/2017C2/h56_220_2017C2_all_1D");
		_hef201745[static_cast<int>(Periods::era2017D)]  = (TH1D *)_file0->Get("Pixel/2017/2017D/h45_220_2017D_all_1D");
		_hef201756[static_cast<int>(Periods::era2017D)]  = (TH1D *)_file0->Get("Pixel/2017/2017D/h56_220_2017D_all_1D");
		_hef201745[static_cast<int>(Periods::era2017E)]  = (TH1D *)_file0->Get("Pixel/2017/2017E/h45_220_2017E_all_1D");
		_hef201756[static_cast<int>(Periods::era2017E)]  = (TH1D *)_file0->Get("Pixel/2017/2017E/h56_220_2017E_all_1D");
		_hef201745[static_cast<int>(Periods::era2017F1)] = (TH1D *)_file0->Get("Pixel/2017/2017F1/h45_220_2017F1_all_1D");
		_hef201756[static_cast<int>(Periods::era2017F1)] = (TH1D *)_file0->Get("Pixel/2017/2017F1/h56_220_2017F1_all_1D");
		_hef201745[static_cast<int>(Periods::era2017F2)] = (TH1D *)_file0->Get("Pixel/2017/2017F2/h45_220_2017F2_all_1D");
		_hef201756[static_cast<int>(Periods::era2017F2)] = (TH1D *)_file0->Get("Pixel/2017/2017F2/h56_220_2017F2_all_1D");
		_hef201745[static_cast<int>(Periods::era2017F3)] = (TH1D *)_file0->Get("Pixel/2017/2017F3/h45_220_2017F3_all_1D");
		_hef201756[static_cast<int>(Periods::era2017F3)] = (TH1D *)_file0->Get("Pixel/2017/2017F3/h56_220_2017F3_all_1D");
		
		_nbins = _hef201756[0]->GetNbinsX();
		_bw    = _hef201756[0]->GetBinWidth(1);
		_xmin  = _hef201756[0]->GetBinCenter(1) - 0.5 * _bw;
		_xmax  = _hef201756[0]->GetBinCenter(_nbins) + 0.5 * _bw;
		
	}
	
	// Destructor 
	~PixEff() {_file0->Close();}
	
	float getEff(float xi, int rpid, int runNumber){
		if( xi < _xmin || xi > _xmax) return 0;
		int ibin = int( (xi-_xmin)/_bw ) + 1;
		if(rpid==23){
			if(runNumber>=297050&& runNumber<=299329) return _hef201745[static_cast<int>(Periods::era2017B)]->GetBinContent(ibin);
			if(runNumber>=299368&& runNumber<=300780) return _hef201745[static_cast<int>(Periods::era2017C1)]->GetBinContent(ibin);
			if(runNumber>=300806&& runNumber<=302029) return _hef201745[static_cast<int>(Periods::era2017C2)]->GetBinContent(ibin);
			if(runNumber>=302031&& runNumber<=302663) return _hef201745[static_cast<int>(Periods::era2017D)]->GetBinContent(ibin);
			if(runNumber>=303825&& runNumber<=304797) return _hef201745[static_cast<int>(Periods::era2017E)]->GetBinContent(ibin);
			if(runNumber>=305044&& runNumber<=305114) return _hef201745[static_cast<int>(Periods::era2017F1)]->GetBinContent(ibin);
			if(runNumber>=305178&& runNumber<=305902) return _hef201745[static_cast<int>(Periods::era2017F2)]->GetBinContent(ibin);
			if(runNumber>=305967&& runNumber<=306460) return _hef201745[static_cast<int>(Periods::era2017F3)]->GetBinContent(ibin);
			return 0;
		}
		else if(rpid==123){
			if(runNumber>=297050&& runNumber<=299329) return _hef201756[static_cast<int>(Periods::era2017B)]->GetBinContent(ibin);
			if(runNumber>=299368&& runNumber<=300780) return _hef201756[static_cast<int>(Periods::era2017C1)]->GetBinContent(ibin);
			if(runNumber>=300806&& runNumber<=302029) return _hef201756[static_cast<int>(Periods::era2017C2)]->GetBinContent(ibin);
			if(runNumber>=302031&& runNumber<=302663) return _hef201756[static_cast<int>(Periods::era2017D)]->GetBinContent(ibin);
			if(runNumber>=303825&& runNumber<=304797) return _hef201756[static_cast<int>(Periods::era2017E)]->GetBinContent(ibin);
			if(runNumber>=305044&& runNumber<=305114) return _hef201756[static_cast<int>(Periods::era2017F1)]->GetBinContent(ibin);
			if(runNumber>=305178&& runNumber<=305902) return _hef201756[static_cast<int>(Periods::era2017F2)]->GetBinContent(ibin);
			if(runNumber>=305967&& runNumber<=306460) return _hef201756[static_cast<int>(Periods::era2017F3)]->GetBinContent(ibin);
			return 0;
		}
		else{ cout << "wrong arm number (expect 23 or 123)"<<endl; return 0;}
	}
	
	
	private: 
	int _nbins; float _bw, _xmin, _xmax;
	TFile * _file0;
	TString eff_file = "/eos/project/c/ctpps/subsystems/Pixel/RPixTracking/pixelEfficiencies_multiRP.root";
	int _year = 2017;
	TH1D *_hef201745[8], *_hef201756[8];
	
	
};

