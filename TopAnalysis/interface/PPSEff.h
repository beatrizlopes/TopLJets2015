// ROOT include
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

// cpp include
#include <iostream>
#include <string>

enum class Periods {
	undefined = 0,
	era2017B  = 1,
	era2017C1 = 2,
	era2017C2 = 3,
	era2017D  = 4,
	era2017E  = 5,
	era2017F1 = 6,
	era2017F2 = 7,
	era2017F3 = 8,
	Count     = 9
};



class PPSEff
{

	public: 
	
	//Constructors
	PPSEff(){} 
	PPSEff(std::string filename){
        init(filename);  
	}
	
	// Destructor 
	~PPSEff() { if(RecoInit || EffInit) _file0->Close();}
	
	int GetRunIndx(unsigned int runNumber){
		if(runNumber>=297050&& runNumber<=299329) return static_cast<int>(Periods::era2017B);
		if(runNumber>=299368&& runNumber<=300780) return static_cast<int>(Periods::era2017C1);
		if(runNumber>=300806&& runNumber<=302029) return static_cast<int>(Periods::era2017C2);
		if(runNumber>=302031&& runNumber<=302663) return static_cast<int>(Periods::era2017D);
		if(runNumber>=303825&& runNumber<=304797) return static_cast<int>(Periods::era2017E);
		if(runNumber>=305044&& runNumber<=305114) return static_cast<int>(Periods::era2017F1);
		if(runNumber>=305178&& runNumber<=305902) return static_cast<int>(Periods::era2017F2);
		if(runNumber>=305967&& runNumber<=306460) return static_cast<int>(Periods::era2017F3);
		//if(runNumber>=306936&& runNumber<=307082) return static_cast<int>(Periods::era2017H);
		return static_cast<int>(Periods::undefined);
	}

	TString era[static_cast<int>(Periods::Count)] = {
		"UNDEF","2017B","2017C1","2017C2","2017D","2017E","2017F1","2017F2","2017F3"
	};
	
	float getEff(float x, float y, int arm, unsigned int runNumber){
		if(!EffInit){ std::cout << "ERROR: call getEff() w/o proper initialization" << std::endl; return 0;}
		if(arm==0){
			int ibin = _heff201745[0]->FindBin(x,y);
			return _heff201745[GetRunIndx(runNumber)]->GetBinContent(ibin);
		}
		else if(arm==1){
			int ibin = _heff201756[0]->FindBin(x,y);
			return _heff201756[GetRunIndx(runNumber)]->GetBinContent(ibin);
		}
		else{ std::cout <<" wrong arm number (expect 0 or 1), return 0...\n"; return 0;}
	}

	float getEffErr(float x, float y, int arm, unsigned int runNumber){ // returns stat + sys error
		if(!EffInit){ std::cout << "ERROR: call getEffErr() w/o proper initialization" << std::endl; return 0;}
		float val = 0 , err = 0;
		if(arm==0){
			int ibin = _heff201745[0]->FindBin(x,y);
            val = _heff201745[GetRunIndx(runNumber)]->GetBinContent(ibin);
			err = _heff201745[GetRunIndx(runNumber)]->GetBinError(ibin);
			return sqrt( err*err + (val*addUncFlat)*(val*addUncFlat) );
		}
		else if(arm==1){
			int ibin = _heff201756[0]->FindBin(x,y);
			val = _heff201745[GetRunIndx(runNumber)]->GetBinContent(ibin);
			err = _heff201745[GetRunIndx(runNumber)]->GetBinError(ibin);
			return sqrt( err*err + (val*addUncFlat)*(val*addUncFlat) );
		}
		else{ std::cout <<" wrong arm number (expect 0 or 1), return 0...\n"; return 0;}
	}
	
	float getRelEffErrSq(float x, float y, int arm, unsigned int runNumber){ // returns squared relative (stat + sys) error
		if(!EffInit){ std::cout << "ERROR: call getEffErr() w/o proper initialization" << std::endl; return 0;}
		float val = 0 , err = 0;
		if(arm==0){
			int ibin = _heff201745[0]->FindBin(x,y);
			val = _heff201745[GetRunIndx(runNumber)]->GetBinContent(ibin);
			err = _heff201745[GetRunIndx(runNumber)]->GetBinError(ibin);
            return (val>0) ? ((err/val)*(err/val) + addUncFlat*addUncFlat) : 0;
		}
		else if(arm==1){
			int ibin = _heff201756[0]->FindBin(x,y);
			val = _heff201756[GetRunIndx(runNumber)]->GetBinContent(ibin);
			err = _heff201756[GetRunIndx(runNumber)]->GetBinError(ibin);
            return (val>0) ? ((err/val)*(err/val) + addUncFlat*addUncFlat) : 0;
		}
		else{ std::cout <<" wrong arm number (expect 0 or 1), return 0...\n"; return 0;}
	}	
	
	float getRecoErr(float xi, int arm, unsigned int runNumber){
		if(!RecoInit){ std::cout << "ERROR: call getRecoErr() w/o proper initialization" << std::endl; return 0;}
		if(arm==0){
			if(runNumber>=297050&& runNumber<=302663) return _reco_err45[0]->Eval(xi);
			if(runNumber>=303825&& runNumber<=306460) return _reco_err45[1]->Eval(xi);
			return 0;
		}
		else if(arm==1){
			if(runNumber>=297050&& runNumber<=302663) return _reco_err56[0]->Eval(xi);
			if(runNumber>=303825&& runNumber<=306460) return _reco_err56[1]->Eval(xi);
			return 0;
		}
		else{ std::cout <<" wrong arm number (expect 0 or 1), return 0...\n"; return 0;}
	}
	
	float getXiHigh(int arm, unsigned int runNumber, float xangle){
		// from https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsGettingStarted#Fiducial_cuts
		if(xangle<120 || xangle>150) return 0; // sanity check
		if(arm<100){ // arm = 0
			if(runNumber<303825) return 0.066 + (3.54E-4 * xangle);
			if(runNumber>=303825) return 0.073 + (4.11E-4 * xangle);
			return 0;
		}
		else {
			if(runNumber<303825) return 0.062 + (5.96E-4 * xangle);
			if(runNumber>=303825) return 0.067 + (6.87E-4 * xangle);
			return 0;
		}
	}
	
	float getThXStarHigh(int arm, unsigned int runNumber, float xi, float xangle){
    // from https://twiki.cern.ch/twiki/bin/view/CMS/TaggedProtonsFiducialCuts
        if(xangle<120 || xangle>150) return 999; // sanity check
		float thX_max = 0;
        if(arm<100){ // arm = 0
            if(runNumber<303825) thX_max= -(8.71198E-07*xangle-0.000134726)+((xi<(0.000264704*xangle+0.081951))*-(4.32065E-05*xangle-0.0130746)+(xi>=(0.000264704*xangle+0.081951))*-(0.000183472*xangle-0.0395241))*(xi-(0.000264704*xangle+0.081951));
            else thX_max= -(8.92079E-07*xangle-0.000150214)+((xi<(0.000278622*xangle+0.0964383))*-(3.9541e-05*xangle-0.0115104)+(xi>=(0.000278622*xangle+0.0964383))*-(0.000108249*xangle-0.0249303))*(xi-(0.000278622*xangle+0.0964383));
        }
        else {
            if(runNumber<303825) thX_max= 3.43116E-05+((xi<(0.000626936*xangle+0.061324))*0.00654394+(xi>=(0.000626936*xangle+0.061324))*-(0.000145164*xangle-0.0272919))*(xi-(0.000626936*xangle+0.061324));
            else thX_max= 4.56961E-05+((xi<(0.00075625*xangle+0.0643361))*-(3.01107e-05*xangle-0.00985126)+(xi>=(0.00075625*xangle+0.0643361))*-(8.95437e-05*xangle-0.0169474))*(xi-(0.00075625*xangle+0.0643361));
        }
		return -thX_max;
    }
	
	bool InFiducial(int arm, int det, unsigned int runNumber, float x, float y){
		TString era("");
		if(runNumber>=297050&& runNumber<=299329) era="2017B";
		if(runNumber>=299368&& runNumber<=300780) era="2017C1";
		if(runNumber>=300806&& runNumber<=302029) era="2017C2";
		if(runNumber>=302031&& runNumber<=302663) era="2017D";
		if(runNumber>=303825&& runNumber<=304797) era="2017E";
		if(runNumber>=305044&& runNumber<=305114) era="2017F1";
		if(runNumber>=305178&& runNumber<=305902) era="2017F2";
		if(runNumber>=305967&& runNumber<=306460) era="2017F3";
		if(runNumber>=306936&& runNumber<=307082) era="2017H";

		std::map<std::pair<int, int>, double> XLow_, XHigh_, YLow_, YHigh_;
		fillFiducialCutsVectors(era, XLow_, XHigh_, YLow_, YHigh_);
		std::pair<int, int> rp(arm,det);
		return (x>XLow_[rp]) && (x<XHigh_[rp]) && (y>YLow_[rp]) && (y<YHigh_[rp]);
	}	
	
	private: 
	float addUncFlat=0; // additional flat systematic uncertainty (see L167 for more info)
	const static int Nperiods = static_cast<int>(Periods::Count);
	int _nbins; float _bw, _xmin=0, _xmax=0.2;
	bool RecoInit = false, EffInit = false, is2D = false;
	TFile * _file0 = NULL;
	TH2D *_heff201745[Nperiods], *_heff201756[Nperiods];
	TGraphErrors * _reco_err45[2], * _reco_err56[2];
	void init(std::string filename){
		TString hname;
		if(filename.empty()) {
			std::cout << "Warning: PPSEff::init() No initialization file specified" << std::endl;
			return;
		}
		std::cout << "INFO: PPSEff::init() Reads eff histogram from "<<filename.c_str() << std::endl;
		_file0 = TFile::Open(filename.c_str());
		if(!_file0) std::cout  <<"ERROR: PPSEff::init() no such file "<< filename.c_str() << std::endl;
		
		if(filename.find("reco")!=std::string::npos){
			_reco_err45[0] = (TGraphErrors *)_file0->Get("2017_preTS2/multi rp-0/xi/g_systematics_vs_xi");
			_reco_err45[1] = (TGraphErrors *)_file0->Get("2017_postTS2/multi rp-0/xi/g_systematics_vs_xi");
			_reco_err56[0] = (TGraphErrors *)_file0->Get("2017_preTS2/multi rp-1/xi/g_systematics_vs_xi");
			_reco_err56[1] = (TGraphErrors *)_file0->Get("2017_postTS2/multi rp-1/xi/g_systematics_vs_xi");

			RecoInit=true;
		}
		else if(filename.find("1D2DMultiTrack")!=std::string::npos){
			// load the strips efficiency histograms
			for(int i=1;i<Nperiods;i++){
				hname = Form("Strips/2017/%s/h45_%s_all_2D",era[i].Data(),era[i].Data());
				if(hname.Contains("2017E")) hname+="_restrictedJSON"; // correction discussed in https://indico.cern.ch/event/1022654/
				_heff201745[i]  = (TH2D *)_file0->Get(hname);
				if(!_heff201745[i]) {std::cout << "ERROR init() " << hname<< " not loaded " << std::endl; return;}
				hname = Form("Strips/2017/%s/h56_%s_all_2D",era[i].Data(),era[i].Data());
				if(hname.Contains("2017E")) hname+="_restrictedJSON"; // correction discussed in https://indico.cern.ch/event/1022654/
				_heff201756[i]  = (TH2D *)_file0->Get(hname);
				if(!_heff201756[i]) {std::cout << "ERROR init() " << hname<< " not loaded " << std::endl; return;}				
			}
			_heff201745[0] = (TH2D *)_heff201745[1]->Clone("empty"); _heff201745[0]->Reset();
			_heff201756[0] = (TH2D *)_heff201756[1]->Clone("empty"); _heff201756[0]->Reset();
			EffInit=true;
			// set Flat unc from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TaggedProtonsStripsEfficiencies#2017
			addUncFlat = sqrt(0.01*0.01+0.01*0.01+0.02*0.02);
		}
		else if(filename.find("pixelEfficiencies_multiRP")!=std::string::npos){
			// load the strips efficiency histograms
			for(int i=1;i<Nperiods;i++){
				hname=Form("Pixel/2017/%s/h45_220_%s_all_2D",era[i].Data(),era[i].Data());
				_heff201745[i]  = (TH2D *)_file0->Get(hname);
				if(!_heff201745[i]) {std::cout << "ERROR init() " << hname<< " not loaded " << std::endl; return;}				
				hname = Form("Pixel/2017/%s/h56_220_%s_all_2D",era[i].Data(),era[i].Data());
				_heff201756[i]  = (TH2D *)_file0->Get(hname);
				if(!_heff201756[i]) {std::cout << "ERROR init() " << hname<< " not loaded " << std::endl; return;}								
			}
			_heff201745[0] = (TH2D *)_heff201745[1]->Clone("empty"); _heff201745[0]->Reset();
			_heff201756[0] = (TH2D *)_heff201756[1]->Clone("empty"); _heff201756[0]->Reset();
			EffInit=true;
			addUncFlat = 0;
		}
		else{
			std::cout << "ERROR: PPSEff::init() Wrong initialization file specified - " << filename << std::endl;
		}
		//_nbins = _hef201756[0]->GetNbinsX();
		//_bw    = _hef201756[0]->GetBinWidth(1);
		//_xmin  = _hef201756[0]->GetBinCenter(1) - 0.5 * _bw;
		//_xmax  = _hef201756[0]->GetBinCenter(_nbins) + 0.5 * _bw;	
	}
    
	
    void fillFiducialCutsVectors( // see https://twiki.cern.ch/twiki/bin/viewauth/CMS/TaggedProtonsPixelEfficiencies
		TString era, std::map<std::pair<int, int>, double> &fiducialXLow_,
		std::map<std::pair<int, int>, double> &fiducialXHigh_,
		std::map<std::pair<int, int>, double> &fiducialYLow_,
		std::map<std::pair<int, int>, double> &fiducialYHigh_) {
		if (era == "2017B") {
		fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.098;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.698;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.698;
		} else if ((era == "2017C1")|| (era == "2017C2") || (era == "2017C") || (era == "2017D")) {
		fiducialXLow_[std::pair<int, int>(0, 2)] = 1.860;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.334;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.098;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.698;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.698;
		} else if (era == "2017E") {
		fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -10.098;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.998;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -9.698;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.398;
		} else if ((era == "2017F1")||(era == "2017F2")||(era == "2017F3")||(era == "2017F")||(era == "2017H")) {
		fiducialXLow_[std::pair<int, int>(0, 2)] = 1.995;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.479;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -10.098;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.998;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.422;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -9.698;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.398;
		} else if (era == "2018A") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.710;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.927;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -11.598;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 3.698;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -10.898;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.398;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 18.498;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -11.298;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.098;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.420;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 20.045;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.398;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.098;
		} else if (era == "2018B1") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.850;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.927;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -11.598;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 3.698;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.420;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -10.798;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 4.298;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 18.070;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -11.198;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.098;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.420;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 25.045;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.398;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 5.098;
		} else if (era == "2018B2") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.562;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.640;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.198;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.135;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.798;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.498;
		} else if (era == "2018C") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.564;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.930;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.198;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.698;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.398;
		} else if (era == "2018D1") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.847;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.930;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -11.098;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.098;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.398;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.698;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 4.398;
		} else if (era == "2018D2") {
		fiducialXLow_[std::pair<int, int>(0, 0)] = 2.847;
		fiducialXHigh_[std::pair<int, int>(0, 0)] = 17.931;
		fiducialYLow_[std::pair<int, int>(0, 0)] = -10.598;
		fiducialYHigh_[std::pair<int, int>(0, 0)] = 4.498;
		fiducialXLow_[std::pair<int, int>(0, 2)] = 2.278;
		fiducialXHigh_[std::pair<int, int>(0, 2)] = 24.620;
		fiducialYLow_[std::pair<int, int>(0, 2)] = -11.598;
		fiducialYHigh_[std::pair<int, int>(0, 2)] = 3.398;
		fiducialXLow_[std::pair<int, int>(1, 0)] = 3.000;
		fiducialXHigh_[std::pair<int, int>(1, 0)] = 17.931;
		fiducialYLow_[std::pair<int, int>(1, 0)] = -10.498;
		fiducialYHigh_[std::pair<int, int>(1, 0)] = 4.698;
		fiducialXLow_[std::pair<int, int>(1, 2)] = 2.279;
		fiducialXHigh_[std::pair<int, int>(1, 2)] = 24.760;
		fiducialYLow_[std::pair<int, int>(1, 2)] = -10.598;
		fiducialYHigh_[std::pair<int, int>(1, 2)] = 3.898;
		} else
		  std::cout << "WARNING: Era not recognized!" << std::endl;
		return;
	}	
	
};

