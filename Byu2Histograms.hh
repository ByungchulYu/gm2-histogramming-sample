#ifndef BYU2_HISTOGRAMS_HH
#define BYU2_HISTOGRAMS_HH

#include <iostream>

#include "HistogramBase.hh"

// ROOT libraries.
#include <TTree.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"

static const double ct2us = 1.25/1000;    // Conversion factor from clock tick to microsecond.
static const double cyclotronPeriod = 0.1492; // microseconds

class Byu2Histograms: public HistogramBase {

public:

    // Constructor.
    Byu2Histograms();

    void bookHistograms(int seedIndex, int skimIndex) override;
    void fillSinglesHistograms(PositronData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillDoublesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillTriplesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void fillLostMuonHistograms(LostMuonData& entry, LostMuonInput& lmInput, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) override;
    void writeHistograms(TFile* outputFile, int seedIndex) override;

private:

	double   			t_min; 			     	// 0 us
	double   			t_max;					// 700 us rounding issue가 있어서 뒤에서 재정의됨
	int      			t_n_bins;				// 700/0.1492 = 4691 (0.1492us = bin width)

	double   			E_min;					// 1050 MeV
	double   			E_max;					// 3060 MeV
	int     			E_n_bins;				// 30 = (2010/67)   
	int      			E_bin_width;			// 67 MeV	

	int 				prev_runIndexS_;		// runIndex in singlefill function
	int 				prev_runIndexD_;		// runIndex in doublefill function
	int 				prev_runIndexH_;		// runIndex in higherfill function

	int 				prev_subrunIndexS_;		// subrunIndex in singlefill function
	int 				prev_subrunIndexD_;		// subrunIndex in doublefill function
	int 				prev_subrunIndexH_;		// subrunIndex in higherfill function

    TH2F*    			EvsT_;					// raw ET histogram
	TH2F*    			EvsT_D_;				// double PU histogram : PileupIndex == 2    (+0.5 weight for PU | -0.5 weight for PC)
	TH2F*    			EvsT_H_;				// higher PU histogram : PileupIndex == pu3  (+0.5 weight for PU | -0.5 weight for PC)
	TH2F*				EvsT_PU_;				// total  PU histogram

	TTree*				TREE_ET_aux_;	 		// Tree for subrun level information
	TTree*				TREE_ET_;	    		// Tree for subrun level information
	
	TBranch* 			prev_runindex_S;		// runIndex Branch in singlefill function
	TBranch* 			prev_runindex_D;		// runIndex Branch in doublefill function
	TBranch* 			prev_runindex_H;		// runIndex Branch in higherfill function

	TBranch* 			prev_index_S;			// subrunIndex Branch in singlefill function
	TBranch* 			prev_index_D;			// subrunIndex Branch in doublefill function
	TBranch* 			prev_index_H;			// subrunIndex Branch in higherfill function

	TBranch* 			EvsT_branch;			// singlefill ET histogram
	TBranch* 			EvsT_D_branch;			// doublefill ET histogram	
	TBranch* 			EvsT_H_branch;			// higherfill ET histogram
	TBranch* 			EvsT_PU_branch;			// pileupfill ET histogram

	double				subruntimeindex_;
	TBranch*			subruntime_;
	std::vector<unsigned int>	timestamps_;			// Unixtimestamp vector
	std::vector<double>	ave_time_vec_;			// subrun average time vector
};

#endif