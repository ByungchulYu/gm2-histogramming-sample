#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#ifndef HISTOGRAM_BASE
#define HISTOGRAM_BASE

// =================================================================================================

// encapsulates a skim TTree entry with relevant branches as member variables for a single positron
class PositronData {
      
  public:

    unsigned int gpsInteger;
    double time; // cluster time in clock ticks
    double energy; // cluster energy in MeV

    double x; // calorimeter x position in crystal widths
    double y; // calorimeter y position in crystal widths

    int caloIndex;
    int runIndex;
    int subrunIndex;
    int fillIndex;
    int bunchNumber;

    bool laserInFill; // FALSE -> indicates positron is from non-laser fill

    // std::vector<double>* inFillGain = 0;
    // std::vector<double>* crystalEnergy = 0;

};

// =================================================================================================

// encapsulates a skim TTree entry with relevant branches as member variables for a constructed pileup event
class PileupData {

  public:

    bool laserInFill; // FALSE -> indicates positron is from non-laser fill

    std::vector<int> pileupIndex;
    std::vector<bool> pileupFlagged;

    std::vector<double> pileupTime; // cluster times in clock ticks
    std::vector<double> pileupEnergy; // cluster energies in MeV
    
    std::vector<double> pileupX; // calorimeter x positions in crystal widths
    std::vector<double> pileupY; // calorimeter y positions in crystal widths

    std::vector<int> pileupCaloIndex;
    int runIndex;
    int subrunIndex;
    int	fillIndex;
    int	bunchNumber;

};

// =================================================================================================

// encapsulates external inputs needed for constructing the lost muon histogram
class LostMuonInput {

  public:

    TH1D* timeOfFlight; // expected lost muon times-of-flight for each calorimeter to the next
    std::vector<TH1D*> caloEfficiency; // positron intensity (wiggle histograms) per calorimeter
    double events; // normalization factor for efficiency histograms

};

// =================================================================================================

// encapsulates a skim TTree entry with relevant branches as member variables for a lost muon candidate
class LostMuonData {

  public:

    bool laserInFill;

    int calo1;
    double time1;
    double energy1;
    double x1;
    double y1;

    std::vector<double> calo2times;
    std::vector<double> calo3times;
    std::vector<double> calo4times;

    std::vector<double> calo2energies;
    std::vector<double> calo3energies;
    std::vector<double> calo4energies;

    std::vector<double> calo2x;
    std::vector<double> calo3x;
    std::vector<double> calo4x;

    std::vector<double> calo2y;
    std::vector<double> calo3y;
    std::vector<double> calo4y;

    int fillIndex;
    int subrunIndex;
    int runIndex;
    int bunchNumber;

};

// =================================================================================================

class HistogramBase {

  public:

    virtual ~HistogramBase() {};

    // use this method to initialize all histogram objects
    virtual void bookHistograms(int seedIndex, int skimIndex) = 0;

    // this method will be called once for every entry in the singles TTree
    // PositronData object contains all relevant branches from the singles TTree entry as members
    // e.g. entry.time, entry.energy, ...
    virtual void fillSinglesHistograms(PositronData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) = 0;

    // this method will be called once for every entry in the double-pileup TTree
    // PileupData object contains all relevant branches from the double-pileup TTree entry as members
    // e.g. entry.pileupTime, entry.pileupEnergy, ...
    virtual void fillDoublesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) = 0;

    // this method will be called once for every entry in the triple-pileup TTree
    // PileupData object contains all relevant branches from the triple-pileup TTree entry as members
    // e.g. entry.pileupTime, entry.pileupEnergy, ...
    virtual void fillTriplesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) = 0;

    // this method will be called once for every entry in the lost-muon-candidate TTree
    // LostMuonData object contains all relevant branches from the lost muon TTree entry as members
    // LostMuonInput object contains the expected lost muon times-of-flight per-calorimeter
    virtual void fillLostMuonHistograms(LostMuonData& entry, LostMuonInput& lmInput, double frRandomization, double vwRandomization, int seedIndex, int skimIndex) = 0;

    // use this method to call Write() on all histograms (and delete, if pointers)
    virtual void writeHistograms(TFile* outputFile, int seedIndex) = 0;

};

#endif
