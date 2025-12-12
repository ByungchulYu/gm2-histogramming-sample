#include "HistogramBase.hh"
// #include "CornellHistograms.hh"
#include "Byu2Histograms.hh"
// #include "RatioHistograms.hh"

#include "TTree.h"
#include "TRandom3.h"

#include <map>
#include <set>
#include <iostream>
#include <algorithm>
#include <memory>

// =================================================================================================

static constexpr double frPeriod = 0.1492; // microseconds
static constexpr double vwPeriod = 0.4366; // microseconds

// =================================================================================================

static std::vector<std::string> allowedClassNames = {
  // "CornellHistograms",
  "RatioHistograms",
  "Byu2Histograms"
};

bool validClassName(std::string& className) {
  return std::find(allowedClassNames.begin(), allowedClassNames.end(), className) != allowedClassNames.end();
}

// =================================================================================================

// read command line inputs of the form "./runHistogramming -d dataset -s skimIndex -p skimFilePath -l lostMuonPath -c class1,class2,... -o outputPath"
void parseInputs(int argc, char** argv, std::string& dataset, int& runYear, int& datasetIndex, int& skimIndex, std::string& skimFilePath, std::string& lostMuonPath, std::vector<std::string>& classNames, std::string& outputPath) {

  // list single-char argument keys (':' = value required, '::' = value optional, otherwise flag)
  const char* const options = "d:s:p:l:c:o:";

  // intermediate string to hold argument listing class names
  std::string classNamesArg = "";

  bool done = false;
  while (!done) {

    // returns next option key as char, and reads value as char* into global variable 'optarg'
    const char option = getopt(argc, argv, options);

    switch(option) {
      case -1: // argument key will be -1 when no more arguments found
        done = true;
        break;
      case 'd':
        dataset = optarg;
        break;
      case 's':
        skimIndex = std::atoi(optarg);
        break;
      case 'p':
        skimFilePath = optarg;
        break;
      // case 'l':
      //   lostMuonPath = optarg;
      //   break;
      case 'c':
        classNamesArg = optarg;
        break;
      case 'o':
        outputPath = optarg;
        break;
      default:
        printf("Unrecognized input option '%c'.\n", option);
        std::exit(1);
    }

  }

  // extract the run year and production dataset letters from the dataset name
  runYear = dataset[0] - '0';

  // convert first capital dataset letter from ASCII char value range [65, 90] to [0, 25]
  datasetIndex = dataset[1] - 'A';

  // intermediate variables to parse argument string listing class names
  std::size_t position = 0;
  std::string delimiter = ",";

  // continually find the position of the next delimiter, until no more found
  while ((position = classNamesArg.find(delimiter)) != std::string::npos) {

    // extract the substring up to the next delimiter
    std::string token = classNamesArg.substr(0, position);

    // if the token is a recognized class name, then add it to the list, otherwise exit with error
    if (validClassName(token)) {
      classNames.push_back(token);
    } else {
      printf("HistogramBase subclass '%s' not recognized.\n", token.c_str());
      std::exit(1);
    }

    // remove this token and delimiter, so the string begins with the next token
    classNamesArg.erase(0, position + delimiter.length());

  }

  if (!classNamesArg.empty()) {
    if (validClassName(classNamesArg)) {
      classNames.push_back(classNamesArg);
    } else {
      printf("HistogramBase subclass '%s' not recognized.\n", classNamesArg.c_str());
      std::exit(1);
    }
  }

}

// =================================================================================================

static constexpr int maxDatasetsPerYear = 100;
static constexpr int maxSkimsPerDataset = 20000;
static constexpr int maxSeedsPerSkim = 1000;

static const int maxSeedsPerDataset = maxSkimsPerDataset * maxSeedsPerSkim;
static const int maxSeedsPerYear = maxDatasetsPerYear * maxSeedsPerDataset;

// Computes the number of unique seeds which must be reserved for skim files *before* this one.
// Therefore, the next 100 integers are available for this skim file to use as random seeds.
int getSeedOffset(int runYear, int datasetIndex, int skimIndex) {
  return maxSeedsPerYear * runYear + maxSeedsPerDataset * datasetIndex + maxSeedsPerSkim * skimIndex;
}

// =================================================================================================

long long getUniqueFillIndex(int runIndex, int subrunIndex, int fillIndex) {
  return runIndex * 1000000LL + subrunIndex * 1000LL + fillIndex;
}

// =================================================================================================

int main(int argc, char** argv) {
  // declare variables for inputs: dataset name, skim file index, and list of classes to run
  std::string skimFilePath = "";
  std::string lostMuonPath = "";
  std::string dataset = "";
  int runYear = -1;
  int datasetIndex = -1;
  int skimIndex = -1;
  std::vector<std::string> classNames;
  std::string outputPath = "";

  // parse command line arguments into above variables (modified by reference)
  parseInputs(argc, argv, dataset, runYear, datasetIndex, skimIndex, skimFilePath, lostMuonPath, classNames, outputPath);
  // std::cout << "[Debug] parsed" << std::endl;

  // compute global offset for the batch of 100 unique random seeds this skim file will use
  const int seedOffset = getSeedOffset(runYear, datasetIndex, skimIndex);

  // construct the file paths to the skim file and lost muon file
  // std::string lostMuonPath = Form("/gm2data/cornell/run%d/muonloss_run%d.root", runYear, runYear);
  // skimFilePath = "skimTest.root";
  // std::string lostMuonPath = "lostmuon.root";

  // open the skim file and lost muon file
  TFile* skimFile = new TFile(skimFilePath.c_str(), "READ");
  // TFile* lostMuonFile = new TFile(lostMuonPath.c_str(), "READ");

  // fetch the TTrees from the skim file
  TTree* singlesTree = (TTree*) skimFile -> Get("crystalTreeMaker1EP/ntuple");
  TTree* doublesTree = (TTree*) skimFile -> Get("crystalTreeMaker2EP/ntuple");
  TTree* triplesTree = (TTree*) skimFile -> Get("crystalTreeMaker3EP/ntuple");
  // TTree* lostMuonTree = (TTree*) skimFile -> Get("lostMuonEP/ntuple");

  // preload the TTree entries into vectors in memory
  std::vector<PositronData> positronEntries;
  std::vector<PileupData> doubleEntries;
  std::vector<PileupData> tripleEntries;
  // std::vector<LostMuonData> lostMuonEntries;

  // ===============================================================================================

  // create dummy positron data object to hold data from current TTree entry
  PositronData tempPositronEntry;

  // point the singles TTree branches to the member variables in the PositronData object
  // singlesTree -> SetBranchAddress("laserInFill", &(tempPositronEntry.laserInFill));
  singlesTree -> SetBranchAddress("gpsInteger", &(tempPositronEntry.gpsInteger));
  singlesTree -> SetBranchAddress("time", &(tempPositronEntry.time));
  singlesTree -> SetBranchAddress("energy", &(tempPositronEntry.energy));
  singlesTree -> SetBranchAddress("x", &(tempPositronEntry.x));
  singlesTree -> SetBranchAddress("y", &(tempPositronEntry.y));
  singlesTree -> SetBranchAddress("caloIndex", &(tempPositronEntry.caloIndex));
  singlesTree -> SetBranchAddress("runIndex", &(tempPositronEntry.runIndex));
  singlesTree -> SetBranchAddress("subrunIndex", &(tempPositronEntry.subrunIndex));
  singlesTree -> SetBranchAddress("fillIndex", &(tempPositronEntry.fillIndex));
  singlesTree -> SetBranchAddress("bunchNumber", &(tempPositronEntry.bunchNumber));
  // singlesTree -> SetBranchAddress("inFillGain", &(tempPositronEntry.inFillGain));
  // singlesTree -> SetBranchAddress("crystalEnergy", &(tempPositronEntry.crystalEnergy));

  // std::cout << "[Debug] before the singlesTree loop" << std::endl;
  // loop over the tree
  for (int i = 0; i < singlesTree -> GetEntries(); i++) {
    singlesTree -> GetEntry(i);
    // add a *copy* of the dummy object to the list of positron objects
    positronEntries.push_back(tempPositronEntry);
  }

  // ===============================================================================================

  // create dummy pileup data object to hold data from current TTree entry
  PileupData tempDoubleEntry;

  // temporary pointers-to-vectors to use for SetBranchAddress
  // these vector contents must be be copied into each pileup data object for each entry
  std::vector<int>* tempPileupIndex = 0;
  std::vector<bool>* tempPileupFlagged = 0;
  std::vector<double>* tempPileupTime = 0;
  std::vector<double>* tempPileupEnergy = 0;
  std::vector<double>* tempPileupX = 0;
  std::vector<double>* tempPileupY = 0;
  std::vector<int>* tempPileupCaloIndex = 0;

  // point the doubles TTree branches to the member variables in the PileupData object for doubles
  // vector types must point to the pointers-to-vectors above
  // non-vector types can point directly inside the dummy object, and will be copied
  doublesTree -> SetBranchAddress("pileupIndex", &tempPileupIndex);
  doublesTree -> SetBranchAddress("pileupFlagged", &tempPileupFlagged);
  doublesTree -> SetBranchAddress("pileupTime", &tempPileupTime);
  doublesTree -> SetBranchAddress("pileupEnergy", &tempPileupEnergy);
  doublesTree -> SetBranchAddress("pileupX", &tempPileupX);
  doublesTree -> SetBranchAddress("pileupY", &tempPileupY);
  doublesTree -> SetBranchAddress("pileupCaloIndex", &tempPileupCaloIndex);
  // doublesTree -> SetBranchAddress("laserInFill", &(tempDoubleEntry.laserInFill));
  doublesTree -> SetBranchAddress("runIndex", &(tempDoubleEntry.runIndex));
  doublesTree -> SetBranchAddress("subrunIndex", &(tempDoubleEntry.subrunIndex));
  doublesTree -> SetBranchAddress("fillIndex", &(tempDoubleEntry.fillIndex));
  doublesTree -> SetBranchAddress("bunchNumber", &(tempDoubleEntry.bunchNumber));

  // std::cout << "[Debug] before the doublesTree loop" << std::endl;
  // loop over the tree
  for (int i = 0; i < doublesTree -> GetEntries(); i++) {
    doublesTree -> GetEntry(i);
    // add a *copy* of the dummy object to the list of double-pileup objects
    doubleEntries.push_back(tempDoubleEntry);
    // must explicitly copy temporary pointers-to-vectors into data object's vectors
    // because ROOT will overwrite its internal buffer that pointers-to-vectors point to
    doubleEntries.back().pileupIndex = *tempPileupIndex;
    doubleEntries.back().pileupFlagged = *tempPileupFlagged;
    doubleEntries.back().pileupTime = *tempPileupTime;
    doubleEntries.back().pileupEnergy = *tempPileupEnergy;
    doubleEntries.back().pileupX = *tempPileupX;
    doubleEntries.back().pileupY = *tempPileupY;
    doubleEntries.back().pileupCaloIndex = *tempPileupCaloIndex;
  }
  // ===============================================================================================

  // create dummy pileup data object to hold data from current TTree entry
  PileupData tempTripleEntry;

  // point the triples TTree branches to the member variables in the PileupData object for triples
  // vector types must point to the pointers-to-vectors declared above for double pileup (reused here)
  // non-vector types can point directly inside the dummy object, and will be copied
  triplesTree -> SetBranchAddress("pileupIndex", &tempPileupIndex);
  triplesTree -> SetBranchAddress("pileupFlagged", &tempPileupFlagged);
  triplesTree -> SetBranchAddress("pileupTime", &tempPileupTime);
  triplesTree -> SetBranchAddress("pileupEnergy", &tempPileupEnergy);
  triplesTree -> SetBranchAddress("pileupX", &tempPileupX);
  triplesTree -> SetBranchAddress("pileupY", &tempPileupY);
  triplesTree -> SetBranchAddress("pileupCaloIndex", &tempPileupCaloIndex);
  // triplesTree -> SetBranchAddress("laserInFill", &(tempTripleEntry.laserInFill));
  triplesTree -> SetBranchAddress("runIndex", &(tempTripleEntry.runIndex));
  triplesTree -> SetBranchAddress("subrunIndex", &(tempTripleEntry.subrunIndex));
  triplesTree -> SetBranchAddress("fillIndex", &(tempTripleEntry.fillIndex));
  triplesTree -> SetBranchAddress("bunchNumber", &(tempTripleEntry.bunchNumber));

  // std::cout << "[Debug] before the triplesTree loop" << std::endl;
  // loop over the tree
  for (int i = 0; i < triplesTree -> GetEntries(); i++) {
    triplesTree -> GetEntry(i);
    // add a *copy* of the dummy object to the list of triple-pileup objects
    tripleEntries.push_back(tempTripleEntry);
    // must explicitly copy temporary pointers-to-vectors into data object's vectors
    // because ROOT will overwrite its internal buffer that pointers-to-vectors point to
    tripleEntries.back().pileupIndex = *tempPileupIndex;
    tripleEntries.back().pileupFlagged = *tempPileupFlagged;
    tripleEntries.back().pileupTime = *tempPileupTime;
    tripleEntries.back().pileupEnergy = *tempPileupEnergy;
    tripleEntries.back().pileupX = *tempPileupX;
    tripleEntries.back().pileupY = *tempPileupY;
    tripleEntries.back().pileupCaloIndex = *tempPileupCaloIndex;
  }

  // ===============================================================================================

  // LostMuonData tempLostMuonEntry;

  // temporary pointers-to-vectors to use for SetBranchAddress
  // these vector contents must be be copied into each lost muon data object for each entry
  std::vector<double>* tempCalo2times = 0;
  std::vector<double>* tempCalo3times = 0;
  std::vector<double>* tempCalo4times = 0;
  std::vector<double>* tempCalo2energies = 0;
  std::vector<double>* tempCalo3energies = 0;
  std::vector<double>* tempCalo4energies = 0;
  std::vector<double>* tempCalo2x = 0;
  std::vector<double>* tempCalo3x = 0;
  std::vector<double>* tempCalo4x = 0;
  std::vector<double>* tempCalo2y = 0;
  std::vector<double>* tempCalo3y = 0;
  std::vector<double>* tempCalo4y = 0;

  // // point the lost muon TTree branches to the member variables in the LostMuonData object for lost muon candidates
  // // vector types must point to the pointers-to-vectors above
  // // non-vector types can point directly inside the dummy object, and will be copied
  // lostMuonTree -> SetBranchAddress("calo1", &(tempLostMuonEntry.calo1));
  // lostMuonTree -> SetBranchAddress("time1", &(tempLostMuonEntry.time1));
  // lostMuonTree -> SetBranchAddress("energy1", &(tempLostMuonEntry.energy1));
  // lostMuonTree -> SetBranchAddress("x1", &(tempLostMuonEntry.x1));
  // lostMuonTree -> SetBranchAddress("y1", &(tempLostMuonEntry.y1));
  // lostMuonTree -> SetBranchAddress("calo2times", &tempCalo2times);
  // lostMuonTree -> SetBranchAddress("calo3times", &tempCalo3times);
  // lostMuonTree -> SetBranchAddress("calo4times", &tempCalo4times);
  // lostMuonTree -> SetBranchAddress("calo2energies", &tempCalo2energies);
  // lostMuonTree -> SetBranchAddress("calo3energies", &tempCalo3energies);
  // lostMuonTree -> SetBranchAddress("calo4energies", &tempCalo4energies);
  // lostMuonTree -> SetBranchAddress("calo2x", &tempCalo2x);
  // lostMuonTree -> SetBranchAddress("calo3x", &tempCalo3x);
  // lostMuonTree -> SetBranchAddress("calo4x", &tempCalo4x);
  // lostMuonTree -> SetBranchAddress("calo2y", &tempCalo2y);
  // lostMuonTree -> SetBranchAddress("calo3y", &tempCalo3y);
  // lostMuonTree -> SetBranchAddress("calo4y", &tempCalo4y);
  // lostMuonTree -> SetBranchAddress("fillIndex", &(tempLostMuonEntry.fillIndex));
  // lostMuonTree -> SetBranchAddress("subrunIndex", &(tempLostMuonEntry.subrunIndex));
  // lostMuonTree -> SetBranchAddress("runIndex", &(tempLostMuonEntry.runIndex));
  // lostMuonTree -> SetBranchAddress("bunchNumber", &(tempLostMuonEntry.bunchNumber));

  // // std::cout << "[Debug] before the lostMuonTree loop" << std::endl;
  // // loop over the tree
  // for (int i = 0; i < lostMuonTree -> GetEntries(); i++) {
  //   lostMuonTree -> GetEntry(i);
  //   // add a *copy* of the dummy object to the list of lost muon candidate objects
  //   lostMuonEntries.push_back(tempLostMuonEntry);
  //   // must explicitly copy temporary pointers-to-vectors into data object's vectors
  //   // because ROOT will overwrite its internal buffer that pointers-to-vectors point to
  //   lostMuonEntries.back().calo2times = *tempCalo2times;
  //   lostMuonEntries.back().calo3times = *tempCalo3times;
  //   lostMuonEntries.back().calo4times = *tempCalo4times;
  //   lostMuonEntries.back().calo2energies = *tempCalo2energies;
  //   lostMuonEntries.back().calo3energies = *tempCalo3energies;
  //   lostMuonEntries.back().calo4energies = *tempCalo4energies;
  //   lostMuonEntries.back().calo2x = *tempCalo2x;
  //   lostMuonEntries.back().calo3x = *tempCalo3x;
  //   lostMuonEntries.back().calo4x = *tempCalo4x;
  //   lostMuonEntries.back().calo2y = *tempCalo2y;
  //   lostMuonEntries.back().calo3y = *tempCalo3y;
  //   lostMuonEntries.back().calo4y = *tempCalo4y;
  // }
  
  // // read the lost muon input data (i.e. expected times of flight and positron rates vs. time, per calorimeter)
  // LostMuonInput lostMuonInput;
  // lostMuonInput.timeOfFlight = (TH1D*) lostMuonFile -> Get("timeOfFlight");
  // lostMuonInput.events = ((TH1D*) lostMuonFile -> Get("events")) -> GetBinContent(1);
  // for (int caloIndex = 0; caloIndex < 24; caloIndex++) {
  //   lostMuonInput.caloEfficiency.push_back((TH1D*) lostMuonFile -> Get(Form("caloeff%d", caloIndex + 1)));
  // }

  // ===============================================================================================

  // initialize one output file for each subclass
  std::vector<TFile*> outputFiles;
  for (std::string& className: classNames) {
    outputFiles.push_back(new TFile(Form("%s/%s_dataset%s_skim%05d.root", outputPath.c_str(), className.c_str(), dataset.c_str(), skimIndex), "RECREATE"));
  }

  // initialize instances of each subclass and book their histograms
  std::vector<HistogramBase*> classInstances;
  for (std::string& className: classNames) {
    // if (className == "CornellHistograms") {
    //   classInstances.push_back(new CornellHistograms());
    // } else if (className == "RatioHistograms") {
      // classInstances.push_back(new RatioHistograms());
    // } else if (className == "Byu2Histograms") {
      // classInstances.push_back(new Byu2Histograms());
    // }
    if (className == "Byu2Histograms") {
      classInstances.push_back(new Byu2Histograms());
    }
    classInstances.back() -> bookHistograms(0, skimIndex);
  }

  // keep track of which fills were marked as in-fill laser fills, since lost muon tree is missing this information
  std::set<long long> laserFillIndices;

  // std::cout << "[Debug] before the random seed loop" << std::endl;
  // loop over random seeds
  for (int seedIndex = 0; seedIndex < 1; seedIndex++) {
    // fprintf(stderr, "Creating histograms for seedIndex = %i\n", seedIndex);

    // create maps from unique fill index to fast rotation and vertical waist randomization amounts
    std::map<long long, double> frRandomizationPerFill;
    std::map<long long, double> vwRandomizationPerFill;

    // create random number generator with unique seed for this skim file + seed index combination
    TRandom3* generator = new TRandom3(seedOffset + seedIndex);

    // keep track of the last uniqueFillIndex so that we don't check if randomization map contains each positron's unique fill index
    // this will save time when we're iterating through a sequence of positrons from the same fill, for example
    long long lastUniqueFillIndex = -1;
    
    // std::cout << "Loop over singles" << std::endl;
    // loop over the preloaded positron entries
    for (int i = 0; i < positronEntries.size(); i++) {
      
      PositronData& positronEntry = positronEntries[i];

      // literals need 'LL' to avoid overflows from intermediate types that are too small
      long long uniqueFillIndex = getUniqueFillIndex(positronEntry.runIndex, positronEntry.subrunIndex, positronEntry.fillIndex);

      // if (positronEntry.laserInFill) { // for now, skip entries from laser-fills
      //   laserFillIndices.insert(uniqueFillIndex);
      //   continue;
      // }

      // if the fill index has changed...
      if (lastUniqueFillIndex != uniqueFillIndex){
        // ...and it's definitely not already in the map (in case entries are out of order)...
        if (frRandomizationPerFill.count(uniqueFillIndex) == 0) {
          // ...add new randomization amounts to the maps for this fill
          frRandomizationPerFill[uniqueFillIndex] = ((generator -> Rndm()) - 0.5) * frPeriod;
          vwRandomizationPerFill[uniqueFillIndex] = ((generator -> Rndm()) - 0.5) * vwPeriod;
        }
      }

      // leave randomization amounts at zero for seedIndex == -1 (unrandomized)
      double frRandomization = 0.0;
      double vwRandomization = 0.0;

      // update randomization amounts from map when seedIndex > -1
      if (seedIndex > -1) {
        frRandomization = frRandomizationPerFill[uniqueFillIndex];
        vwRandomization = vwRandomizationPerFill[uniqueFillIndex];
      }

      for (HistogramBase* instance: classInstances) {
        instance -> fillSinglesHistograms(positronEntry, frRandomization, vwRandomization, seedIndex, skimIndex);
      }

    }

    // std::cout << "Loop over doubles" << std::endl;
    // loop over the preloaded double-pileup entries
    for (int i = 0; i < doubleEntries.size(); i++) {
      
      PileupData& doubleEntry = doubleEntries[i];
      // if (doubleEntry.laserInFill) { // for now, skip entries from laser-fills
      //   continue;
      // }

      // literals need 'LL' to avoid overflows from intermediate types that are too small
      long long uniqueFillIndex = getUniqueFillIndex(doubleEntry.runIndex, doubleEntry.subrunIndex, doubleEntry.fillIndex);

      double frRandomization = 0.0;
      double vwRandomization = 0.0;

      // seedIndex -1 is unrandomized; set randomizationTime to 0
      if (seedIndex > -1) {
        frRandomization = frRandomizationPerFill[uniqueFillIndex];
        vwRandomization = vwRandomizationPerFill[uniqueFillIndex];
      }

      for (HistogramBase* instance: classInstances) {
        instance -> fillDoublesHistograms(doubleEntry, frRandomization, vwRandomization, seedIndex, skimIndex);
      }

    }

    // std::cout << "Loop over triples" << std::endl;
    // loop over the preloaded triple-pileup entries
    for (int i = 0; i < tripleEntries.size(); i++) {
      
      PileupData& tripleEntry = tripleEntries[i];
      // if (tripleEntry.laserInFill) { // for now, skip entries from laser-fills
      //   continue;
      // }

      // literals need 'LL' to avoid overflows from intermediate types that are too small
      long long uniqueFillIndex = getUniqueFillIndex(tripleEntry.runIndex, tripleEntry.subrunIndex, tripleEntry.fillIndex);

      double frRandomization = 0.0;
      double vwRandomization = 0.0;

      // seedIndex -1 is unrandomized; set randomizationTime to 0
      if (seedIndex > -1) {
        frRandomization = frRandomizationPerFill[uniqueFillIndex];
        vwRandomization = vwRandomizationPerFill[uniqueFillIndex];
      }

      for (HistogramBase* instance: classInstances) {
        instance -> fillTriplesHistograms(tripleEntry, frRandomization, vwRandomization, seedIndex, skimIndex);
      }

    } // end for loop over triplesTree

    // std::cout << "Loop over lost muons" << std::endl;
    // // loop over the preloaded lost muon candidate entries
    // for (int i = 0; i < lostMuonEntries.size(); i++) {

    //   LostMuonData& lostMuonEntry = lostMuonEntries[i];

    //   // literals need 'LL' to avoid overflows from intermediate types that are too small
    //   long long uniqueFillIndex = getUniqueFillIndex(lostMuonEntry.runIndex, lostMuonEntry.subrunIndex, lostMuonEntry.fillIndex);

    //   if (laserFillIndices.find(uniqueFillIndex) != laserFillIndices.end()) { // for now, skip entries from laser-fills
    //     continue;
    //   }

    //   double frRandomization = 0.0;
    //   double vwRandomization = 0.0;

    //   // seedIndex -1 is unrandomized; set randomizationTime to 0
    //   if (seedIndex > -1) {
    //     frRandomization = frRandomizationPerFill[uniqueFillIndex];
    //     vwRandomization = vwRandomizationPerFill[uniqueFillIndex];
    //   }

    //   // for (HistogramBase* instance: classInstances) {
    //   //   instance -> fillLostMuonHistograms(lostMuonEntry, lostMuonInput, frRandomization, vwRandomization, seedIndex, skimIndex);
    //   // }

    // }

    // std::cout << "Loop over classes" << std::endl;
    // loop over the instances and write histograms to disk for this seed
    for (unsigned int instanceIndex = 0; instanceIndex < classInstances.size(); instanceIndex++) {
      std::string seedLabel = Form("seed%d", seedIndex);
      outputFiles[instanceIndex] -> mkdir(seedLabel.c_str());
      outputFiles[instanceIndex] -> cd(seedLabel.c_str());
      classInstances[instanceIndex] -> writeHistograms(outputFiles[instanceIndex], seedIndex);
    }

  } // end loop over seedIndex

  // close input and output files
  skimFile -> Close();
  // lostMuonFile -> Close();
  for (TFile* outputFile: outputFiles) {
    outputFile -> Close();
  }

}

