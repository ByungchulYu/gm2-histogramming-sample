#include "Byu2Histograms.hh"


template <typename T> bool isElementOf(std::vector<T> array, T element)
{
    for (uint i=0; i<array.size(); i++) {
        if (element == array[i]) {
            return true;
        }
    }
    return false;
}

// Constructor.
Byu2Histograms::Byu2Histograms()
{
    // Initialization of the energy time bin information
	t_min                = 0; 								// in us
	t_max                = 700; 							// in us
	t_n_bins             = (t_max - t_min) / 0.1492;  		// 4691
	t_max                = t_min + t_n_bins * 0.1492; 		// 699.8972 us

	E_min                = 1050;							// in MeV
	E_max                = 3060;							// in MeV
	E_bin_width          = 67;								// in MeV
	E_n_bins             = (E_max - E_min) / E_bin_width;   // 30
	E_max                = E_min + E_n_bins * E_bin_width;  // 3060 MeV

    // Initialization of the runIndex and subrunIndex
	prev_runIndexS_	    =-1;
	prev_runIndexD_	    =-1;
	prev_runIndexH_	    =-1;

    prev_subrunIndexS_ = -1;
    prev_subrunIndexD_ = -1;
    prev_subrunIndexH_ = -1;

    
    // Initialization of the tree and branch
    TREE_ET_aux_            = new TTree("ET", "ET");    

    prev_runindex_S     = TREE_ET_aux_->Branch("prev_runIndexS_", &prev_runIndexS_, "prev_runIndexS_/I");
    prev_runindex_D     = TREE_ET_aux_->Branch("prev_runIndexD_", &prev_runIndexD_, "prev_runIndexD_/I");
    prev_runindex_H     = TREE_ET_aux_->Branch("prev_runIndexH_", &prev_runIndexH_, "prev_runIndexH_/I");

    prev_index_S        = TREE_ET_aux_->Branch("prev_subrunIndexS_", &prev_subrunIndexS_, "prev_subrunIndexS_/I");
    prev_index_D        = TREE_ET_aux_->Branch("prev_subrunIndexD_", &prev_subrunIndexD_, "prev_subrunIndexD_/I");
    prev_index_H        = TREE_ET_aux_->Branch("prev_subrunIndexH_", &prev_subrunIndexH_, "prev_subrunIndexH_/I");
    
    // initialization of subruntime
    subruntimeindex_    = 0;
    subruntime_         = TREE_ET_aux_->Branch("subruntimeindex_", &subruntimeindex_, "subruntimeindex_/D");
}


void Byu2Histograms::bookHistograms(int seedIndex, int skimIndex)
{

    // Histograms initialization 
	EvsT_                = new TH2F("EvsT_",    "Energy vs Time             ;Time [us]; Energy [MeV]", t_n_bins, t_min, t_max, E_n_bins, E_min, E_max);
	EvsT_D_	             = new TH2F("EvsT_D_",  "Energy vs Time (Double PU) ;Time [us]; Energy [MeV]", t_n_bins, t_min, t_max, E_n_bins, E_min, E_max);
	EvsT_H_	             = new TH2F("EvsT_H_",  "Energy vs Time (Higher PU) ;Time [us]; Energy [MeV]", t_n_bins, t_min, t_max, E_n_bins, E_min, E_max);
	EvsT_PU_	         = new TH2F("EvsT_PU_", "Energy vs Time (Total PU)  ;Time [us]; Energy [MeV]", t_n_bins, t_min, t_max, E_n_bins, E_min, E_max);

    // Branch initialization (This must be after the histogram initialization)
    EvsT_branch     = TREE_ET_aux_->Branch("EvsT_",    "TH2F",  &EvsT_);
    EvsT_D_branch   = TREE_ET_aux_->Branch("EvsT_D_",  "TH2F",  &EvsT_D_);
    EvsT_H_branch   = TREE_ET_aux_->Branch("EvsT_H_",  "TH2F",  &EvsT_H_);
    EvsT_PU_branch  = TREE_ET_aux_->Branch("EvsT_PU_", "TH2F",  &EvsT_PU_);

}




void Byu2Histograms::fillSinglesHistograms(PositronData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex)
{

    
    // Call energy and time for each entry
    double energy = entry.energy;
    double convertedTime = entry.time * ct2us + frRandomization;
    int caloIndex = entry.caloIndex;
    timestamps_.push_back(entry.gpsInteger);
    // Exclude calorimeter 18 in Run2F dataset       //////////////////////////Don't Forget////////////////////////
    // if (caloIndex == 18) {
    //     return;
    // }   

    // Fill EvsT histogram, runIndex, and subrunIndex branches in this if statement 
    if ((entry.runIndex != prev_runIndexS_ || entry.subrunIndex != prev_subrunIndexS_) && prev_subrunIndexS_ != -1) {
        double gpstime_size = timestamps_.size();
        double average_time = 0;
        for (unsigned int i = 0; i < timestamps_.size(); i++){
		    average_time += timestamps_[i] / gpstime_size;
	    }
	    subruntimeindex_ = average_time;
        subruntime_->Fill();
        timestamps_.clear();
        
        EvsT_->SetTitle(Form("EvsT_subrun%d", entry.subrunIndex));
        prev_index_S->Fill();
        prev_runindex_S->Fill();
        EvsT_branch->Fill();
        EvsT_->Reset();
    }

    // Fill clusters (entry in PositronData) in EvsT histogram (assign runIndex and subrunIndex after each cluster is filled)
    EvsT_->Fill(convertedTime, energy);
    prev_subrunIndexS_      = entry.subrunIndex; // Update previousSubrunIndex
    prev_runIndexS_         = entry.runIndex;
    subruntimeindex_        = entry.gpsInteger;
}




void Byu2Histograms::fillDoublesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex)
{

    // Fill EvsT histogram, runIndex, and subrunIndex branches in this if statement 
    if ((entry.runIndex != prev_runIndexD_ || entry.subrunIndex != prev_subrunIndexD_) && prev_subrunIndexD_ != -1) {
        EvsT_D_->SetTitle(Form("EvsT_D_subrun%d", entry.subrunIndex));
        prev_index_D->Fill();
        prev_runindex_D->Fill();
        EvsT_D_branch->Fill();
        EvsT_D_->Reset();
    }

    // Fill clusters (entry in PileupData) in EvsT_D histogram with proper weights (assign runIndex and subrunIndex after each cluster is filled)
    double weight = 0;
    for (uint i=0; i<entry.pileupIndex.size(); i++) {

        double energy = entry.pileupEnergy.at(i);
        double convertedTime = entry.pileupTime.at(i) * ct2us + frRandomization + 0.5 * cyclotronPeriod;
        int caloIndex = entry.pileupCaloIndex.at(i);

        // Exclude calorimeter 18 in Run2F dataset       //////////////////////////Don't Forget////////////////////////
        // if (caloIndex == 18) {
        //     return;
        // }    


        if (entry.pileupIndex.at(i) == 2) {
            weight = 0.5;
            EvsT_D_->Fill(convertedTime, energy, weight);

        } 
        else {
            weight = -0.5;
            EvsT_D_->Fill(convertedTime, energy, weight);

        }

    } 
    prev_subrunIndexD_      = entry.subrunIndex; // Update previousSubrunIndex
    prev_runIndexD_         = entry.runIndex;
        

}



void Byu2Histograms::fillTriplesHistograms(PileupData& entry, double frRandomization, double vwRandomization, int seedIndex, int skimIndex)
{
    
    // Fill EvsT histogram, runIndex, and subrunIndex branches in this if statement 
    if ((entry.runIndex != prev_runIndexH_ || entry.subrunIndex != prev_subrunIndexH_) && prev_subrunIndexH_ != -1) {
        EvsT_H_->SetTitle(Form("EvsT_H_subrun%d", entry.subrunIndex));
        prev_index_H->Fill();
        prev_runindex_H->Fill();
        EvsT_H_branch->Fill();
        EvsT_H_->Reset();
    }


    // Fill clusters (entry in PileupData) in EvsT_H histogram with proper weights (assign runIndex and subrunIndex after each cluster is filled)
    std::vector<int> pu3DoublesIndices = {0, 1, 2, 6, 9, 12};     // Analogous to doubles in double-pileup events.
    std::vector<int> pu3SinglesIndices = {3, 4, 5, 7, 8, 10, 11}; // Analogous to singles in double-pileup events.
    double weight = 0;
    for (uint i=0; i<entry.pileupIndex.size(); i++) {
        double energy = entry.pileupEnergy.at(i);
        double convertedTime = entry.pileupTime.at(i) * ct2us + frRandomization + cyclotronPeriod;
        int caloIndex = entry.pileupCaloIndex.at(i);

        // Exclude calorimeter 18 in Run2F dataset       //////////////////////////Don't Forget////////////////////////
        // if (caloIndex == 18) {
        //     return;
        // }   
        
        if (isElementOf(pu3DoublesIndices, entry.pileupIndex.at(i))) {
            weight = 0.5;
            EvsT_H_->Fill(convertedTime, energy, weight);

        } else if (isElementOf(pu3SinglesIndices, entry.pileupIndex.at(i))) {
            weight = -0.5;
            EvsT_H_->Fill(convertedTime, energy, weight);

        } else {
            continue;
        }

    }
    prev_subrunIndexH_ = entry.subrunIndex; // Update previousSubrunIndex
    prev_runIndexH_         = entry.runIndex;

}

void Byu2Histograms::fillLostMuonHistograms(LostMuonData& entry, LostMuonInput& lmInput, double frRandomization, double vwRandomization, int seedIndex, int skimIndex)
{
    
}

void Byu2Histograms::writeHistograms(TFile* outputFile, int seedIndex)
{   
    
    // 각 fillsingles/double/triple Histograms 함수들에서 마지막으로 들어온 subrun에 대해서는 각 함수에 있는 if문의 조건이 성립되지 않아서 EvsT_ EvsT_D_ EvsT_H_가 tree에 fill이 안되었다.
    // 이곳 writeHistgrams에서 저 히스토그램들을 각각의 branch에 fill을 해준다.
    prev_runindex_S->Fill();
    prev_index_S->Fill();
    EvsT_branch->Fill();
    
    prev_runindex_D->Fill();
    prev_index_D->Fill();
    EvsT_D_branch->Fill();

    prev_runindex_H->Fill();
    prev_index_H->Fill();
    EvsT_H_branch->Fill();

    double gpstime_size = timestamps_.size();
    double average_time = 0;
    for (unsigned int i = 0; i < timestamps_.size(); i++){
        average_time += timestamps_[i] / gpstime_size;
    }
    subruntimeindex_ = average_time;
    subruntime_->Fill();

    // tree에 fill을 하지 않고 branch마다 fill을 따로 하였기 때문에 tree의 entry는 수동으로 아래와 같이 직접 정해주어야 한다.
    TREE_ET_aux_->SetEntries(prev_index_S->GetEntries());
    // TREE_ET_aux_->SetEntries(10);
    
   
    // Fill PU EvsT histogram
    TH2F* double_h = new TH2F();
    TREE_ET_aux_->SetBranchAddress("EvsT_D_", &double_h);
    TH2F* higher_h = new TH2F();
    TREE_ET_aux_->SetBranchAddress("EvsT_H_", &higher_h);
    for (int i = 0; i < TREE_ET_aux_->GetEntries(); i++){
        TREE_ET_aux_->GetEntry(i);
        EvsT_PU_->Reset();
        EvsT_PU_->Add(double_h, 1);
        EvsT_PU_->Add(higher_h, 1);
        EvsT_PU_branch->Fill();
    }

    


    // TREE_subrun_->Fill();
    // TREE_ET_aux_->Write();


    TREE_ET_aux_->SetBranchStatus("*", 0);
    TREE_ET_aux_->SetBranchStatus("prev_runIndexS_", 1);
    TREE_ET_aux_->SetBranchStatus("prev_subrunIndexS_", 1);
    TREE_ET_aux_->SetBranchStatus("subruntimeindex_", 1);
    TREE_ET_aux_->SetBranchStatus("EvsT_", 1);
    TREE_ET_aux_->SetBranchStatus("EvsT_PU_", 1);

    TREE_ET_ = TREE_ET_aux_->CloneTree();

    TREE_ET_->Write();

}

