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
    /*
    fillSinglesHistograms는 positron cluster (entry) 에 대해 for문을 이용하여 EvsT_ 히스토그램에 fill을 해주는 함수이다. (반복문은 HistogramBase.hh에 이미 할당되었다)
    1. 처음에 prev_subrunIndexS_ 이 constructor에서 0으로 initialization이 되어 있기 때문에 if문에 들어가지 않고, 바로 EvsT에 첫 cluster (entry)가 fill이 된다.
    2. 이후 첫번째 prev_subrunIndex에 해당 entry의 subrunIndex가 assign이 된 후 다시 if문으로 가지만, entry.subrunIndex와 prev_subrunIndex가 같기에 if문으로 들어가지지 않는다.
    3. 결국 하나의 subrunIndex에 있는 모든 cluster들의 fill이 다 끝날때까지 if문으로 들어가지지 않고, 그 다음 subrunIndex가 prev_subrunIndex로 assign이 된 후에 if문으로 들어가진다.
    4. 즉, 동일한 subrunIndex를 가지고 있는 모든 entry들에 대해서 EvsT_에 fill을 다하고 다시 다음 subrunIndex가 되었을 때 if문으로 들어가지게 된다.
    5. if문에 들어오면, 하나의 subrunIndex에 대하여 모든 entry들에 대해 fill이 된 EvsT_ 2F 히소토그램이 Tree에 fill이되고 동시에 prev_subrunIndex도 Tree에 fill이 되며, EvsT_은 Reset이 된다. 
    그 다음부터는 위의 logic을 반복하며 동일한 subrunIndex를 갖는 모든 cluster들에 대해 EvsT_ 2F 히스토그램에 fill이 되고 이런 EvsT_히스토그램은 subrun의 개수만큼 Tree에 fill이 된다. 
    */
    
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
        /*
        fillDoubleHistograms는 double pileup cluster (entry) 에 대해 for 문을 이용하여 EvsT_D_ 히스토그램에 fill을 해주는 함수이다. (반복문은 HistogramBase.hh에 이미 할당되었다)
        유의할 점은, fillsinglesHistogram과는 다르게, entry 하나가 positron 한 개 가 아니라 pileup reconstruction의 묶음에 해당된다. 이 묶음은 pilupIndex가 0,1,2로 되어 있다.
        이 묶음은 pilupIndex가 0,1,2로 되어 있다. 0과 1이 아주 짧은 시간 동안에 동시에 검출된 positron에 해당되고, 2는 이들을 하나로 취급했을 때의 double pileup에 해당된다.
        1. 각 entry마다 (0,1,2) 3개의 index에 대한 정보를 갖는 pileup reconstruction 1묶음이 대응된다. (종종 index가 2인 녀석이 1개 이상이라서 0과 1이 각각 n개씩 포함된 묶음도 있다.)
        2. 그러므로 for문을 통해서 각 묶음에 대한 정보를 불러와서 적합한 weight를 주어 EvsT_D_에 fill을 한다. 
        3. 1번과 2번의 과정이 잘 이해가 되었다면, fillsinglesHistogram에서 한 것 처럼 subrunIndex가 바뀔 때 마다 Tree에 Fill을 해주는 if문을 for문 보다 먼저 위치시킨다.
        */


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

