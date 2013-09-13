#include <Track.h>
#include <GFRaveVertex.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>




int main() {


  genfit::Track* aTrackPtr(NULL);
  genfit::GFRaveVertex* aVertexPtr(NULL);

  TFile* trackFile = TFile::Open("tracks.root", "READ");
  TTree* trackTree = (TTree*)trackFile->Get("trackTree");
  TClonesArray* trackArray = new TClonesArray("genfit::Track");
  trackTree->SetBranchAddress("trackBranch", &trackArray);

  trackTree->Print();

  //TFile* vertexFile = TFile::Open("vertices.root", "READ");
  //TTree* vertexTree = (TTree*)trackFile->Get("vertexTree");
  TClonesArray* vertexArray = new TClonesArray("genfit::GFRaveVertex");
  trackTree->SetBranchAddress("vertexBranch", &vertexArray);

  //vertexTree->Print();

  for (Long_t i = 0; i < trackTree->GetEntries(); ++i) {
    trackTree->GetEntry(i);

    std::cout << "trackArray nr of entries: " << trackArray->GetEntries() << "\n";

    for (Long_t j = 0; j < trackArray->GetEntriesFast(); ++j) {
      std::cout << "track uniqueID: " << static_cast<genfit::Track*>(trackArray->At(j))->GetUniqueID() <<
          " (" << static_cast<genfit::Track*>(trackArray->At(j))->GetUniqueID() - 16777216 << ")\n";
    }

    /*for (Long_t j = 0; j < trackArray->GetEntriesFast(); ++j) {
      aTrackPtr = static_cast<genfit::Track*>(trackArray->At(j));
    }*/
  /*}

  for (Long_t i = 0; i < vertexTree->GetEntries(); ++i) {*/
    //vertexTree->GetEntry(i);

    for (Long_t j = 0; j < vertexArray->GetEntriesFast(); ++j) {


      aVertexPtr = (genfit::GFRaveVertex*)(vertexArray->At(j));
      //aVertexPtr->Print();

      for (unsigned int k=0; k<aVertexPtr->getNTracks(); ++k) {
        std::cout << "track parameters uniqueID: " << aVertexPtr->getParameters(k)->GetUniqueID() << "\n";
      }

      // when the track branch from the tracks.root file is loaded, the TRefs to the tracks
      // in the GFRaveTrackParameters are again pointing to them.
      for (unsigned int k = 0; k<aVertexPtr->getNTracks(); ++k) {
        if (aVertexPtr->getParameters(i)->hasTrack()) {
          std::cout << "track parameters have track \n";
        }
        else {
          std::cout << "track parameters have NO track <--------------------------------- \n";
        }
      }


    }

  }


}


