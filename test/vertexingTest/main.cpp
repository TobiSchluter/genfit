#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <GFRaveVertexFactory.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>




int main() {

  gRandom->SetSeed(14);

  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 15.)); // 15 kGauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  // init vertex factory
  genfit::GFRaveVertexFactory vertexFactory;
  vertexFactory.setMethod("kalman-smoothing:1");


  std::vector<genfit::Track*> tracks;

  // main loop
  for (unsigned int iEvent=0; iEvent<100; ++iEvent){

    // clean up
    for (unsigned int i=0; i<tracks.size(); ++i){
      delete tracks[i];
    }
    tracks.clear();

    unsigned int nTracks = gRandom->Uniform(2, 10);
    tracks.reserve(nTracks);

    for (unsigned int iTrack=0; iTrack<nTracks; ++iTrack){

      // true start values
      TVector3 pos(0, 0, 0);
      TVector3 mom(1.,0,0);
      mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
      mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
      mom.SetMag(gRandom->Uniform(0.2, 1.));


      // helix track model
      const int pdg = 13;               // particle pdg code
      const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
      genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
      measurementCreator.setTrackModel(helix);


      unsigned int nMeasurements = gRandom->Uniform(5, 15);


      // smeared start values
      const bool smearPosMom = true;     // init the Reps with smeared pos and mom
      const double posSmear = 0.1;     // cm
      const double momSmear = 3. /180.*TMath::Pi();     // rad
      const double momMagSmear = 0.1;   // relative

      TVector3 posM(pos);
      TVector3 momM(mom);
      if (smearPosMom) {
        posM.SetX(gRandom->Gaus(posM.X(),posSmear));
        posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
        posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

        momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
        momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
        momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
      }
      // approximate covariance
      TMatrixDSym covM(6);
      double resolution = 0.01;
      for (int i = 0; i < 3; ++i)
        covM(i,i) = resolution*resolution;
      for (int i = 3; i < 6; ++i)
        covM(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);


      // trackrep
      genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

      // smeared start state
      genfit::MeasuredStateOnPlane stateSmeared(rep);
      rep->setPosMomCov(stateSmeared, posM, momM, covM);


      // create track
      TVectorD seedState(6);
      TMatrixDSym seedCov(6);
      rep->get6DStateCov(stateSmeared, seedState, seedCov);
      tracks.push_back(new genfit::Track(rep, seedState, seedCov));


      // create random measurement types
      std::vector<genfit::eMeasurementType> measurementTypes;
      for (unsigned int i = 0; i < nMeasurements; ++i)
        measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(7)));


      // create smeared measurements and add to track
      try{
        for (unsigned int i=1; i<measurementTypes.size(); ++i){
          genfit::AbsMeasurement* measurement = measurementCreator.create(measurementTypes[i], i*5.);
          tracks.back()->insertPoint(new genfit::TrackPoint(measurement, tracks.back()));
        }
      }
      catch(genfit::Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        std::cerr << e.what();
        continue; // here is a memleak!
      }

      //check
      assert(tracks.back()->checkConsistency());

      // do the fit
      try{
        fitter->processTrack(tracks.back());
      }
      catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception, next track" << std::endl;
        continue;
      }

      //check
      assert(tracks.back()->checkConsistency());

    } // end loop over tracks

    // vertexing
    std::vector <  genfit::GFRaveVertex* > vertices;
    vertexFactory.findVertices(&vertices, tracks);





    if (iEvent < 1000) {
      // add tracks to event display
      display->addEvent(tracks);
    }

  }// end loop over events

  delete fitter;

  // open event display
  display->open();

}


