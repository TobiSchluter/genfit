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

  const unsigned int nEvents = 100;
  const unsigned int nMeasurements = 11;
  const double BField = 15.;       // kGauss
  const double momentum = 0.4;     // GeV

  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = true;     // init the Reps with smeared pos and mom
  const double posSmear = 0.01;     // cm
  const double momSmear = 5. /180.*TMath::Pi();     // rad
  const double momMagSmear = 0.2;   // relative
  const double zSmearFac = 2;



  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();



  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // main loop
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){

    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0,0);
    mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
    mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
    mom.SetMag(momentum);


    // helix track model
    const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    measurementCreator.setTrackModel(helix);

    // smeared start values
    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),zSmearFac*posSmear));

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
    genfit::Track fitTrack(rep, seedState, seedCov);


    // create random measurement types
    std::vector<genfit::eMeasurementType> measurementTypes;
    for (unsigned int i = 0; i < nMeasurements; ++i)
      measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(7)));


    // create smeared measurements and add to track
    try{
      for (unsigned int i=0; i<measurementTypes.size(); ++i){
        genfit::AbsMeasurement* measurement = measurementCreator.create(measurementTypes[i], i*5.);
        fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
      }
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
      continue; // here is a memleak!
    }

    //check
    assert(fitTrack.checkConsistency());

    // do the fit
    try{
      fitter->processTrack(&fitTrack);
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }

    //check
    assert(fitTrack.checkConsistency());


    if (iEvent < 1000) {
      // add track to event display
      display->addEvent(&fitTrack);
    }



  }// end loop over events

  delete fitter;

  // open event display
  display->open();

}


