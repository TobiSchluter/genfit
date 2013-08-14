
#include "EventDisplay.h"

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#include <AbsMeasurement.h>
#include <PlanarMeasurement.h>
#include <ProlateSpacepointMeasurement.h>
#include <SpacepointMeasurement.h>
#include <WireMeasurement.h>
#include <WirePointMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <Tools.h>
#include <KalmanFitterInfo.h>

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveStraightLineSet.h>
#include <TEveTriangleSet.h>
#include <TDecompSVD.h>
#include <TGButton.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoSphere.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixDEigen.h>
#include <TROOT.h>
#include <TVector2.h>
#include <TVectorD.h>
#include <TSystem.h>

ClassImp(genfit::EventDisplay)

namespace genfit {


EventDisplay* EventDisplay::eventDisplay_ = NULL;

EventDisplay::EventDisplay() {

  if((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication))) {
    std::cout << "In EventDisplay ctor: gApplication not found, creating..." << std::flush;
    new TApplication("ROOT_application", 0, 0);
    std::cout << "done!" << std::endl;
  }
  if(!gEve) {
    std::cout << "In EventDisplay ctor: gEve not found, creating..." << std::flush;
    TEveManager::Create();
    std::cout << "done!" << std::endl;
  }

  eventId_ = 0;
  setOptions();
  setErrScale();

}

void EventDisplay::setOptions(std::string opts) { option_ = opts; }

void EventDisplay::setErrScale(double errScale) { errorScale_ = errScale; }

double EventDisplay::getErrScale() { return errorScale_; }

EventDisplay* EventDisplay::getInstance() {

  if(eventDisplay_ == NULL) {
    eventDisplay_ = new EventDisplay();
  }
  return eventDisplay_;

}

EventDisplay::~EventDisplay() { reset(); }

void EventDisplay::reset() {

  for(unsigned int i = 0; i < events_.size(); i++) {

    for(unsigned int j = 0; j < events_[i]->size(); j++) {

      delete events_[i]->at(j);

    }
    delete events_[i];
  }

  events_.clear();
}

void EventDisplay::addEvent(std::vector<Track*>& evts) {

  std::vector<Track*>* vec = new std::vector<Track*>;

  for(unsigned int i = 0; i < evts.size(); i++) {

    vec->push_back(new Track(*(evts[i])));

  }

  events_.push_back(vec);

}

void EventDisplay::next(unsigned int stp) {

  gotoEvent(eventId_ + stp);

}

void EventDisplay::prev(unsigned int stp) {

  if(events_.size() == 0) return;
  if(eventId_ < (int)stp) {
    gotoEvent(0);
  } else {
    gotoEvent(eventId_ - stp);
  }

}

int EventDisplay::getNEvents() { return events_.size(); }

void EventDisplay::gotoEvent(unsigned int id) {

  if(id >= events_.size()) id = events_.size() - 1;

  eventId_ = id;

  std::cout << "At event " << id << std::endl;
  if (gEve->GetCurrentEvent()) {
    gEve->GetCurrentEvent()->DestroyElements();
  }
  double old_error_scale = errorScale_;
  drawEvent(eventId_);
  if(old_error_scale != errorScale_) {
    if (gEve->GetCurrentEvent()) {
      gEve->GetCurrentEvent()->DestroyElements();
    }
    drawEvent(eventId_); // if autoscaling changed the error, draw again.
  }
  errorScale_ = old_error_scale;

}

void EventDisplay::open() {

  std::cout << "EventDisplay::open(); " << events_.size() << " events loaded" << std::endl;

  bool drawSilent = false;
  bool drawGeometry = false;

  // parse the global options
  for(size_t i = 0; i < option_.length(); i++) {
    if(option_[i] == 'X') drawSilent = true;
    if(option_[i] == 'G') drawGeometry = true;
  }

  // draw the geometry, does not really work yet. If it's fixed, the docu in the header file should be changed.
  if(drawGeometry) {
    TGeoNode* top_node = gGeoManager->GetTopNode();
    assert(top_node != NULL);

    //Set transparency & color of geometry
    TObjArray* volumes = gGeoManager->GetListOfVolumes();
    for(int i = 0; i < volumes->GetEntriesFast(); i++) {
      TGeoVolume* volume = dynamic_cast<TGeoVolume*>(volumes->At(i));
      assert(volume != NULL);
      volume->SetLineColor(12);
      volume->SetTransparency(50);
    }

    TEveGeoTopNode* eve_top_node = new TEveGeoTopNode(gGeoManager, top_node);
    eve_top_node->IncDenyDestroy();
    gEve->AddGlobalElement(eve_top_node);
  }

  if(getNEvents() > 0) {
    double old_error_scale = errorScale_;
    drawEvent(0);
    if(old_error_scale != errorScale_) gotoEvent(0); // if autoscaling changed the error, draw again.
    errorScale_ = old_error_scale;
  }


  if(!drawSilent) {
    makeGui();
    gApplication->Run(kTRUE);
  }

  std::cout << "opened" << std::endl;

}

void EventDisplay::drawEvent(unsigned int id) {

  std::cout << "EventDisplay::drawEvent(" << id << ")" << std::endl;

  // parse the option string ------------------------------------------------------------------------
  bool drawAutoScale = false;
  bool drawBackward = false;
  bool drawDetectors = false;
  bool drawErrors = false;
  bool drawForward = false;
  bool drawHits = false;
  bool drawScaleMan = false;
  bool drawTrackMarkers = false;
  bool drawPlanes = false;
  bool drawTrack = false;

  if(option_ != "") {
    for(size_t i = 0; i < option_.length(); i++) {
      if(option_[i] == 'A') drawAutoScale = true;
      if(option_[i] == 'B') drawBackward = true;
      if(option_[i] == 'D') drawDetectors = true;
      if(option_[i] == 'E') drawErrors = true;
      if(option_[i] == 'F') drawForward = true;
      if(option_[i] == 'H') drawHits = true;
      if(option_[i] == 'M') drawTrackMarkers = true;
      if(option_[i] == 'P') drawPlanes = true;
      if(option_[i] == 'S') drawScaleMan = true;
      if(option_[i] == 'T') drawTrack = true;
    }
  }
  // finished parsing the option string -------------------------------------------------------------



  for(unsigned int i = 0; i < events_[id]->size(); i++) { // loop over all tracks in an event

    Track* track = events_[id]->at(i);
    if (! track->checkConsistency()){
      std::cerr<<"track is not consistent"<<std::endl;
      continue;
    }

    AbsTrackRep* rep(track->getCardinalRep());

    //track->Print();
    track->Print("C");
    track->getFitStatus(rep)->Print();

    unsigned int numhits = track->getNumPointsWithMeasurement();

    KalmanFitterInfo* fi;
    KalmanFitterInfo* prevFi = 0;
    const MeasuredStateOnPlane* fittedState(NULL);
    const MeasuredStateOnPlane* prevFittedState(NULL);

    for(unsigned int j = 0; j < numhits; j++) { // loop over all hits in the track

      TrackPoint* tp = track->getPointWithMeasurement(j);
      if (! tp->hasRawMeasurements()) {
        std::cerr<<"trackPoint has no raw measurements"<<std::endl;
        continue;
      }

      if (! tp->hasFitterInfo(rep)) {
        std::cerr<<"trackPoint has no fitterInfo for rep"<<std::endl;
        continue;
      }

      // get the fitter infos ------------------------------------------------------------------
      AbsFitterInfo* fitterInfo = tp->getFitterInfo(rep);
      const AbsMeasurement* m = tp->getRawMeasurement();  // FIXME draw all measurements, not only 1st

      fi = dynamic_cast<KalmanFitterInfo*>(fitterInfo);
      if(fi == NULL) {
        std::cerr<<"can only display KalmanFitterInfo"<<std::endl;
        continue;
      }
      if (! fi->hasPredictionsAndUpdates()) {
        std::cerr<<"KalmanFitterInfo does not have all predictions and updates"<<std::endl;
        continue;
      }
      try {
        fittedState = &(fi->getFittedState(true));
      }
      catch (Exception& e) {
        std::cerr << e.what();
        std::cerr<<"can not get fitted state"<<std::endl;
        continue;
      }

      TVector3 track_pos = rep->getPos(*fittedState);

      double charge = rep->getCharge(*fittedState);

      const MeasurementOnPlane* mop = fi->getMeasurementOnPlane(); // FIXME draw all measurements, not only 1st
      const TVectorT<double>& hit_coords = mop->getState();
      const TMatrixTSym<double>& hit_cov = mop->getCov();

      // finished getting the hit infos -----------------------------------------------------

      // sort hit infos into variables ------------------------------------------------------
      TVector3 o = fittedState->getPlane()->getO();
      TVector3 u = fittedState->getPlane()->getU();
      TVector3 v = fittedState->getPlane()->getV();

      bool planar_hit = false;
      bool planar_pixel_hit = false;
      bool space_hit = false;
      bool wire_hit = false;
      bool wirepoint_hit = false;
      double_t hit_u = 0;
      double_t hit_v = 0;
      double_t plane_size = 4;

      int hit_coords_dim = hit_coords.GetNrows();

      if(dynamic_cast<const PlanarMeasurement*>(m) != NULL) {
        planar_hit = true;
        if(hit_coords_dim == 1) {
          hit_u = hit_coords(0);
        } else if(hit_coords_dim == 2) {
          planar_pixel_hit = true;
          hit_u = hit_coords(0);
          hit_v = hit_coords(1);
        }
      } else if (dynamic_cast<const SpacepointMeasurement*>(m) != NULL) {
        space_hit = true;
        plane_size = 4;
      } else if (dynamic_cast<const WireMeasurement*>(m) != NULL) {
        wire_hit = true;
        hit_u = fabs(hit_coords(0));
        hit_v = v*(track_pos-o); // move the covariance tube so that the track goes through it
        plane_size = 4;
        if (dynamic_cast<const WirePointMeasurement*>(m) != NULL) {
          wirepoint_hit = true;
          hit_v = hit_coords(1);
        }
      } else {
        std::cout << "Track " << i << ", Hit " << j << ": Unknown measurement type: skipping hit!" << std::endl;
        break;
      }

      if(plane_size < 4) plane_size = 4;
      // finished setting variables ---------------------------------------------------------

      // draw planes if corresponding option is set -----------------------------------------
      if(drawPlanes || (drawDetectors && planar_hit)) {
        TVector3 move(0,0,0);
        if (wire_hit) move = v*(v*(track_pos-o)); // move the plane along the wire until the track goes through it
        TEveBox* box = boxCreator(o + move, u, v, plane_size, plane_size, 0.01);
        if (drawDetectors && planar_hit) {
          box->SetMainColor(kCyan);
        } else {
          box->SetMainColor(kGray);
        }
        box->SetMainTransparency(50);
        gEve->AddElement(box);
      }
      // finished drawing planes ------------------------------------------------------------

      // draw track if corresponding option is set ------------------------------------------
      struct makeLinesClass {
      void operator()(const StateOnPlane* prevState, const StateOnPlane* state, const AbsTrackRep* rep,
          const Color_t& color, const Style_t& style, bool drawMarkers, bool drawErrors, double lineWidth = 2, int markerPos = 1)
        {
          TVector3 pos, dir, oldPos, oldDir;
          rep->getPosDir(*state, pos, dir);
          rep->getPosDir(*prevState, oldPos, oldDir);

          double distA = (pos-oldPos).Mag();
          double distB = distA;
          if ((pos-oldPos)*oldDir < 0)
            distA *= -1.;
          if ((pos-oldPos)*dir < 0)
            distB *= -1.;
          TVector3 intermediate1 = oldPos + 0.3 * distA * oldDir;
          TVector3 intermediate2 = pos - 0.3 * distB * dir;
          TEveStraightLineSet* ls = new TEveStraightLineSet;
          ls->AddLine(oldPos(0), oldPos(1), oldPos(2), intermediate1(0), intermediate1(1), intermediate1(2));
          ls->AddLine(intermediate1(0), intermediate1(1), intermediate1(2), intermediate2(0), intermediate2(1), intermediate2(2));
          ls->AddLine(intermediate2(0), intermediate2(1), intermediate2(2), pos(0), pos(1), pos(2));
          ls->SetLineColor(color);
          ls->SetLineStyle(style);
          ls->SetLineWidth(lineWidth);
          if (drawMarkers) {
            if (markerPos == 0)
              ls->AddMarker(oldPos(0), oldPos(1), oldPos(2));
            else
              ls->AddMarker(pos(0), pos(1), pos(2));
          }

          if (lineWidth > 0)
            gEve->AddElement(ls);


          if (drawErrors) {
            const MeasuredStateOnPlane* measuredState;
            if (markerPos == 0)
              measuredState = dynamic_cast<const MeasuredStateOnPlane*>(prevState);
            else
              measuredState = dynamic_cast<const MeasuredStateOnPlane*>(state);

            if (measuredState != NULL) {

              // step for evaluate at a distance from the original plane
              TVector3 eval;
              if (markerPos == 0)
                eval = 0.2 * distA * oldDir;
              else
                eval = -0.2 * distB * dir;


              // get cov at first plane
              TMatrixDSym cov;
              TVector3 position, direction;
              rep->getPosMomCov(*measuredState, position, direction, cov);

              // get eigenvalues & -vectors
              TMatrixDEigen eigen_values(cov.GetSub(0,2, 0,2));
              TMatrixT<double> ev = eigen_values.GetEigenValues();
              TMatrixT<double> eVec = eigen_values.GetEigenVectors();
              TVector3 eVec1, eVec2;
              // limit
              static const double maxErr = 1000.;
              double ev0 = std::min(ev(0,0), maxErr);
              double ev1 = std::min(ev(1,1), maxErr);
              double ev2 = std::min(ev(2,2), maxErr);

              // get two largest eigenvalues/-vectors
              if (ev0 < ev1 && ev0 < ev2) {
                eVec1.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
                eVec1 *= sqrt(ev1);
                eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
                eVec2 *= sqrt(ev2);
              }
              else if (ev1 < ev0 && ev1 < ev2) {
                eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
                eVec1 *= sqrt(ev0);
                eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
                eVec2 *= sqrt(ev2);
              }
              else {
                eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
                eVec1 *= sqrt(ev0);
                eVec2.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
                eVec2 *= sqrt(ev1);
              }

              if (eVec1.Cross(eVec2)*eval < 0)
                eVec2 *= -1;
              //assert(eVec1.Cross(eVec2)*eval > 0);

              const TVector3 oldEVec1(eVec1);
              const TVector3 oldEVec2(eVec2);

              const int nEdges = 24;
              std::vector<TVector3> vertices;

              vertices.push_back(position);

              // vertices at plane
              for (int i=0; i<nEdges; ++i) {
                const double angle = 2*TMath::Pi()/nEdges * i;
                vertices.push_back(position + cos(angle)*eVec1 + sin(angle)*eVec2);
              }



              DetPlane* newPlane = new DetPlane(*(measuredState->getPlane()));
              newPlane->setO(position + eval);

              MeasuredStateOnPlane stateCopy(*measuredState);
              try{
                rep->extrapolateToPlane(stateCopy, SharedPlanePtr(newPlane));
              }
              catch(Exception& e){
                std::cerr<<e.what();
                return;
              }

              // get cov at 2nd plane
              rep->getPosMomCov(stateCopy, position, direction, cov);

              // get eigenvalues & -vectors
              TMatrixDEigen eigen_values2(cov.GetSub(0,2, 0,2));
              ev = eigen_values2.GetEigenValues();
              eVec = eigen_values2.GetEigenVectors();
              // limit
              ev0 = std::min(ev(0,0), maxErr);
              ev1 = std::min(ev(1,1), maxErr);
              ev2 = std::min(ev(2,2), maxErr);

              // get two largest eigenvalues/-vectors
              if (ev0 < ev1 && ev0 < ev2) {
                eVec1.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
                eVec1 *= sqrt(ev1);
                eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
                eVec2 *= sqrt(ev2);
              }
              else if (ev1 < ev0 && ev1 < ev2) {
                eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
                eVec1 *= sqrt(ev0);
                eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
                eVec2 *= sqrt(ev2);
              }
              else {
                eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
                eVec1 *= sqrt(ev0);
                eVec2.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
                eVec2 *= sqrt(ev1);
              }

              if (eVec1.Cross(eVec2)*eval < 0)
                eVec2 *= -1;
              //assert(eVec1.Cross(eVec2)*eval > 0);

              if (oldEVec1*eVec1 < 0) {
                eVec1 *= -1;
                eVec2 *= -1;
              }

              // vertices at 2nd plane
              double angle0 = eVec1.Angle(oldEVec1);
              if (eVec1*(eval.Cross(oldEVec1)) < 0)
                angle0 *= -1;
              for (int i=0; i<nEdges; ++i) {
                const double angle = 2*TMath::Pi()/nEdges * i - angle0;
                vertices.push_back(position + cos(angle)*eVec1 + sin(angle)*eVec2);
              }

              vertices.push_back(position);


              TEveTriangleSet* error_shape = new TEveTriangleSet(vertices.size(), nEdges*2);
              for(unsigned int k = 0; k < vertices.size(); ++k) {
                error_shape->SetVertex(k, vertices[k].X(), vertices[k].Y(), vertices[k].Z());
              }

              assert(vertices.size() == 2*nEdges+2);

              int iTri(0);
              for (int i=0; i<nEdges; ++i) {
                //error_shape->SetTriangle(iTri++,  0,             i+1,        (i+1)%nEdges+1);
                error_shape->SetTriangle(iTri++,  i+1,           i+1+nEdges, (i+1)%nEdges+1);
                error_shape->SetTriangle(iTri++, (i+1)%nEdges+1, i+1+nEdges, (i+1)%nEdges+1+nEdges);
                //error_shape->SetTriangle(iTri++,  2*nEdges+1,    i+1+nEdges, (i+1)%nEdges+1+nEdges);
              }

              //assert(iTri == nEdges*4);

              error_shape->SetMainColor(color);
              error_shape->SetMainTransparency(25);
              gEve->AddElement(error_shape);
            }
          }

        }
      } makeLines;

      if (j > 0) {
        if(drawTrack) {
          makeLines(prevFittedState, fittedState, rep, charge > 0 ? kRed : kBlue, 1, drawTrackMarkers, drawErrors, 3);
          if (drawErrors) { // make sure to draw errors in both directions
            makeLines(prevFittedState, fittedState, rep, charge > 0 ? kRed : kBlue, 1, false, drawErrors, 0, 0);
          }
        }
        if (drawForward)
          makeLines(prevFi->getForwardUpdate(), fi->getForwardPrediction(), rep, charge > 0 ? kMagenta : kCyan, 1, drawTrackMarkers, drawErrors, 1, 0);
        if (drawBackward)
          makeLines(prevFi->getBackwardPrediction(), fi->getBackwardUpdate(), rep, charge > 0 ? kYellow : kMagenta, 1, drawTrackMarkers, drawErrors, 1);
        // draw reference track if corresponding option is set ------------------------------------------
        if(drawTrack && fi->hasReferenceState() && prevFi->hasReferenceState())
          makeLines(prevFi->getReferenceState(), fi->getReferenceState(), rep, charge > 0 ? kRed + 2 : kBlue + 2, 2, drawTrackMarkers, false, 3);
      }

      // draw detectors if option is set, only important for wire hits ----------------------
      if(drawDetectors) {

        if(wire_hit) {
          TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
          det_shape->IncDenyDestroy();
          det_shape->SetShape(new TGeoTube(std::max(0., (double)(hit_u-0.0105/2.)), hit_u+0.0105/2., plane_size));

          TVector3 norm = u.Cross(v);
          TGeoRotation det_rot("det_rot", (u.Theta()*180)/TMath::Pi(), (u.Phi()*180)/TMath::Pi(),
              (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi(),
              (v.Theta()*180)/TMath::Pi(), (v.Phi()*180)/TMath::Pi()); // move the tube to the right place and rotate it correctly
          TVector3 move = v*(v*(track_pos-o)); // move the tube along the wire until the track goes through it
          TGeoCombiTrans det_trans(o(0) + move.X(),
                                   o(1) + move.Y(),
                                   o(2) + move.Z(),
                                   &det_rot);
          det_shape->SetTransMatrix(det_trans);
          det_shape->SetMainColor(kCyan);
          det_shape->SetMainTransparency(25);
          if((drawHits && (hit_u+0.0105/2 > 0)) || !drawHits) {
            gEve->AddElement(det_shape);
          }
        }

      }
      // finished drawing detectors ---------------------------------------------------------

      if(drawHits) {

        // draw planar hits, with distinction between strip and pixel hits ----------------
        if(planar_hit) {
          if(!planar_pixel_hit) {
            TEveBox* hit_box;
            hit_box = boxCreator((o + hit_u*u), u, v, errorScale_*std::sqrt(hit_cov(0,0)), plane_size, 0.0105);
            hit_box->SetMainColor(kYellow);
            hit_box->SetMainTransparency(0);
            gEve->AddElement(hit_box);
          } else {
            // calculate eigenvalues to draw error-ellipse ----------------------------
            TMatrixDEigen eigen_values(hit_cov);
            TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
            cov_shape->IncDenyDestroy();
            TMatrixT<double> ev = eigen_values.GetEigenValues();
            TMatrixT<double> eVec = eigen_values.GetEigenVectors();
            double pseudo_res_0 = errorScale_*std::sqrt(ev(0,0));
            double pseudo_res_1 = errorScale_*std::sqrt(ev(1,1));
            // finished calcluating, got the values -----------------------------------

            // do autoscaling if necessary --------------------------------------------
            if(drawAutoScale) {
              double min_cov = std::min(pseudo_res_0,pseudo_res_1);
              if(min_cov < 1e-5) {
                std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
              } else {
                if(min_cov < 0.049) {
                  double cor = 0.05 / min_cov;
                  std::cout << "Track " << i << ", Hit " << j << ": Pixel covariance too small, rescaling by " << cor;
                  errorScale_ *= cor;
                  pseudo_res_0 *= cor;
                  pseudo_res_1 *= cor;
                  std::cout << " to " << errorScale_ << std::endl;
                }
              }
            }
            // finished autoscaling ---------------------------------------------------

            // calculate the semiaxis of the error ellipse ----------------------------
            cov_shape->SetShape(new TGeoEltu(pseudo_res_0, pseudo_res_1, 0.0105));
            TVector3 pix_pos = o + hit_u*u + hit_v*v;
            TVector3 u_semiaxis = (pix_pos + eVec(0,0)*u + eVec(1,0)*v)-pix_pos;
            TVector3 v_semiaxis = (pix_pos + eVec(0,1)*u + eVec(1,1)*v)-pix_pos;
            TVector3 norm = u.Cross(v);
            // finished calculating ---------------------------------------------------

            // rotate and translate everything correctly ------------------------------
            TGeoRotation det_rot("det_rot", (u_semiaxis.Theta()*180)/TMath::Pi(), (u_semiaxis.Phi()*180)/TMath::Pi(),
                (v_semiaxis.Theta()*180)/TMath::Pi(), (v_semiaxis.Phi()*180)/TMath::Pi(),
                (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi());
            TGeoCombiTrans det_trans(pix_pos(0),pix_pos(1),pix_pos(2), &det_rot);
            cov_shape->SetTransMatrix(det_trans);
            // finished rotating and translating --------------------------------------

            cov_shape->SetMainColor(kYellow);
            cov_shape->SetMainTransparency(0);
            gEve->AddElement(cov_shape);
          }
        }
        // finished drawing planar hits ---------------------------------------------------

        // draw spacepoint hits -----------------------------------------------------------
        if(space_hit) {

          // get eigenvalues of covariance to know how to draw the ellipsoid ------------
          TMatrixDEigen eigen_values(m->getRawHitCov());
          TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
          cov_shape->IncDenyDestroy();
          cov_shape->SetShape(new TGeoSphere(0.,1.));
          TMatrixT<double> ev = eigen_values.GetEigenValues();
          TMatrixT<double> eVec = eigen_values.GetEigenVectors();
          TVector3 eVec1(eVec(0,0),eVec(1,0),eVec(2,0));
          TVector3 eVec2(eVec(0,1),eVec(1,1),eVec(2,1));
          TVector3 eVec3(eVec(0,2),eVec(1,2),eVec(2,2));
          TVector3 norm = u.Cross(v);
          // got everything we need -----------------------------------------------------


          TGeoRotation det_rot("det_rot", (eVec1.Theta()*180)/TMath::Pi(), (eVec1.Phi()*180)/TMath::Pi(),
              (eVec2.Theta()*180)/TMath::Pi(), (eVec2.Phi()*180)/TMath::Pi(),
              (eVec3.Theta()*180)/TMath::Pi(), (eVec3.Phi()*180)/TMath::Pi()); // the rotation is already clear

          // set the scaled eigenvalues -------------------------------------------------
          double pseudo_res_0 = errorScale_*std::sqrt(ev(0,0));
          double pseudo_res_1 = errorScale_*std::sqrt(ev(1,1));
          double pseudo_res_2 = errorScale_*std::sqrt(ev(2,2));
          if(drawScaleMan) { // override again if necessary
            pseudo_res_0 = errorScale_*0.5;
            pseudo_res_1 = errorScale_*0.5;
            pseudo_res_2 = errorScale_*0.5;
          }
          // finished scaling -----------------------------------------------------------

          // autoscale if necessary -----------------------------------------------------
          if(drawAutoScale) {
            double min_cov = std::min(pseudo_res_0,std::min(pseudo_res_1,pseudo_res_2));
            if(min_cov < 1e-5) {
              std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
            } else {
              if(min_cov <= 0.149) {
                double cor = 0.15 / min_cov;
                std::cout << "Track " << i << ", Hit " << j << ": Space hit covariance too small, rescaling by " << cor;
                errorScale_ *= cor;
                pseudo_res_0 *= cor;
                pseudo_res_1 *= cor;
                pseudo_res_2 *= cor;
                std::cout << " to " << errorScale_ << std::endl;

              }
            }
          }
          // finished autoscaling -------------------------------------------------------

          // rotate and translate -------------------------------------------------------
          TGeoGenTrans det_trans(o(0),o(1),o(2),
                                 //std::sqrt(pseudo_res_0/pseudo_res_1/pseudo_res_2), std::sqrt(pseudo_res_1/pseudo_res_0/pseudo_res_2), std::sqrt(pseudo_res_2/pseudo_res_0/pseudo_res_1), // this workaround is necessary due to the "normalization" performed in  TGeoGenTrans::SetScale
                                 //1/(pseudo_res_0),1/(pseudo_res_1),1/(pseudo_res_2),
                                 pseudo_res_0, pseudo_res_1, pseudo_res_2,
                                 &det_rot);
          cov_shape->SetTransMatrix(det_trans);
          // finished rotating and translating ------------------------------------------

          cov_shape->SetMainColor(kYellow);
          cov_shape->SetMainTransparency(10);
          gEve->AddElement(cov_shape);
        }
        // finished drawing spacepoint hits -----------------------------------------------

        // draw wire hits -----------------------------------------------------------------
        if(wire_hit) {
          TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
          cov_shape->IncDenyDestroy();
          double pseudo_res_0 = errorScale_*std::sqrt(hit_cov(0,0));
          double pseudo_res_1 = plane_size;
          if (wirepoint_hit) pseudo_res_1 = errorScale_*std::sqrt(hit_cov(1,1));

          // autoscale if necessary -----------------------------------------------------
          if(drawAutoScale) {
            if(pseudo_res_0 < 1e-5) {
              std::cout << "Track " << i << ", Hit " << j << ": Invalid wire resolution (< 1e-5), autoscaling not possible!" << std::endl;
            } else {
              if(pseudo_res_0 < 0.0049) {
                double cor = 0.005 / pseudo_res_0;
                std::cout << "Track " << i << ", Hit " << j << ": Wire covariance too small, rescaling by " << cor;
                errorScale_ *= cor;
                pseudo_res_0 *= cor;
                std::cout << " to " << errorScale_ << std::endl;
              }
            }

            if(wirepoint_hit && pseudo_res_1 < 1e-5) {
              std::cout << "Track " << i << ", Hit " << j << ": Invalid wire resolution (< 1e-5), autoscaling not possible!" << std::endl;
            } else {
              if(pseudo_res_1 < 0.0049) {
                double cor = 0.005 / pseudo_res_1;
                std::cout << "Track " << i << ", Hit " << j << ": Wire covariance too small, rescaling by " << cor;
                errorScale_ *= cor;
                pseudo_res_1 *= cor;
                std::cout << " to " << errorScale_ << std::endl;
              }
            }
          }
          // finished autoscaling -------------------------------------------------------

          cov_shape->SetShape(new TGeoTube(std::max(0., (double)(hit_u - pseudo_res_0)), hit_u + pseudo_res_0, pseudo_res_1));
          TVector3 norm = u.Cross(v);

          // rotate and translate -------------------------------------------------------
          TGeoRotation det_rot("det_rot", (u.Theta()*180)/TMath::Pi(), (u.Phi()*180)/TMath::Pi(),
              (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi(),
              (v.Theta()*180)/TMath::Pi(), (v.Phi()*180)/TMath::Pi());
          TGeoCombiTrans det_trans(o(0) + hit_v*v.X(),
                                   o(1) + hit_v*v.Y(),
                                   o(2) + hit_v*v.Z(),
                                   &det_rot);
          cov_shape->SetTransMatrix(det_trans);
          // finished rotating and translating ------------------------------------------

          cov_shape->SetMainColor(kYellow);
          cov_shape->SetMainTransparency(50);
          gEve->AddElement(cov_shape);
        }
        // finished drawing wire hits -----------------------------------------------------
      }

      prevFi = fi;
      prevFittedState = fittedState;

    }

  }

  gEve->Redraw3D(kTRUE);

}




TEveBox* EventDisplay::boxCreator(TVector3 o, TVector3 u, TVector3 v, float ud, float vd, float depth) {

  TEveBox* box = new TEveBox("detPlane_shape");
  float vertices[24];

  TVector3 norm = u.Cross(v);
  u *= (0.5*ud);
  v *= (0.5*vd);
  norm *= (0.5*depth);

  vertices[0] = o(0) - u(0) - v(0) - norm(0);
  vertices[1] = o(1) - u(1) - v(1) - norm(1);
  vertices[2] = o(2) - u(2) - v(2) - norm(2);

  vertices[3] = o(0) + u(0) - v(0) - norm(0);
  vertices[4] = o(1) + u(1) - v(1) - norm(1);
  vertices[5] = o(2) + u(2) - v(2) - norm(2);

  vertices[6] = o(0) + u(0) - v(0) + norm(0);
  vertices[7] = o(1) + u(1) - v(1) + norm(1);
  vertices[8] = o(2) + u(2) - v(2) + norm(2);

  vertices[9] = o(0) - u(0) - v(0) + norm(0);
  vertices[10] = o(1) - u(1) - v(1) + norm(1);
  vertices[11] = o(2) - u(2) - v(2) + norm(2);

  vertices[12] = o(0) - u(0) + v(0) - norm(0);
  vertices[13] = o(1) - u(1) + v(1) - norm(1);
  vertices[14] = o(2) - u(2) + v(2) - norm(2);

  vertices[15] = o(0) + u(0) + v(0) - norm(0);
  vertices[16] = o(1) + u(1) + v(1) - norm(1);
  vertices[17] = o(2) + u(2) + v(2) - norm(2);

  vertices[18] = o(0) + u(0) + v(0) + norm(0);
  vertices[19] = o(1) + u(1) + v(1) + norm(1);
  vertices[20] = o(2) + u(2) + v(2) + norm(2);

  vertices[21] = o(0) - u(0) + v(0) + norm(0);
  vertices[22] = o(1) - u(1) + v(1) + norm(1);
  vertices[23] = o(2) - u(2) + v(2) + norm(2);


  for(int k = 0; k < 24; k += 3) box->SetVertex((k/3), vertices[k], vertices[k+1], vertices[k+2]);

  return box;

}

void EventDisplay::makeGui() {

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);

  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
  {
    TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    TGPictureButton* b = 0;
    EventDisplay*  fh = EventDisplay::getInstance();

    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
    hf->AddFrame(b);
    b->Connect("Pressed()", "genfit::EventDisplay", fh, "prev(UInt_t)");

    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
    hf->AddFrame(b);
    b->Connect("Clicked()", "genfit::EventDisplay", fh, "next(UInt_t)");
  }
  frmMain->AddFrame(hf);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);
}

}
