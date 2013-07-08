/* Copyright 2013, Ludwig-Maximilians Universität München,
   Authors: Tobias Schlüter & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_DAF_h
#define genfit_DAF_h

#include "AbsKalmanFitter.h"
#include <vector>
#include <map>

namespace genfit {

/** @brief Determinstic Annealing Filter (DAF) implementation.
 *
 * @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 * @author Karl Bicker (Technische Universit&auml;t M&uuml;nchen)
 *
 * The DAF is an iterative Kalman filter with annealing. It is capable of
 * fitting tracks which are contaminated with noise hits. The algorithm is
 * taken from the references R. Fruehwirth & A. Strandlie, Computer Physics
 * Communications 120 (1999) 197-214 and CERN thesis: Dissertation by Matthias
 * Winkler.
 *
 * The weights which were assigned to the hits by the DAF are accessible by using the
 * bookkeeping object of the fitted track. The weight is stored under the key "dafWeight".
 * So to retrieve for example the weight of hit 10, fitted with track representation 2,
 * use GFTrack::getBK(2)->getNumber("dafWeight", 10,  double& wght).
 */

class DAF : public AbsKalmanFitter {

  DAF();
  DAF(AbsKalmanFitter* kalman);
  ~DAF() {};

  /** @brief Process a track using the DAF.
   */
  void processTrack(Track* tr, const AbsTrackRep* rep);

  /** @brief Set the probability cut for the weight calculation for the hits.
   *
   * By default the cut values for measurements of dimensionality from 1 to 5 are calculated.
   * If you what to have cut values for an arbitrary measurement dimensionality use
   * addProbCut(double prob_cut, int maxDim);
   */
  void setProbCut(const double prob_cut);

  /** Set the probability cut for the weight calculation for the hits for a specific measurement dimensionality*/
  void addProbCut(const double prob_cut, const int measDim);

  /** @brief Configure the annealing scheme.
   *
   * In the current implementation you need to provide at least one temperatures
   * and not more then ten temperatures.
   */
  void setBetas(double b1,double b2=-1, double b3=-1., double b4=-1., double b5=-1., double b6=-1., double b7=-1., double b8=-1., double b9=-1., double b10=-1.);

  void resolveWireHitAmbi(bool resolve = true);


 private:

  /** @brief check if convergence is met so the main iteration (over the betas) can be left
   * the convergence criteria is the largest change in the weights
   */
  bool isConvergent(const std::vector<std::vector<double> >& oldWeights, AbsTrackRep* rep) const;

  /** @brief Calculate the weights for the next fitting pass.
    */
  std::vector<std::vector<double> > calcWeights(Track* trk, double beta);

  std::map<AbsTrackRep*, std::vector<std::vector<double> > > weights_;
  std::vector<double> betas_;
  std::map<int,double>  chi2Cuts_;
  boost::scoped_ptr<AbsKalmanFitter> kalman_;

  /** The maximal number of iterations in the main DAF loop. If the weights do not converge the loop will end after c_maxIter iterations*/
  const static int c_maxIter = 10;
};

}  /* End of namespace genfit */
/** @} */

#endif //genfit_DAF_h
