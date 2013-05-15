/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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


#ifndef genfit_FieldManager_h
#define genfit_FieldManager_h

#include "AbsBField.h"
#include <iostream>
#include <stdexcept>
#include <string>

namespace genfit {

/** @brief Singleton which provides access to magnetic field for track representations
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */

class FieldManager {

 public:

  AbsBField* getField(){
    if(field_==NULL){
      std::cerr << "FieldManager hasn't been initialized with a correct AbsBField pointer!" << std::endl;
      std::string msg("FieldManager hasn't been initialized with a correct AbsBField pointer!");
      std::runtime_error err(msg);
      throw err;
    }
    return field_;
  }

  static TVector3 getFieldVal(const TVector3& position){
    if(instance_==NULL){
      std::cerr << "FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!" << std::endl;
      std::string msg("FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!");
      std::runtime_error err(msg);
      throw err;
    }
    if(field_==NULL){
      std::cerr << "FieldManager hasn't been initialized with a correct AbsBField pointer!" << std::endl;
      std::string msg("FieldManager hasn't been initialized with a correct AbsBField pointer!");
      std::runtime_error err(msg);
      throw err;
    }
    return field_->get(position);
  }

  static void getFieldVal(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz){
    if(instance_==NULL){
      std::cerr << "FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!" << std::endl;
      std::string msg("FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!");
      std::runtime_error err(msg);
      throw err;
    }
    if(field_==NULL){
      std::cerr << "FieldManager hasn't been initialized with a correct AbsBField pointer!" << std::endl;
      std::string msg("FieldManager hasn't been initialized with a correct AbsBField pointer!");
      std::runtime_error err(msg);
      throw err;
    }
    return field_->get(posX, posY, posZ, Bx, By, Bz);
  }

  //! set the magnetic field here. Magnetic field classes must be derived from AbsBField.
  void init(AbsBField* b) {
    field_=b;
  }

  static FieldManager* getInstance(){
    if(instance_==NULL) {
      instance_ = new FieldManager();
    }
    return instance_;
  }


 private:

  FieldManager(){}
  static FieldManager* instance_;
  static AbsBField* field_;

};

} /* End of namespace genfit */

#endif // genfit_FieldManager_h
