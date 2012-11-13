/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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

#include "GFBookkeeping.h"
#include "GFException.h"
#include <TString.h>


GFBookkeeping::GFBookkeeping(const GFBookkeeping& bk) {
  fNhits = bk.fNhits;
  fVectors = bk.fVectors;
  fMatrices = bk.fMatrices;
  fSymMatrices = bk.fSymMatrices;
  fNumbers = bk.fNumbers; 
  fPlanes = bk.fPlanes; 
  fFailedHits = bk.fFailedHits;
}

// Use a custom streamer, the auto-generated one is prohibitively slower (root 5.34).
void GFBookkeeping::Streamer(TBuffer &R__b)
{
   // Stream an object of class GFBookkeeping.
   if (R__b.IsReading()) {

     //     Version_t R__v = R__b.ReadVersion(); 
     TObject::Streamer(R__b); 

     clearAll();

     R__b >> fNhits;

     TString s;
     std::string key;
     unsigned int nkeys;
     TVectorT<double> vec;
     TMatrixT<double> mat;
     TMatrixTSym<double> symmat;
     GFDetPlane pl;

     {//reading vectors
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         s.Streamer(R__b);
         key = s.Data();
         bookVectors(key);
         for(int j=0;j<fNhits;++j){
           vec.Streamer(R__b);
           setVector(key,j,vec);
         }
       }
     }//done reading vectors
     {//reading matrices
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         s.Streamer(R__b);
         key = s.Data();
         bookMatrices(key);
         for(int j=0;j<fNhits;++j){
           mat.Streamer(R__b);
           setMatrix(key,j,mat);
         }
       }
     }//done reading matrices
     {//reading symmetric matrices
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         s.Streamer(R__b);
         key = s.Data();
         bookSymMatrices(key);
         for(int j=0;j<fNhits;++j){
           symmat.Streamer(R__b);
           setSymMatrix(key,j,symmat);
         }
       }
     }//done reading matrices
     {//reading planes
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         s.Streamer(R__b);
         key = s.Data();
         bookGFDetPlanes(key);
         for(int j=0;j<fNhits;++j){
           pl.Streamer(R__b);
           setDetPlane(key,j,pl);
         }
       }
     }//done reading planes
     {//reading numbers
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         s.Streamer(R__b);
         key = s.Data();
         bookNumbers(key);
         for(int j=0;j<fNhits;++j){
           double d;
           R__b >> d;
           setNumber(key,j,d);
         }
       }
     }//done reading numbers
     {//read failed hits
       clearFailedHits();
       unsigned int nFailedHits;
       R__b >> nFailedHits;
       fFailedHits.reserve(nFailedHits);
       for(unsigned int i=0;i<nFailedHits;++i){
         unsigned int aFailedHit;
         R__b >> aFailedHit;
         fFailedHits[i] = aFailedHit;
       }
     }//done reading failed hits
   } else {
     //     R__b.WriteVersion(GFBookkeeping::IsA()); 
     TObject::Streamer(R__b); 

     //write number of hits
     R__b << fNhits;

     std::vector<std::string> keys;
     {//save vectors
       keys = getVectorKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         TString s(keys.at(i));
         s.Streamer(R__b);
         for(int j=0;j<fNhits;++j){
           ((fVectors[keys.at(i)])[j]).Streamer(R__b);
         }
       }
     }
     {//save matrices
       keys = getMatrixKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         TString s(keys.at(i));
         s.Streamer(R__b);
         for(int j=0;j<fNhits;++j){
           ((fMatrices[keys.at(i)])[j]).Streamer(R__b);
         }
       }
     }
     {//save symmetric matrices
       keys = getSymMatrixKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         TString s(keys.at(i));
         s.Streamer(R__b);
         for(int j=0;j<fNhits;++j){
           ((fSymMatrices[keys.at(i)])[j]).Streamer(R__b);
         }
       }
     }
     keys.clear();
     {//save GFDetPlanes
       keys = getGFDetPlaneKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         TString s(keys.at(i));
         s.Streamer(R__b);
         for(int j=0;j<fNhits;++j){
           ((fPlanes[keys.at(i)])[j]).Streamer(R__b);
         }
       }
     }//done saving GFDetPlanes
     keys.clear();
     {//save numbers    
       keys = getNumberKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         TString s(keys.at(i));
         s.Streamer(R__b);
         for(int j=0;j<fNhits;++j){
	   R__b << (fNumbers[keys.at(i)])[j];
         }
       }
     }//done saving numbers 
     {//save failedHits
       R__b << ((unsigned int) fFailedHits.size());
       for(unsigned int i=0;i<fFailedHits.size();++i){
         R__b << fFailedHits.at(i);
       }
     }//done saving failed Hits    
   }
}


void GFBookkeeping::bookVectors(const std::string& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<std::string, std::vector<TVectorT<double> > >::const_iterator it;
  it = fVectors.find(key);
  if(it != fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookVectors()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fVectors[key].resize(fNhits);
}

void GFBookkeeping::bookMatrices(const std::string& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<std::string, std::vector<TMatrixT<double> > >::const_iterator it;
  it = fMatrices.find(key);
  if(it != fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookMatrices()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fMatrices[key].resize(fNhits);
}

void GFBookkeeping::bookSymMatrices(const std::string& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<std::string, std::vector<TMatrixTSym<double> > >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it != fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookSymMatrices()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fSymMatrices[key].resize(fNhits);
}

void GFBookkeeping::bookGFDetPlanes(const std::string& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<std::string, std::vector<GFDetPlane > >::const_iterator it;
  it = fPlanes.find(key);
  if(it != fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookGFDetPlanes()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fPlanes[key].resize(fNhits);
}

//val is default set to 0.
void GFBookkeeping::bookNumbers(const std::string& key,double val){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<std::string, std::vector<double> >::const_iterator it;
  it = fNumbers.find(key);
  if(it != fNumbers.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookNumbers()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fNumbers[key].resize(fNhits);
  for(int i=0;i<fNhits;++i){
    (fNumbers[key])[i] = val;
  }
  
}

void GFBookkeeping::setVector(const std::string& key, unsigned int index,
			      const TVectorT<double>& vec){
  std::map<std::string, std::vector<TVectorT<double> > >::const_iterator it;
  it = fVectors.find(key);
  if(it == fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setMatrix()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fVectors[key])[index].ResizeTo(vec);
  (fVectors[key])[index] = vec;
}

void GFBookkeeping::setMatrix(const std::string& key, unsigned int index,
			    const TMatrixT<double>& mat){
  std::map<std::string, std::vector<TMatrixT<double> > >::const_iterator it;
  it = fMatrices.find(key);
  if(it == fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setMatrix()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fMatrices[key])[index].ResizeTo(mat);
  (fMatrices[key])[index] = mat;
}

void GFBookkeeping::setSymMatrix(const std::string& key, unsigned int index,
				 const TMatrixTSym<double>& mat){
  std::map<std::string, std::vector<TMatrixTSym<double> > >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it == fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setMatrix()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fSymMatrices[key])[index].ResizeTo(mat);
  (fSymMatrices[key])[index] = mat;
}

void GFBookkeeping::setDetPlane(const std::string& key, unsigned int index,
				  const GFDetPlane& pl){
  std::map<std::string, std::vector<GFDetPlane> >::const_iterator it;
  it = fPlanes.find(key);
  if(it == fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fPlanes[key])[index] = pl;
}
void GFBookkeeping::setNumber(const std::string& key, unsigned int index,
			      const double& num){
  if(fNumbers[key].size() == 0){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fNumbers[key])[index] = num;
}

bool GFBookkeeping::getVector(const std::string& key,
			    unsigned int index,
			    TVectorT<double>& vec) const {
  std::map<std::string, std::vector<TVectorT<double> > >::const_iterator it;
  it = fVectors.find(key);
  if(it == fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  vec.ResizeTo(((*it).second)[index]);
  vec = ((*it).second)[index];
  return true;
}

bool GFBookkeeping::getMatrix(const std::string& key,
			    unsigned int index,
			    TMatrixT<double>& mat) const {
  std::map<std::string, std::vector<TMatrixT<double> > >::const_iterator it;
  it = fMatrices.find(key);
  if(it == fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  mat.ResizeTo(((*it).second)[index]);
  mat = ((*it).second)[index];
  return true;
}

bool GFBookkeeping::getSymMatrix(const std::string& key,
			    unsigned int index,
			    TMatrixTSym<double>& mat) const {
  std::map<std::string, std::vector<TMatrixTSym<double> > >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it == fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  mat.ResizeTo(((*it).second)[index]);
  mat = ((*it).second)[index];
  return true;
}

bool GFBookkeeping::getDetPlane(const std::string& key,
			      unsigned int index,
			      GFDetPlane& pl) const {
  std::map<std::string, std::vector<GFDetPlane> >::const_iterator it;
  it = fPlanes.find(key);
  if(it == fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  pl = ((*it).second)[index];
  return true;
}

bool GFBookkeeping::getNumber(const std::string& key,
			    unsigned int index,
			    double& num) const {
  std::map<std::string, std::vector<double> >::const_iterator it;
  it = fNumbers.find(key);
  if(it == fNumbers.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  num = ((*it).second)[index];
  return true;
}

void GFBookkeeping::addFailedHit(unsigned int id){
  fFailedHits.push_back( id );
}

unsigned int GFBookkeeping::getNumFailed(){
  return fFailedHits.size();
}

unsigned int GFBookkeeping::hitFailed(unsigned int id){
  unsigned int retVal = 0;
  for(unsigned int i=0; i<fFailedHits.size(); ++i){
    if(fFailedHits.at(i) == id){
      ++retVal;
    }
  }
  return retVal;
}

void GFBookkeeping::clearFailedHits(){
  fFailedHits.clear();
}

void GFBookkeeping::reset() {
  std::vector<std::string> vecKeys = getVectorKeys();
  std::vector<std::string> matKeys = getMatrixKeys();
  std::vector<std::string> symMatKeys = getSymMatrixKeys();
  std::vector<std::string> planeKeys = getGFDetPlaneKeys();
  std::vector<std::string> numKeys = getNumberKeys();

  clearAll();
  clearFailedHits();

  for(unsigned int i=0;i<vecKeys.size();++i){
    bookVectors(vecKeys.at(i));
  }
  for(unsigned int i=0;i<matKeys.size();++i){
    bookMatrices(matKeys.at(i));
  }
  for(unsigned int i=0;i<symMatKeys.size();++i){
    bookSymMatrices(symMatKeys.at(i));
  }
  for(unsigned int i=0;i<planeKeys.size();++i){
    bookGFDetPlanes(planeKeys.at(i));
  }
  for(unsigned int i=0;i<numKeys.size();++i){
    bookNumbers(numKeys.at(i));
  }

}

void GFBookkeeping::clearAll(){
  fVectors.clear();
  fMatrices.clear();
  fSymMatrices.clear();
  fPlanes.clear();
  fNumbers.clear();
}

std::vector< std::string > GFBookkeeping::getVectorKeys() const {
  std::vector< std::string > keys;
  keys.reserve(fVectors.size());
  std::map<std::string, std::vector<TVectorT<double> > >::const_iterator it;
  for(it=fVectors.begin();it!=fVectors.end();it++){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}
std::vector< std::string > GFBookkeeping::getMatrixKeys() const {
  std::vector< std::string > keys;
  keys.reserve(fMatrices.size());
  std::map<std::string, std::vector<TMatrixT<double> > >::const_iterator it;
  for(it=fMatrices.begin();it!=fMatrices.end();it++){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}
std::vector< std::string > GFBookkeeping::getSymMatrixKeys() const {
  std::vector< std::string > keys;
  keys.reserve(fSymMatrices.size());
  std::map<std::string, std::vector<TMatrixTSym<double> > >::const_iterator it;
  for(it=fSymMatrices.begin();it!=fSymMatrices.end();it++){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}
std::vector< std::string > GFBookkeeping::getGFDetPlaneKeys() const {
  std::vector< std::string > keys;
  keys.reserve(fPlanes.size());
  std::map<std::string, std::vector<GFDetPlane> >::const_iterator it;
  for(it=fPlanes.begin();it!=fPlanes.end();it++){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}
std::vector< std::string > GFBookkeeping::getNumberKeys() const {
  std::vector< std::string > keys;
  keys.reserve(fNumbers.size());
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it=fNumbers.begin();it!=fNumbers.end();it++){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}


void GFBookkeeping::Print(const Option_t* option) const {
  std::cout << "=============GFBookkeeping::print()==============\n";
  std::cout << "GFBookkeeping has " << fNhits << " hits.\n";
  std::cout << "-----printing all vectors:------\n";
  std::vector<std::string> keys = getVectorKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys.at(i) << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      TVectorT<double> m;
      getVector(keys.at(i),j,m);
      m.Print(option);
    }
  }
  std::cout << "-----printing all matrices:------\n";
  keys = getMatrixKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys.at(i) << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      TMatrixT<double> m;
      getMatrix(keys.at(i),j,m);
      m.Print(option);
    }
  }
  std::cout << "-----printing all symmetric matrices:------\n";
  keys = getSymMatrixKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys.at(i) << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      TMatrixTSym<double> m;
      getSymMatrix(keys.at(i),j,m);
      m.Print(option);
    }
  }
  std::cout << "-----printing all GFDetPlanes:------\n";
  keys = getGFDetPlaneKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys.at(i) << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      GFDetPlane p;
      getDetPlane(keys.at(i),j,p);
      p.Print(option);
    }
  }
  std::cout << "-----printing all numbers:------\n";
  keys = getNumberKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys.at(i) << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      double n(-1111.);
      getNumber(keys.at(i),j,n);
      std::cout << n << std::endl;
    }
  }
  std::cout << "-----failed hits:------\n";
  for(unsigned int i=0;i<fFailedHits.size();++i){
    std::cout << fFailedHits.at(i) << " ";
  }
  std::cout << "==========GFBookkeeping::print() - Done==========" << std::endl;
}


ClassImp(GFBookkeeping)
