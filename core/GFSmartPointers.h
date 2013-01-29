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

#ifndef GFSMARTPOINTERS_H
#define GFSMARTPOINTERS_H

#if defined (__CINT__)
#  define UNIQUE_PTR(T) T*
#  define SHARED_PTR(T) T*
#else
#  ifdef __GLIBCXX__
#    include <tr1/memory>
#  else
#    ifdef __IBMCPP__
#      define __IBMCPP_TR1__
#    endif
#    include <memory>
#  endif
#  define UNIQUE_PTR(T) std::tr1::unique_ptr<T>
#  define SHARED_PTR(T) std::tr1::shared_ptr<T>
#endif

/** @brief ROOT CINT compatible C++ 11 smart pointers
 *
 */

#endif
