/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006, Uche Mennel and Marzio Sala
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301, USA.
 */

#ifndef _INTMAP_H_
#define _INTMAP_H_

#include <map>
#include <iostream>

//! A class to generate printable maps

/*! For every type inserted into a Teuchos::Parameterlist, there must exist a method to print this type
    in order to print the whole ParameterList.
*/
class IntMap : public std::map<int, int> {
  
  typedef IntMap::iterator Iterator;
  
public:
  typedef IntMap::const_iterator ConstIterator;

  //! Printing method for this map
  std::ostream& print(std::ostream& os) const;
  
};

//! Printing operator for type IntMap
inline std::ostream& operator<<(std::ostream& os, const IntMap& l)
{
  return l.print(os);
}

#endif
