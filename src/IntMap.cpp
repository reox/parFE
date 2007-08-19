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

#include "IntMap.h"

std::ostream& IntMap::print(std::ostream& os) const
{
  os << "IntMap" << std::endl;
  if (begin() == end()) 
  {
    os << "[empty list]" << std::endl;
  }
  else 
    for (ConstIterator i = begin(); i != end(); ++i) 
    {
      os << (*i).second << " <-> " << (*i).first << std::endl;
    }
  return os;
}
