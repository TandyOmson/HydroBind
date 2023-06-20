//    Copright (C) 1999-2021, Bernd Gaertner
//    November 12, 2021
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//    Contact:
//    --------
//    Bernd Gaertner
//    Institute of Theoretical Computer Science 
//    ETH Zuerich
//    CAB G31.1
//    CH-8092 Zuerich, Switzerland
//    http://www.inf.ethz.ch/personal/gaertner

// =========================================================================
// Container for Miniball
// Finds the smallest sphere that encloses a set of points in 3D
//  - Input: a list of tuples, each containing three doubles
// Called from Python using pybind11
// =========================================================================

#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include "Miniball.hpp"

// Include pybind11 headers
#include<pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;

using std::vector;
using std::tuple;

typedef double mytype; // coordinate type for Miniball

// Create a type for a list of tuples, each containing three doubles
typedef vector< tuple<double,double,double> > tuple_list;
// Create a type for a list of vectors of doubles
typedef std::list<std::vector<mytype> > vector_list;


// Create a function that takes a list of tuples (each with three doubles) and returns a list of vectors (three doubles)
vector_list vectors_from_tuples (tuple_list &initlist) { 
  int d = 3;
  int n = initlist.size();
  std::list<std::vector<mytype> > lp;
  for (int i=0; i<n; ++i) {     
    std::vector<mytype> p(d);
    p[0] = std::get<0>(initlist.at(i));
    p[1] = std::get<1>(initlist.at(i));
    p[2] = std::get<2>(initlist.at(i));
    lp.push_back(p);
  }
  return lp;
};

// Effectively the main function for Miniball, to be called from Python
// Takes a set of coordinates, and calculates the smallest sphere that encloses all of them
std::vector<mytype> get_smallest_enclosing_sphere (tuple_list &initlist)
{
  int d = 3;
  typedef double mytype;             // coordinate type
  // Declare a list of vectors of doubles
  std::list<std::vector<mytype> > lp;
  // Convert the input list of tuples to a list of vectors and assign it to lp
  lp = vectors_from_tuples(initlist);

  // define the types of iterators through the points and their coordinates
  // ----------------------------------------------------------------------
  typedef std::list<std::vector<mytype> >::const_iterator PointIterator; 
  typedef std::vector<mytype>::const_iterator CoordIterator;

  // create an instance of Miniball, run it
  // ------------------------------
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
    MB;
  // MB object, d dimensions, lp.begin() and lp.end() as iterators  
  MB mb (d, lp.begin(), lp.end());
  
  // output results
  // --------------
  // Declare a vector of doubles to store the center and squared radius of the sphere
  std::vector<mytype>  results;
  const mytype* center = mb.center(); 
  for(int i=0; i<d; ++i, ++center){
    results.push_back(*center);
  }
  results.push_back(mb.squared_radius());

  return results;
}

// Define the module name
PYBIND11_MODULE(miniball_container, m) {
  // Docstring
  m.doc() = "pybind11 plugin for Miniball";
  // Define the function name and attach C++ function to be called from Python
  m.def("get_smallest_enclosing_sphere", &get_smallest_enclosing_sphere, "A function that takes a list of tuples (each with three doubles) and returns a list of vectors (three doubles)");
}

// To compile: g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)
// Creates a .so file (shared object file that contains compiled code and data that can be dynamically linked to an executable at runtime)
// Importable in python:
// from miniball_container import get_smallest_enclosing_sphere
