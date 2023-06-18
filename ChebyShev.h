/*  
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   ChebyShev.h
 * Copyright (C) 2018
 * Jenny Farmer jfarmer6@uncc.edu
 * Donald Jacobs djacobs1@uncc.edu
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in 
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CHEBYSHEV_HPP
#define	CHEBYSHEV_HPP

#include <vector>

using namespace std;
class ChebyShev {
public:
    ChebyShev();
    ChebyShev(const ChebyShev& orig);
    virtual ~ChebyShev();
    void initialize(double dzLocal[], int sizeLocal);
    void initializeDx(double dzLocal[], int sizeLocal);   

    vector < vector < double > >  getAllTerms(unsigned mode);
    vector < vector < double > >  getAllTermsDx(unsigned mode);
    
private:         
    int size;    
    double * dz;
    vector < vector<double> > termsT;
    vector < vector<double> > termsQ;
    
    vector <double> addMode(int mode);
    vector <double> addModeDx(int mode);
    
};

#endif	/* CHEBYSHEV_HPP */

