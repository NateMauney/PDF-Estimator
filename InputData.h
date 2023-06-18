/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   inputData.h
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

#ifndef INPUTDATA_HPP
#define	INPUTDATA_HPP

#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include "OutputControl.h"
#include "InputParameters.h"
#include "ChebyShev.h"

using namespace std;

class InputData {
public:    
    virtual ~InputData();
    double * inverse;
    double * doubleInverse;
    double * transformedZeroOne;
    double * dz;
    double * dzWeight1;
    double * dzWeight2;
    double * dzWeight3;
    int N;
    int nPoints;
    double minimumRaw;
    double maximumRaw;
    double minimumCalc;
    double maximumCalc;
    int nRightOutliers;
    int nLeftOutliers;
    int nPointsAdjust;
    ChebyShev cheby;
    
    vector <int> smoothWindow;
    vector <double> smoothSize;
        
    InputData() {};
    InputData(const InputParameters& input);
//    InputData(const InputData& orig);
    bool readData();
    void setData(vector <double> & data);
    vector<int> realIdx;
    bool processData();
    OutputControl out;   
private: 
    InputParameters input; 
    
    bool leftOutliers;
    bool rightOutliers;
           
    vector <double> rawData;
    vector <double> transformedData;
    vector <double> tempData;
    bool transformData();
    void setAdaptiveDz();
//    void setUniformDz();
    void identifyOutliers(); 
};

#endif	/* INPUTDATA_HPP */

