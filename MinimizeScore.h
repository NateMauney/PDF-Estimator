/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   MinimizeScore.h
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

#ifndef MINIMIZESCORE_HPP
#define	MINIMIZESCORE_HPP

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>

#include "ChebyShev.h"
#include "InputData.h"
#include "ScoreQZ.h"
#include "OutputControl.h"
#include <limits>
#include <climits>

#ifdef outputR
#include "Rmath.h"
#endif

//#define clock

using namespace std;

class MinimizeScore {
public:
    MinimizeScore();
    virtual ~MinimizeScore();
    bool minimize(const InputParameters& input, const InputData& data);
    vector <double> getLagrange();    
    
    OutputControl out;
    
    double bestThreshold;
    double * bestRandom;
    
    int    mode;
    float duration;
    double  normalize;
    int     N;    
    int     initPartitionSize;
    int     partitionSize;
    int     targetPartition;
private:
    int     nPoints;
    int maxLagrange;
    int seed;
    bool useLast;
    double y2;    
    ChebyShev cheby;
    double * inverse;
    double * z;
    double * dzWeight1;
    double * dzWeight2;
    double * dzWeight3; 
    double * trialRandom;
    double bestScore;
    double * bestLagrange;
    double * rawDataPartition;
    vector < vector < double > > T;
    vector < vector < double > > Tdx;
    vector <int> smoothWindow;
    vector <double> smoothSize;
    double smoothError;
    bool smooth;
   
    
    void funnelDiffusion(double * original, double * updated, int arraySize, double currentSigmaMu);
    void funnelDiffusion(double * original, double * updated, int arraySize, double currentSigmaMu, int startIndex);
    double funnelDiffusion(double original, double currentSigmaMu);
    double random(double m, double s);
    double ranX();
    void map(double r[], double cdf[], double rawDataPartition[], int partitionSize);
    void calculatePDF (double cdf[], double lagrange[], int modes);
    void calculatePDFAdaptive (double cdf[], double lagrange[], int modes);
    void calculatePDF (double cdf[], double power);
};

#endif	/* MINIMIZESCORE_HPP */

