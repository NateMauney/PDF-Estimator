/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   Score.h
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

#ifndef SCORE_HPP
#define	SCORE_HPP

#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>

using namespace std;

class Score {
public:        
    Score() {};
    Score(const Score& orig) {};
    virtual ~Score();
    
    virtual double calculateScorePartition(double r[], int p) {return 0;}
    virtual double calculateScore(double r[], int N) {return 0;}
    virtual vector <int>  setIndices (int N, int p, double * data, bool index) {return indices;}
    virtual vector <int> getIndices (int N, int p, double * data) {return indices;}
    
    double getLikelihood() {return likelihood;};
//    void setVarianceMin(bool qzVar) {minimizeVariance = qzVar;}
    double getConfidence(double score);    
    double QZVariance;
        
    double getTargetScore(double SURD);  
    
protected:
    vector <double> scores;
    vector <double> SURDs;
    double likelihood;
    vector <int> indices;
            
    void getValues();
    
};

#endif	/* SCORE_HPP */

