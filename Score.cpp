/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   Score.cpp
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


#include "Score.h"

Score::~Score() {
}


 
double Score::getTargetScore(double SURD) {
 
    vector<double>::iterator it;
    it = lower_bound (SURDs.begin(), SURDs.end(), SURD/100);
    unsigned index = it - SURDs.begin();
    
    if (index == SURDs.size()) {
        return scores[index - 1];
    } 
    if (index == 0) {
        return scores[0];
    } 
    
    double E1 = scores[index - 1];
    double E2 = scores[index];
    double P1 = SURDs[index - 1];
    double P2 = SURDs[index];
    double E = E1 + (SURD/100 - P1)*(E2 - E1)/(P2 - P1);
    return E;  
    
}

double Score::getConfidence(double score) {
 
    vector<double>::iterator it;
    it = lower_bound (scores.begin(), scores.end(), score);
    unsigned index = it - scores.begin();
    
    if (index == scores.size()) {
        return SURDs[index - 1];
    } 
    if (index == 0) {
        return SURDs[0];
    } 
    
    double E1 = scores[index - 1];
    double E2 = scores[index];
    double P1 = SURDs[index - 1];
    double P2 = SURDs[index];
    double P = P1 + (score - E1)*(P2 - P1)/(E2 - E1);
    return P * 100;  
    
}
