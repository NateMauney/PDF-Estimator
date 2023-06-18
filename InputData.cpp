/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   inputData.cpp
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

#include <numeric>

#include "InputData.h"
#include "WriteResults.h"


InputData::InputData(const InputParameters& input) {
    this->input = input;
        
    N = 0;
    nPoints = 0;
    minimumRaw = 0;
    maximumRaw = 0;
    minimumCalc = 0;
    maximumCalc = 0;
    nPointsAdjust = 0;
    
    nRightOutliers = 0;
    nLeftOutliers = 0;    
    leftOutliers = false;
    rightOutliers = false;
    
}
//InputData::InputData(const InputData& orig) {
//}

InputData::~InputData() {    
    delete [] doubleInverse;
    delete [] transformedZeroOne;
    delete [] inverse;
    delete [] dz;
    delete [] dzWeight1;
    delete [] dzWeight2;
    delete [] dzWeight3;
}

bool InputData::readData() {
    
    ifstream fin;
    string line;
    fin.open((input.inputPath + input.inputFile).c_str());

    if(!fin.is_open()){
        out.error("Failed to open data file " + input.inputFile);
        return false;
    }
	    
    while (getline(fin, line)) {
        double test = atof(line.c_str());
        if (test == 0) {
            test = 0;
        }
        rawData.push_back(test);
    }
    if (rawData.size() == 0) {
        out.error("No data in " + input.inputFile);
        return false;        
    }
   
    fin.close();       
    sort(rawData.begin(), rawData.end());
    return processData();
}
    
void InputData::setData(vector<double> & data) {       
    rawData.clear();
    rawData.reserve(data.size());
    rawData = data;
    sort(rawData.begin(), rawData.end());
    
}

bool InputData::processData() {   
    nPoints = input.integrationPoints;
    if (nPoints == -1) {
        nPoints = (int) (200 + rawData.size()/200.0);
        if (nPoints > 1500) nPoints = 1500;
    }
    
    minimumRaw = rawData[0];
    maximumRaw = rawData[rawData.size() - 1];
    if (minimumRaw == maximumRaw) {        
        out.error("All input data has the same value ", minimumRaw);
        return false;        
    }
    
    int nValues = rawData.size();
    if (input.upperBoundSpecified) {  
            maximumCalc = input.upperBound;
    } else {
        double max = rawData[nValues - 1];
        maximumCalc = max + (max - rawData[rawData.size() - 5]);
        if (maximumCalc < max) {
            maximumCalc = max;
        }
    }
    if (input.lowerBoundSpecified) {  
        minimumCalc = input.lowerBound;
    } else {
        double min = rawData[0];
        minimumCalc = min + (min - rawData[4]);
        if (minimumCalc > min) {
            minimumCalc = min;
        }
    }
    
    if (input.outlierCutoff > 0) {        
        identifyOutliers();
    }
    
    if (!transformData()) {
        return false;
    }
    setAdaptiveDz();  
    cheby.initialize(doubleInverse, 2*nPointsAdjust-1);
    cheby.initializeDx(doubleInverse, 2*nPointsAdjust-1);
    return true;
}


  void InputData::identifyOutliers() {
         
        double q1 = 0;
        double q3 = 0;      
        
        int nValues = rawData.size();    
        
        int middle = (int) (nValues/2);
        int quarter = (int) (middle/2);
        
        if (nValues%2 == 0) {
            if (middle%2 == 0) {
                q1 = (rawData[(int)quarter - 1] + rawData[(int)quarter])/2;
                q3 = (rawData[(int)quarter + (int)middle - 1] + rawData[(int)quarter + (int)middle])/2;
            } else {
                q1 = rawData[(int)quarter];
                q3 = rawData[(int)quarter + (int)middle];
            }
        } else {
            if (middle%2 == 0) {
                q1 = (rawData[(int)quarter - 1] + rawData[(int)quarter])/2;
                q3 = (rawData[(int)quarter + (int)middle] + rawData[(int)quarter + (int)middle + 1])/2;
            } else {
                q1 = rawData[(int)quarter];
                q3 = rawData[(int)quarter + (int) middle];
            }
        }               
        double iqr = input.outlierCutoff*(q3 - q1);
        double leftOutlier = q1 - iqr;
        double rightOutlier = q3 + iqr;    
                
        if (maximumCalc > rightOutlier) {
            maximumCalc = rightOutlier;
            rightOutliers = true;                
        }
        if (minimumCalc < leftOutlier) {
            minimumCalc = leftOutlier;
            leftOutliers = true;                
        }
  }
    


bool InputData::transformData() {    
    
    int nValues = rawData.size();       
            
    int nAllData = 0;
    for (vector<double>::iterator iter = rawData.begin(); iter != rawData.end(); ++iter) {
        nAllData++;
        if (*iter >= minimumCalc) {
            if (*iter <= maximumCalc) {
                tempData.push_back(*iter);
                realIdx.push_back(nAllData);
            } else {
                nRightOutliers++;
            }
        } else {
            nLeftOutliers++;
        }            
    }
        
    nValues = tempData.size(); 
    if (nValues == 0) {
        out.error("No data within specified boundaries");
        return false;
    }
    
    
    transformedData.clear();      
    transformedData.reserve(nValues);
    transformedZeroOne = new double[nValues];
    int count = 0;
    for (vector<double>::iterator iter = tempData.begin(); iter != tempData.end(); ++iter) {
        transformedData.push_back((2*(*iter) - maximumCalc - minimumCalc)/(maximumCalc - minimumCalc));
        transformedZeroOne[count] = (transformedData[count] + 1)/2.0;
        count++;    
    }
    return true;
}


void InputData::setAdaptiveDz() {
   
    vector <double> inverseVector;   
    N = transformedData.size();      
    double last = 0;
    double next;
    bool breakOut = false;
    double dzMax = 1.0/(nPoints - 1);
    double maxSmoothWindow = 10;
          
    int skip = (int) (N/(nPoints - 1));
    if (skip==0) skip = 1;
           
    inverseVector.push_back(0);    
    for (int b = 0; b <= (N + skip); b+=skip) {
        if (b >= (N-1)) {
            next = 1;
            breakOut = true;
        }
        else {
            next = transformedZeroOne[b];
        }
        double test = next - last;
        double difference = fabs(test);
        if (difference > dzMax) {
            double steps =  difference / dzMax;
            int iSteps = (int) steps;
            int wStepSize = iSteps + 1;
            double add = difference / (iSteps + 1);
            if (iSteps > maxSmoothWindow) {
                double windowSteps = (iSteps + 1) / maxSmoothWindow;
                int wSteps = ceil(windowSteps);
                wStepSize = (int) (iSteps * 1.0 + 1) / wSteps;
            }
            double sum = 0;
            int windows = 0;
            for (int k = 0; k < (iSteps + 1); k++) {
                inverseVector.push_back(inverseVector[inverseVector.size() - 1] + add);
                sum += add;
                windows ++;
                if ((windows) > wStepSize) {
                    smoothWindow.push_back((windows) * 2);
                    smoothSize.push_back(sum);
                    sum = 0;
                    windows = 0;
                }
            }
            if (windows > 0) {
                smoothWindow.push_back((windows) * 2);
                smoothSize.push_back(sum);
            }
        } else { 
            inverseVector.push_back(next);
            smoothWindow.push_back(2);
            smoothSize.push_back(difference);
        }
        if (breakOut) break;
        last = next; 
    }            
    
    inverseVector[inverseVector.size() - 1] = 1;
    int inverseSize = inverseVector.size();
    inverse = new double[inverseSize];
    doubleInverse = new double[2 * inverseSize - 1];
    dz = new double[2 * inverseSize - 2];    
    dzWeight1 = new double[inverseSize];  
    dzWeight2 = new double[inverseSize];  
    dzWeight3 = new double[inverseSize]; 
    
    int count = 1;
    sort(inverseVector.begin(), inverseVector.end());
    inverse[0] = inverseVector[0];
    doubleInverse[0] = inverse[0];
    double delta = 10e-10;
    for (int j = 1; j < inverseSize; j++) {
        inverse[j] = inverseVector[j];
        doubleInverse[count] = (inverse[j-1] + inverse[j])/2.0;
        doubleInverse[count+1] = inverse[j];
        dz[count-1] = doubleInverse[count] - doubleInverse[count-1];
        dz[count] = doubleInverse[count+1] - doubleInverse[count];
        if (dz[count] == 0) dz[count] = delta;
        if (dz[count - 1] == 0) dz[count - 1] = delta;
        double hph = dz[count] + dz[count - 1];
        dzWeight1[j] = (pow(dz[count], 3) + pow(dz[count - 1], 3) + 3 * dz[count] * dz[count - 1] * hph ) / (6 * dz[count] * dz[count - 1]);
        dzWeight2[j] = (2 * pow(dz[count - 1], 3) - pow(dz[count], 3) + 3 * dz[count] * pow(dz[count - 1], 2)) / (6 * dz[count - 1] * hph);
        dzWeight3[j] = (2 * pow(dz[count], 3) - pow(dz[count - 1], 3) + 3 * dz[count - 1] * pow(dz[count], 2)) / (6 * dz[count] * hph);         
        
        count += 2;
    }        
    nPointsAdjust = inverseSize;        
}
    
