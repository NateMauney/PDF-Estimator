/* 
 * File:   callPDF.cpp
 * Author: jenny
 * 
 * Created on September 18, 2018, 1:20 PM
 */

#include "callPDF.h"

callPDF::callPDF() {
}

callPDF::callPDF(const callPDF& orig) {
}

callPDF::~callPDF() {
}

void callPDF::makeCall(double * sampleData, int sampleLength, double * estimationPoints, int estimationLength, int isSpecifyPoints, double low, double high, int isLow, int isHigh, double target, int points, int lagrangeMin, int lagrangeMax, int outlierCutoff, int debug, int smooth) {
    InputParameters input;
    out.debug = debug;
    input.out.debug = debug;    
    input.outlierCutoff = outlierCutoff;
    input.smooth = smooth;
    input.minLagrange = lagrangeMin;
    input.maxLagrange = lagrangeMax;
    input.SURDTarget = target;
    input.adaptive = false;
    if (isSpecifyPoints) {
        input.estimatePoints = true;
        vector <double> estimationData;
        for (int i = 0; i < estimationLength; i++) {
            estimationData.push_back(estimationPoints[i]);
        }
        input.setEstimationPoints(estimationData);
    }
    if (isLow) {
        input.lowerBoundSpecified = true;
        input.lowerBound = low;
    } else {
         input.lowerBoundSpecified = false;
    }
    if (isHigh) {
        input.upperBoundSpecified = true;
        input.upperBound = high;
    } else {
         input.upperBoundSpecified = false;
    }
    
    input.integrationPoints = points - 1;
   
    MinimizeScore minimumPDF = MinimizeScore();
    minimumPDF.out.debug = debug;
    InputData data = InputData(input);  
    data.out.debug = debug;
    vector <double> inputData;
    for (int i = 0; i < sampleLength; i++) {
        inputData.push_back(sampleData[i]);
    }
    input.writeHeader = false;
    input.writeFile = false;
    data.setData(inputData);     
    if (data.processData()) {             
        WriteResults write;       
        write.out.debug = debug;
        solutionFailed = minimumPDF.minimize(input, data);
        solutionThreshold = minimumPDF.bestThreshold;
        write.createSolution(input, data, minimumPDF);  
        this->N = data.N;
        write.createQQ(minimumPDF.bestRandom, data.N);   
        Vsqr = write.SQR;
        Vcdf = write.CDF;
        Vpdf = write.PDF;
        VpdfPoints = write.PDFPoints;
        Vx   = write.x;
        Vlagrange = write.L;
        Vr = write.R;
    } else {
        solutionFailed = true;
    }
    
}

