
#include <vector>
#include "callPDF.h"
#include "Variable.h"
#include "JointProbability.h"

extern "C" { 
    void estimatePDFmv(double *sampleData, int *sampleLength, int *variableLength, int *dimension, int *debug,
            double *x, double *pdf){ 
        
        int nSamples = sampleLength[0];
        int nVariables = variableLength[0];
        int matrixSize = nSamples * nVariables;
        int pdfSize = dimension[0];
        
        vector <Variable> variables;
        variables.reserve(nVariables);
        InputParameters input;
        input.out.debug = debug[0];    
        
        int iVariable = 0;
        vector <double> samples;
        samples.reserve(nSamples);
        int vSample = 0;
        for (int i = 0; i < matrixSize; i++) {
            samples.push_back(sampleData[i]);
            if (++vSample == nSamples) {
                ostringstream vString; 
                vString << iVariable++; 
                Variable variable = Variable(input, samples, vString.str(), false);
                variable.out.debug = debug[0];
                variables.push_back(variable);
                samples.clear();
                vSample = 0;
            }
        }           
   
        JointProbability jp = JointProbability(variables, nSamples, pdfSize);
        jp.out.debug = debug[0];
        jp.calculate();
    
        vector <double> Vpdf = jp.getJP();    
        vector <double> Vx = jp.getRange();     
        
        for (unsigned i = 0; i < Vx.size(); i++) {
            x[i] = Vx[i];
        }
        for (unsigned i = 0; i < Vpdf.size(); i++) {
            pdf[i] = Vpdf[i];
        }
        return;    
    }   
} 
