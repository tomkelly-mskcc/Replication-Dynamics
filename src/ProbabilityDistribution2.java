/* 
 * The MIT License
 *
 * Copyright 2018 Thomas Kelly.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package replicationdynamics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 *Probability distribution for initiator distribution based on exponential function of AT content
 * Dependencies - AtFunction and parameter set (with parameters for exponential function)
 * plus an attenuator function to take into account effects of transcription etc.
 * @author tkelly
 */
public class ProbabilityDistribution2 {
    private final List<Double> cdf;  //cummulative distribution function getter returns unmodifiable form
    private final List<Double> pmf; //probability mass function getter returns unmodifiable form 
    private final ParameterSet parameterSet;  //immutable
    private final int seqLength;

    /**
     * Constructor for probability distribution derived from normalized exponential function of AT 
     * content weighted by attenuator function based on transcription blocks
     * @param atFunction List of number of ATs in initiator sites as a function of position in the genome
     * @param paramSet this parameter set defining parameters for exponential function
     * @param transcriptionBlockList List of transcription blocks for transcription attenuator function
     */
    public ProbabilityDistribution2(List<Integer> atFunction, ParameterSet paramSet, List<SeqInterval> transcriptionBlockList) {
        
        //Initialize fields
        parameterSet = paramSet;
        int initiatorSiteLength = parameterSet.getInitiatorSiteLength();  
        seqLength = atFunction.size() + initiatorSiteLength - 1;  //CHECK
 
        //Get exponential table for initiator site length -             
        if (Double.isInfinite(atFunction.size() * Math.exp(parameterSet.getExponentialPowerFactor()))) {throw new ArithmeticException("kvalue * sequence length is too big - result may be infinite");} //check for overflow
         AtExponentialTable expTable = new AtExponentialTable(initiatorSiteLength, parameterSet.getExponentialPowerFactor());  //make a table to elimate duplication in calculating exponeital functin of AT content

        //get Attenuator Function for weighting
        List<Double> attenuatorFunction = new AttenuatorDistribution(transcriptionBlockList, parameterSet, seqLength).getAttenuatorDistribution();
         
        //create probability mass function from exponential(at function) and attenuator functon
        List<Double> tempPmf = atFunction.stream().map(e-> parameterSet.getExponentialCoefficient() * expTable.getAtExponential(e)).collect(Collectors.toList());
        final double total = tempPmf.stream().mapToDouble(e-> e).sum();
        for (int i = 0; i < tempPmf.size(); i++) {tempPmf.set(i, tempPmf.get(i) * attenuatorFunction.get(i) / total);} //normalize and multiply by attenuator function
        final double total1 = tempPmf.stream().mapToDouble(e-> e).sum();
        pmf = tempPmf.stream().map(e-> e / total1).collect(Collectors.toList()); //renormalize
        //printProbabilityDistributionInBins();
        
        //create cummulative distribution function 
        List<Double> tempCDF1 = new ArrayList<>();
        tempCDF1.add(pmf.get(0));
        for (int i = 1; i < pmf.size(); i++) {tempCDF1.add(tempCDF1.get(i - 1) + pmf.get(i));}
        cdf = tempCDF1;
        //System.out.println(paramSet);
        //printCDf(300);
        }    

    /**
     * constructor for a uniform probability distribution
     * @param seqLen
     * @param paramSet 
     */
    public ProbabilityDistribution2 (int seqLen, ParameterSet paramSet) {
        //Initialize fields
        seqLength = seqLen;
        parameterSet = paramSet;
        int initiatorSiteLength = parameterSet.getInitiatorSiteLength();
        int distributionSize = seqLength - initiatorSiteLength + 1; 
                
        List<Double> tempPmf = new ArrayList<>();
        for (int i = 0; i < distributionSize; i++) {
            tempPmf.add(1.0 / distributionSize);
        }
        pmf = tempPmf;
        List<Double> tempCdf = new ArrayList<>();
        tempCdf.add(pmf.get(0));
        for (int i = 1; i < distributionSize; i++) {
            tempCdf.add(tempCdf.get(i - 1) + pmf.get(i));
        }
        cdf = tempCdf;
 }
    
    /**
    * returns the initiator site length from parameter set
    * @return initiator site length
    */
    int getInitiatorSiteLength() {
        return parameterSet.getInitiatorSiteLength();
    }

    /**
     * Returns CDF probability distribution 
     * @return CDF
     */
    public List<Double> getCdf() {
        return Collections.unmodifiableList(new ArrayList<>(cdf));  //Defensive unmodifiable copy
    }

    /**
     * returns this parameter set
     * @return this parameter set
     */
    ParameterSet getParameterSet() {
        return parameterSet;
    }
    
    /**
     * Returns length of sequence - convenience method carries length forward
     * @return length of sequence
     */
    int getSeqLength() {
        return seqLength;
    }
    
     void printCDf (int binSize) {
        //System.out.println(parameterSet);
        int numberBins = cdf.size() / binSize;
        for (int i = 0; i < numberBins; i++) {
            int offset = binSize * i + binSize / 2;
            System.out.println(i * ParameterSet.PU_BIN_SIZE + 150 + "\t" + cdf.get(offset));
        }
    }
}