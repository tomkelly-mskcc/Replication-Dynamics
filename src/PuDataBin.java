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

/**
 * A collection of data in a bin in a PU experiment and methods to return functions of the data
 * immutable
 * @author tkelly
 */
public class PuDataBin {
    
    private final double normDC;  //normalized to total raw counts of delta on watson
    private final double normDW; //normalized to total raw counts of delta on crick
    private final double normEC; //normalized to total raw counts of epsilon on crick
    private final double normEW; //normalized to total raw counts of epsilon on watson
    private final SeqInterval seqInterval; //The seq interval for this bin
    private final int centerOfBin; //position of bin center
    private final double fractionCrick;  //fraction of total counts for this bin on crick strand
    private final double fractionWatson;  //fraction of total counts for this bin on watson strand
    
    /**
     * constructs a collection of data from a bin in a PU Experiment
     * @param normDC normalized DC counts
     * @param normDW normalized DW counts
     * @param normEC normalized EC counts
     * @param normEW normalized EW counts
     */
    PuDataBin(SeqInterval interval, double normDC, double normDW, double normEC, double normEW, double totalCrickCounts, double totalWatsonCounts) {
        this.normDC = normDC;
        this.normDW = normDW;
        this.normEC = normEC;
        this.normEW = normEW;
        seqInterval = new SeqInterval(interval.start(), interval.end());  //Defensive copy
        centerOfBin = seqInterval.center();
        this.fractionCrick = totalCrickCounts / (totalCrickCounts + totalWatsonCounts);
        this.fractionWatson = totalWatsonCounts / (totalWatsonCounts + totalCrickCounts);
    }
    
    //methods to return normalized counts and various ratios of them
    double getNormDC () {return normDC;}
    double getNormDW () {return normDW;}
    double getNormEC () {return normEC;}
    double getNormEW () {return normEW;}
    double getRatioDC() {return normDC/(normDC + normEC);}  //Crick right
    double getRatioDW() {return normDW/(normDW + normEW);}  //Watson right
    double getRatioEC() {return normEC/(normEC + normDC);}  //Crick left
    double getRatioEW() {return normEW/(normEW + normDW);}  //watson left
    
    /**
     * Average ratios representing average fraction rightward moving forks in a bin
     * Two averages (Watson and Crick strands) are weighted according to total number of counts 
     * on each strand
     * @return average fraction rightward moving forks
     */
    double getRatioDCEW() {  
        double weightedDcRatio;
        double weightedEwRatio;
        if ((normDC + normEC) == 0) {weightedDcRatio = Double.NaN;}  //Deal with zero denominators
        else {weightedDcRatio = normDC * fractionCrick / (normDC + normEC);}
        if (normEW + normDW == 0) {weightedEwRatio = Double.NaN;}
        else {weightedEwRatio = normEW * fractionWatson/(normEW + normDW);}
        return weightedDcRatio + weightedEwRatio;
    }
    /**
     * Average ratios representing average fraction leftward moving forks in a bin
     * Two averages (Watson and Crick strands) are weighted according to total number of counts 
     * on each strand
     * @return average fraction leftward moving forks (1 - fraction rightward moving forks)
     */
    
    double getRatioDWEC() {
        double weightedDwRatio;
        double weightedEcRatio;
        if ((normDC + normEC) == 0) {weightedDwRatio = Double.NaN;}
        else {weightedDwRatio = normDW * fractionWatson / (normDW + normEW);}
        if (normEW + normDW == 0) {weightedEcRatio = Double.NaN;}
        else {weightedEcRatio = normEC * fractionCrick / (normEC + normDC);}
        return weightedDwRatio + weightedEcRatio;
    }
    
    /**
     * Returns the SeqInterval for this bin
     * @return the SeqInverval for this bin
     */
    SeqInterval getSeqInterval() {return seqInterval;}
    
    /**
     * Returns center of bin
     * @return center position of bin
     */
    int getCenterOfBin() {return centerOfBin;}
    
}
