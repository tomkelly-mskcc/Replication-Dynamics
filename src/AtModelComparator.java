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
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.stream.Collectors;

/**
 *Compares observed and predicted right fork frequency distributions
 * @author tkelly
 */
class AtModelComparator {
    private final List<Double> rtForkDistributionObserved;
    private final List<Double> rtForkDistributionPredicted;
    private final double meanSquaredDifference;
    private final  PuExperiment puExpt;
    

    /**
     * Creates instance of comparator - computes mean squared difference
     * between observed and predicted right fork frequency distributions
     * predicted and observed have to start at same bin in DNA sequence underlying pu experiment
     * @param expt a PuSeq experiment
     * @param puStartingBin index of first bin in PuSeq experiment to be used for comparison
     * @param rtForkPrediction right fork frequency prediction for each seq bin in DNA sequence
     */
    public AtModelComparator(PuExperiment expt, int puStartingBin, List<Double> rtForkPrediction) {
        
        rtForkDistributionPredicted = rtForkPrediction;
        puExpt = expt;
        int numberBins = rtForkDistributionPredicted.size();
        List<Double> tempRtForkObserved =  new ArrayList<>(); 

       List<PuDataBin> puDataBinList = expt.getPuDataArray();
        for (int i = puStartingBin; i < puStartingBin + numberBins; i++) {
            tempRtForkObserved.add(puDataBinList.get(i).getRatioDCEW());
        }
        rtForkDistributionObserved = Collections.unmodifiableList(tempRtForkObserved);
        
        
        double tempSquaredDiff = 0;
        for (int i = 0; i < rtForkDistributionPredicted.size(); i++) {
            double diff = rtForkDistributionPredicted.get(i) - rtForkDistributionObserved.get(i);
            if (! Double.isNaN(diff)) {tempSquaredDiff += (diff * diff);}
        }
        
        meanSquaredDifference = tempSquaredDiff / rtForkDistributionPredicted.size();
        
    }

    /**
     * returns mean squared difference between observed and predicted average right fork frequencies in bins
     * @return mean squared difference between observed and predicted average right fork frequencies in bins
     */
    public double getMeanSquaredDifference() {
        return meanSquaredDifference;
    }
    
    
    /**
     * Prints predicted and observed right fork frequencies in bins across genome - if observed frequency 
     * is NaN de to low raw bin counts leading to division by zero prints a "" for both predicted and observed
     */
    public void printRtForkFrequencies() {
        System.out.println("RIGHT FORK FREQUENCIES");
        System.out.println("Position" + "\t" + "Predicted" + "\t" + "Observed");
        for (int i = 0; i < rtForkDistributionPredicted.size(); i++) {
            double position = puExpt.getPuBinSize() / 2 + puExpt.getPuBinSize() * i;
            if (Double.isNaN(rtForkDistributionObserved.get(i))) {
                System.out.println("" + "\t" + "");
            }
            else {System.out.println(position + "\t" + rtForkDistributionPredicted.get(i) + "\t" + rtForkDistributionObserved.get(i));}
        }
    }
    
    /**
     * Returns standard deviation of the right fork frequencies in bins across genome
     * @return standard deviation predicted right fork frequencies in bin across genome
     */
    double getPredictedStdDeviation() {
        DoubleSummaryStatistics stats = rtForkDistributionPredicted.stream().collect(Collectors.summarizingDouble(e -> e*e));
        return stats.getAverage() - .25;
    }
    
    /**
     * Returns standard deviation of observed right fork frequency - extra
     * steps required to deal with NaN due to low raw counts in bin
     * @return standard deviation observed right fork frequencies in bin across genome
     */
    double getObservedStdDeviation() {
        double sumSquares = 0;
        int numberValues = 0;
        for (int i = 0; i < rtForkDistributionObserved.size(); i++) {
            double temp = rtForkDistributionObserved.get(i);
            if (! Double.isNaN(temp)) {
                sumSquares += temp*temp;
                numberValues++;
            }
        }
        return sumSquares / numberValues - .25;
    }

    /**
     * returns list of observed average rightward fork frequencies for bins
     * @return rtForkDistributionObserved
     */
    List<Double> getRtForkDistributionObserved() {
        return rtForkDistributionObserved;
    }

    /** 
     * returns list of predicted average rightward fork frequencies for bins
     * @return rtForkDistributionPredicted
     */
    List<Double> getRtForkDistributionPredicted() {
        return rtForkDistributionPredicted;
    }
    
    
}
