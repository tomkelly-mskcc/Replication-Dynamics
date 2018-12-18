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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Immutable Class represents PUSeq Data from Carr PU seq experiments
 * Bins are numbered 0 to PUExperiment.numberBins - data normalized to total counts of each type  No correction for 50:50
 * @author tkelly
 */
class PuExperiment {
    
    private final int MAX_NUMBER_BINS = 16000;
    final static int BIN_SIZE_300 = 300;
    final static int BIN_SIZE_600 = 600;
    final static int BIN_SIZE_900 = 900;
    
    private final int numberBins;
    private final SeqBin puBin; //The seqBin used for PU Experiment (300, 0) for Carr data 
    private final int puBinSize;
    
    //raw bin counts
    private final int[] rawDC = new int[MAX_NUMBER_BINS];
    private final int[] rawDW = new int[MAX_NUMBER_BINS];
    private final int[] rawEC = new int[MAX_NUMBER_BINS];
    private final int[] rawEW = new int[MAX_NUMBER_BINS];
    
    private final List<PuDataBin> puDataArray;  //unmodifianle array
    private List<Double> observedRtForkFrequencyDistribution = null;
    //private final List<Double> puDeltaDWECArray;  //unmodifiable array  DEPRECATED
    
    
    /**
     * Creates instance of PUExperiment from PUSeq data file
     * only complete data bins are included (determined from underlying linked sequence)
     * @param p path to file containing PUSeq data - comma delimited four raw bin counts for DC,DW,EC,EW
     * @param linkedSeq  the DNA sequence linked to this PU experiment
     * @param bin the SeqBin associated with this experiment
     */
    PuExperiment(Path p, SeqBin bin) {
        String rawLine;
        BufferedReader inputStream;
        String[] binCounts; //array of strings giving raw counts for DC [0],DW [1], EC [2], EW [3]
        int binIndex = 0; //initialize bin index to reaading pu data file
        this.puBin = new SeqBin(bin.binSize(), bin.offset());  //Defensive copy
        this.puBinSize = bin.binSize();  //Convenience bin size for pu expt
        
        //totals raw bincnounts for normalization
        int totalDC,totalDW, totalEC, totalEW;
        totalDC = totalDW = totalEC = totalEW = 0;
      
        //Get Path for PU analysis file, open a buffered reader, and read raw counts for DC, DW, EC,EW for each bin
        try {
            inputStream = Files.newBufferedReader(p);
            while((rawLine = inputStream.readLine()) != null) {
                binCounts = rawLine.split(",");
                rawDC[binIndex] = Integer.parseInt(binCounts[0]);
                totalDC += rawDC[binIndex];
                rawDW[binIndex] = Integer.parseInt(binCounts[1]);
                totalDW += rawDW[binIndex];
                rawEC[binIndex] = Integer.parseInt(binCounts[2]);
                totalEC += rawEC[binIndex];
                rawEW[binIndex] = Integer.parseInt(binCounts[3]);
                totalEW += rawEW[binIndex];
                binIndex++;
            }
        } catch(IOException x) {
                System.err.format("Problem opening data file: %s%n", x);                    
        }
        this.numberBins = binIndex;  //set number of bins in pu data

        
        //normalize and store the raw bin counts
        List<PuDataBin> tmpDataArray = new ArrayList<>();
        for (int i = 0; i < numberBins; i++) {
                tmpDataArray.add(new PuDataBin(
                new SeqInterval(i * puBinSize, (i + 1) * puBinSize - 1),
                ((double)rawDC[i])/(totalDC),
                ((double)rawDW[i])/(totalDW),
                ((double)rawEC[i])/(totalEC),
                ((double)rawEW[i])/(totalEW),
                (double)rawDC[i] + (double)rawEC[i],
                (double)rawDW[i] + (double)rawEW[i]));            
       }
        puDataArray = Collections.unmodifiableList(tmpDataArray);
    }
   
    
    /**
     * returns the number of bins in the data set
     * @return number of bins
     */
    int numberBins() {
        return numberBins;
    }
    
    /**
     * Returns the SeqBin for this 
     * @return SeqBin for PUSeqData Object
     */
    SeqBin puSeqBin () {
        return puBin;
    }
    
    /**
     * Returns pu bin size
     * @return size of bin in Pu experiment in bp (300bp for Daigaku Carr Experiment)
     */
    public int getPuBinSize() {
        return puBinSize;
    }
  
    /**
     * Returns the raw bin counts for pol delta on Crick strand
     * @param index of bin 0 to number of bins -1
     * @return 
     */
    int rawDC (int index) {
        checkBinIndex(index);
        return rawDC[index];
    }
    
     /**
     * Returns the raw bin counts for pol delta on Watson strand
     * @param index of bin 0 to number of bins -1
     * @return 
     */
    int rawDW (int index) {
        checkBinIndex(index);
        return rawDW[index];
    }
    
     /**
     * Returns the raw bin counts for pol epsilon on Crick strand
     * @param index of bin 0 to number of bins -1
     * @return 
     */
    int rawEC (int index) {
        checkBinIndex(index);
        return rawEC[index];
    }
    
     /**
     * Returns the raw bin counts for pol epsilon on Watson strand
     * @param index of bin 0 to number of bins -1
     * @return 
     */
    int rawEW (int index) {
        checkBinIndex(index);
        return rawEW[index];
    }
    
    
    /**
     * returns a new ArrayList of PuDataBins containing normalized data for each PuSeq data bin
     * @return ref to the array
     */
    List<PuDataBin> getPuDataArray() {
        return puDataArray;
    }
    
    
    /**
     * returns  ArrayList of the differences between ratioDWEC in adjacent bins
     * @return List of deltas of this experiment, an unmodifiable list
     */
    /*List<Double> getDeltaDWECArray () {
        return puDeltaDWECArray;
    }*/
       
    /** 
     * Throws illegal argument exception if bin index is out of bounds
     * @param index of bin
     */
    private void checkBinIndex(int index) {
        if (index < 0 || index > numberBins - 1) {
            throw new IllegalArgumentException("index to PU data bin is out of bounds 0 to number of bins - 1");
        }
    }
    
    /**
     * Returns observed right fork frequency distribution
     */
    void getRtForkFreqDistributionObserved() {
        observedRtForkFrequencyDistribution =  getPuDataArray().stream().map(e -> e.getRatioDCEW()).collect(Collectors.toList());
    }
   
   /**
    * Returns right fork frequency of a bin
    * @param binIndex the index of bin in the molecule
    * @return right fork frequency of a bin
    */
    double getRtForkFreqOfBin(int binIndex) {
       if (observedRtForkFrequencyDistribution == null) {getRtForkFreqDistributionObserved();}
       return observedRtForkFrequencyDistribution.get(binIndex);
   }
   
   /**
    * returns the RFdelta of a bin
    * @param binIndex index of bin
    * @return  RFdelta of bin
    */
    public double getDeltaOfBin(int binIndex) {
       if (binIndex == 0) {return 0;}
       return getRtForkFreqOfBin(binIndex) - getRtForkFreqOfBin(binIndex - 1);
   }
}
