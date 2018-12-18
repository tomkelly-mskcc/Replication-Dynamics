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
 * A set of parameters for modeling preRC distribution and preRC firing
 * immutable
 * @author tkelly
 */
public class ParameterSet {
    public static final String GENOME_NAME;  //Name of the genome/chromosome to be replicated
    public static final String HOME_DIR;  //The home directory
    public static final String COMBING_DATA_FOLDER; //Folder containing Kaykov and Nurse 2015 dataset reformatted by TKelly
    public final static String FA_SEQ_FILE;  //Seqeunce file of chromosome II in FastA format
    public static final String BIN_COUNT_FILE; //File of bin counts - raw Pu-seq data from Daigaku et al 2015
    public static final String CHROMOSOME_II_GENBANK;  //Genbank file for chromosome II
    public final static String RIF1_INTERVAL_FILE;
    public final static int SEQUENCE_START;  //needs to be on bin boundary
    public final static int SEQUENCE_END;  //needs to be on bin boundary
    public final static int SEQUENCE_LENGTH; //Length of DNA seqeunce of chromosome II
    public final static double TIME_RECORDING_INTERVAL; //time in min interval between recording intermediate results from synthesis engine
    public final static double TIME_RECORDING_START; //Time to start reocrding - i.e. time of first record
    public final static int PU_BIN_SIZE; //size of bin in PuSeq experiment 300bp
    public final static int NUMBER_BINS;  //number of bins in seqeunce
    public final static double HISTOGRAM_BIN_SIZE; //bin size for certain histograms
    public final static int NUMBER_INTERVALS_FOR_CCDF; //number of intervals in complementary CDF for combing experiment
    public final static int DISTANCE_INTERVAL_FOR_CCDF; //interval for plotting inter-centroid distances
    
    //Definitiaitons of these parameters given in DefaultParameterRecipes
    private final int initiatorSiteLength; 
    private final int numberPreRCs; 
    private final int numberCells;
    private final double exponentialCoefficient;
    private final double exponentialPowerFactor;
    private final double maxFiringProbabilityPerMin;
    private final int elongationRate;
    private final double timeConstant;
    private final double minPerCycle;
    private final double attenuationFactor;
    private final double unusedParameter;
    
    
    static {
    GENOME_NAME = "S_POMBE_CHROMOSOME_II";
    HOME_DIR = System.getProperty("user.home");
    FA_SEQ_FILE = HOME_DIR + "/Dropbox/Sync/NetBeansProjects/FilesForTesting/pombe_sequenceII_test.fa";
    COMBING_DATA_FOLDER = HOME_DIR + "/Dropbox/Replication Modeling/Origin Project/csvFiles";
    BIN_COUNT_FILE = HOME_DIR + "/Dropbox/Sync/Replication Modeling/Carr Stuff/2015 Daigaku &Carr Nat S Mol Bio/Chromosome_II_bin_counts.csv";
    CHROMOSOME_II_GENBANK = HOME_DIR + "/Dropbox/Sync/Replication Modeling/Carr Stuff/pombe 294v2.23/Schizosaccharomyces_pombe.ASM294v2.23.II.genbank";
    RIF1_INTERVAL_FILE = HOME_DIR + "/Dropbox/Sync/NetBeansProjects/FilesForTesting/Rif1Intervals.txt";
    SEQUENCE_START = 0;  //needs to be on bin boundary
    SEQUENCE_END = 4539599;  //4539599 this is seqlength - 1 seqLength - seqlength needs to be  integral number of bins 12999900
    SEQUENCE_LENGTH = SEQUENCE_END - SEQUENCE_START + 1;
    TIME_RECORDING_INTERVAL = 0.25;  //interval for recording stgates as fraction replicated - note: MUST BE BIG ENOUGH SO TOTAL FORK MOVEMENT DOESN'T EXCEED IT
    TIME_RECORDING_START =.25;
    PU_BIN_SIZE = 300;
    NUMBER_BINS = SEQUENCE_LENGTH / PU_BIN_SIZE;
    HISTOGRAM_BIN_SIZE = .02;
    NUMBER_INTERVALS_FOR_CCDF = 100;
    DISTANCE_INTERVAL_FOR_CCDF = 10; //in kb
    }
    
    /**
     * //Creates new parameter set for DNA replication
     * @param initSiteLen length of window for initiator site
     * @param numbPreRCs number pre-RCs
     * @param numCells number of cells (molecules replicated)
     * @param exponentialCoefficient deprecated
     * @param exponentialPowerFactor constant for exponential relating AT content to probability
     * @param maxFiringRate maximum rate of preRC firing
     * @param elongRate velocity of forks in nucleotides per minute
     * @param timeConst rate of increase of firing rate 
     * @param cycleTime cycle time of synthesis engine
     * @param attenuator value to set probability in transcription units
     * @param unused for future use
     */
    public ParameterSet(int initSiteLen, int numbPreRCs,int numCells, double exponentialCoefficient, 
            double exponentialPowerFactor, double maxFiringRate, int elongRate, double timeConst, 
            double cycleTime, double attenuator, double unused) {
        this.initiatorSiteLength = initSiteLen;
        this.numberPreRCs = numbPreRCs;
        this.numberCells = numCells;
        this.exponentialCoefficient = exponentialCoefficient;
        this.exponentialPowerFactor = exponentialPowerFactor;
        this.maxFiringProbabilityPerMin = maxFiringRate;
        this.elongationRate = elongRate;
        this.timeConstant = timeConst;
        this.minPerCycle = cycleTime;
        this.attenuationFactor = attenuator;
        this.unusedParameter = unused;
    }

    /**
     * Returns elongation rate in nucleotides per min
     * @return elongation rate
     */
    int getElongationRate() {
        return elongationRate;
    }

    /**
     * Returns window for initiator site
     * @return length of initiator window
     */
    int getInitiatorSiteLength() { 
        return initiatorSiteLength;
    }

    /**
     * Returns number of preRCs
     * @return number of preRCs
     */
    int getNumberInitiators() {
        return numberPreRCs;
    }

    /**
     * number of molecules replicated
     * @return number of replicating molecules with given parameter set
     */
    int getNumberCells() {
        return numberCells;
    }

    /**
     * Returns constant for exponential function relating AT content to probability 
     * @return exponential constant
     */
    double getExponentialPowerFactor() {
        return exponentialPowerFactor;
    }

    /**
     * Returns maximum firing rate
     * @return maximum firing rate
     */
    double getMaxFiringProbabilityPerMin() {
        return maxFiringProbabilityPerMin;
    }

    /**
     * Returns cycle time for synthesis engine
     * @return synthesis cycle time
     */
    double getMinPerCycle() {
        return minPerCycle;
    }

    /**
     * Deprecated
     * @return 
     */
    double getExponentialCoefficient() {
        return exponentialCoefficient;
    }

    /**
     * Returns the rate of increase of firing rate
     * @return rate of increase of firing rate
     */
    double getTimeConstant() {
        return timeConstant;
    }

    /**
     * Returns factor to set probability of preRC in transcription unit
     * @return attenuation factor for transcription unit
     */
    double getAttenuationFactor() {
        return attenuationFactor;
    }

    /**
     * for future use
     * @return 
     */
    double getUnusedParameter() {
        return unusedParameter;
    }

    
    
    /**
     * String giving values of key parameters
     * @return 
     */
    @Override
    public String toString() {
        String str = String.format("REPLICATION PARAMETERS\n");
        str = str.concat(String.format("Genome\t%s\nLength\t%d\nStart\t%d\nEnd\t%d\nPuBinSize\t%d\nNumberBins\t%d\n"
                + "SiteLength\t%d\nNumberPreRCs\t%d\nNumberCells\t%d\nExpCoeff\t%f\nPowerFactor\t%f\n"
                + "FiringProb\t%f\nElongRate\t%d\nTimeCons\t%f\nCycleTime\t%f\nAttenuator%f\n\n",
                GENOME_NAME,SEQUENCE_LENGTH,SEQUENCE_START,SEQUENCE_END,PU_BIN_SIZE,NUMBER_BINS,
                getInitiatorSiteLength(),getNumberInitiators(), getNumberCells(),
                getExponentialCoefficient(), getExponentialPowerFactor(), getMaxFiringProbabilityPerMin(), getElongationRate(), 
                getTimeConstant(), getMinPerCycle(), getAttenuationFactor()));
        return str;
    }
    
    
    
    
}
