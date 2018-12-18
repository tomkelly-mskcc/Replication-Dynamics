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

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

/**
 * Main routine for simulating DNA replication to compare with PuSeq dataset
 * @author tkelly
 */
public class PuExperimentSimulation {
    
    public enum PuTasks   { PRINT_SEQUENCE,
                            PRINT_PARAMETER_SET,
                            PRINT_CDF,
                            PRINT_MEAN_SQUARED_DEVIATION_RIGHT_FORK_FREQUENCY,
                            PRINT_PREDICTED_AND_OBSERVED_RIGHT_FORK_FREQUENCIES,
                            PRINT_INITIATION_FREQUENCIES,
                            PRINT_TERMINATION_FREQUENCIES,
                            PRINT_MEDIAN_REPLICATION_TIMES_FOR_BINS_IN_GENOME,
                            PRINT_SPHASE_DURATION_FOR_MOLECULES_IN_POPULATION,
                            PRINT_FRACTION_REPLICATED_OF_BINS_AT_GIVEN_TIME,
                            PRINT_REPLICATION_TIME_OF_BIN_FOR_MOLECULES_IN_POPULATION,
                            PRINT_RF_DELTA_LIST
                        }
    /**************************
    //SET TASK FOR MAIN ROUTINE
    **************************/
    static PuTasks task = PuTasks.PRINT_RF_DELTA_LIST;
    static double timeOfReplication = 10; //For routine that gives fraction replicated at a given time - this is given time
    static int binNumber = 10587; //Bin number for routine to get distribution of replication times of bin over molecule population
    
    private static PuExperiment pu; //Pu-seq experiment dataset
    private static SequenceDNA seq; //the DNa sequence of the chromosome/genome
    private static SequenceDNAAnnotationGenbank annotation;  //annotation file for chromosome/genome
    private static List<Transcript> transcriptList;  //a list of transcripts in chrmomosome/genome
    private static AtModelComparator modelComparator; //compares observed and predicted right fork frequencies
    private static ProbabilityDistribution2 probDistribution; //the probability distribution for pre-RCs
    private static ReplicatingMoleculePopulation population; // the population of replicating molecules with a given set of parameters

    public static void main(String[] args) throws IOException {
        
        // set parameter set list
        ParameterSetListMaker maker = new ParameterSetListMaker(DefaultParameterRecipes.getDEFAULT_PARAMETER_RECIPES());
        List<ParameterSet> psList = maker.getParameterSetList();
        
        //Get sequenceDNA a singleton chromsomome II currently
        Path p1 = Paths.get(ParameterSet.FA_SEQ_FILE);
        seq = new SequenceDNA(p1, new SeqInterval(ParameterSet.SEQUENCE_START, ParameterSet.SEQUENCE_END));
        
        //Get PUExpt chromosome II - a singleton
        Path p = Paths.get(ParameterSet.BIN_COUNT_FILE);
        pu = new PuExperiment(p, new SeqBin(ParameterSet.PU_BIN_SIZE,0));  //the pu experiment under study - also sets bin size for analysis
    
        //Get genbank annotation file of chromosome II - a singleton
        Path p2 = Paths.get(ParameterSet.CHROMOSOME_II_GENBANK);
        annotation = new SequenceDNAAnnotationGenbank(p2);
      
        //get transcript list - a singleton
        transcriptList = annotation.getTranscriptsByType(type-> type == Transcript.TRANSCRIPT_TYPE.MRNA || type == Transcript.TRANSCRIPT_TYPE.SNRNA ||
                type == Transcript.TRANSCRIPT_TYPE.RRNA || type == Transcript.TRANSCRIPT_TYPE.SNORNA || type == Transcript.TRANSCRIPT_TYPE.TRNA);
        
        
        //RUN PARAMETER SETS FOR DYNAMIC SIMULATION
            
        //Iterate over all parameter sets
        for (ParameterSet paramSet : psList) {                        

            //Create probability distribution based on transcript exclusion and AT content
            probDistribution = new ProbabilityDistribution2(seq.atFunction(paramSet.getInitiatorSiteLength()), paramSet, 
                    annotation.getTranscriptBlocks(transcriptList));

            //Replicate molecules
            population = new ReplicatingMoleculePopulation(probDistribution, paramSet);

           //Create new model comparator to compare pu data with predicted right fork frequency distribution from population of replicating molecules
            modelComparator = new AtModelComparator(pu, ParameterSet.SEQUENCE_START / ParameterSet.PU_BIN_SIZE, population.getPopulationAveRtForkFreqInBins());
            
            //print list of interTranscripts with probabilities
                /*System.out.println("index\tstart\tend\tAT Content\tRif1 index\tprob");    
                InterTranscriptSites sites = new InterTranscriptSites(seq, annotation.getTranscriptBlocks(transcriptList), rif1List, probDistribution);
                sites.getInterTranscriptList().forEach(System.out::println);
                System.out.println("");*/

            switch (task) {
                case PRINT_PARAMETER_SET:
                    System.out.println(paramSet); 
                    break;
                case PRINT_SEQUENCE:
                    System.out.println(seq.intervalToString(new SeqInterval(0, seq.getLength() - 1)));
                    break;
                case PRINT_CDF:
                    System.out.println(paramSet); 
                    System.out.println("MIDPOINT\tCDF");
                    probDistribution.printCDf(ParameterSet.PU_BIN_SIZE);
                    break;
                case PRINT_MEAN_SQUARED_DEVIATION_RIGHT_FORK_FREQUENCY:
                    System.out.println("MSD\tNumPreRCs\tExpConst\tAttenuator\tInitLength\tMaxRate\tStdDevPredicted\tStdDevObserved");
                    System.out.println(modelComparator.getMeanSquaredDifference() + "\t" + paramSet.getNumberInitiators() + "\t" + paramSet.getExponentialPowerFactor() + "\t" + 
                    paramSet.getAttenuationFactor() + "\t" + paramSet.getInitiatorSiteLength() 
                    + "\t" + paramSet.getMaxFiringProbabilityPerMin() + "\t" + modelComparator.getPredictedStdDeviation() + "\t" + modelComparator.getObservedStdDeviation());
                    break;
                case PRINT_PREDICTED_AND_OBSERVED_RIGHT_FORK_FREQUENCIES:
                    printPredictedAndObservedRtForkFrequencies(paramSet);
                    break;
                case PRINT_INITIATION_FREQUENCIES:
                    printInitiationFrequencies(paramSet);
                    break;
                case PRINT_TERMINATION_FREQUENCIES:
                    printTerminationFrequencies(paramSet);
                    break;
                case PRINT_MEDIAN_REPLICATION_TIMES_FOR_BINS_IN_GENOME:
                    System.out.println(paramSet);
                    System.out.println("Position\tMedianTime");
                    List<Double> replicationTimes = population.getPopulationMedianTimeOfReplicationInBins(ParameterSet.PU_BIN_SIZE, ParameterSet.NUMBER_BINS);
                    for (int i = 0; i < replicationTimes.size(); i++) {
                        System.out.println(i * ParameterSet.PU_BIN_SIZE + ParameterSet.PU_BIN_SIZE / 2 + "\t" + replicationTimes.get(i));
                    }
                    break;
                case PRINT_SPHASE_DURATION_FOR_MOLECULES_IN_POPULATION:
                    System.out.println(paramSet);
                    population.getElapsedTimeList().stream().forEach(f-> System.out.println(f));
                    break;
                case PRINT_FRACTION_REPLICATED_OF_BINS_AT_GIVEN_TIME:
                    System.out.println(paramSet);
                    System.out.println("ELAPSED TIME: " + timeOfReplication + " min");
                    System.out.println("Position\tFxReplicated");
                    List<Double> fxReplicated = population.getPopulationAveFxReplicatedInBins(timeOfReplication, ParameterSet.PU_BIN_SIZE, ParameterSet.NUMBER_BINS, paramSet.getMinPerCycle());
                    for (int i = 0; i < fxReplicated.size(); i++) {
                        System.out.println(i * ParameterSet.PU_BIN_SIZE + ParameterSet.PU_BIN_SIZE / 2 + "\t" + fxReplicated.get(i));
                    }
                    break;
                case PRINT_REPLICATION_TIME_OF_BIN_FOR_MOLECULES_IN_POPULATION:
                    System.out.println(paramSet);
                    System.out.println("BIN NUMBER: " + binNumber);
                    for (int i = 0; i < paramSet.getNumberCells(); i++) {
                        System.out.println(population.getPopulationTimesOfReplicationOfBin(binNumber)[i]);
                    }
                    break;
                case PRINT_RF_DELTA_LIST:
                    PrintDeltaAtContentAndNumberTranscribedPositionsOfBins();
            }
        }    
    }
    
    /**Output predicted and observed right fork frequencies
     * @param pset the current parameter set
     */
     public static void printPredictedAndObservedRtForkFrequencies(ParameterSet pset) {
        System.out.println(pset);
        modelComparator.printRtForkFrequencies();
    }
    
    /**
     * Prints termination frequencies in bins
     * @param pset parameter set
     */
    public static void printTerminationFrequencies (ParameterSet pset) {
        System.out.println(pset);
        System.out.println("TERMINATION FREQUENCIES IN BINS");
        System.out.println("midpoint\tfrequency");
        List<Double> popTerminationsList = population.getPopulationAveTerminationFrequencyInBins();
        for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
            System.out.println(i *  ParameterSet.PU_BIN_SIZE + ParameterSet.PU_BIN_SIZE / 2 + "\t" + popTerminationsList.get(i));
        }
    }  
     
    /**
     * Prints Initiation Frequencies
     * @param pset parameter set
     */
    public static void printInitiationFrequencies (ParameterSet pset) {
        System.out.println(pset);
        System.out.println("INITIATION FREQUENCIES IN BINS");
        System.out.println("midpoint\tfrequency");
        List<Double> popInitiationsList = population.getPopulationAverageInitiationFrequencyInBins();
        for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
            System.out.println(i *  ParameterSet.PU_BIN_SIZE + ParameterSet.PU_BIN_SIZE / 2 + "\t" + popInitiationsList.get(i));
        }
    }    
     
     public static void PrintDeltaAtContentAndNumberTranscribedPositionsOfBins() {
         int[] transcribedNucleotidesInBinsList = annotation.getNumberTranscribedNucleotidesInPuBins(transcriptList);
         System.out.println("DELTA LIST");
         System.out.println("Midpoint\tAT Content/tDelta/tNucsInTranscripta");
         for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
             System.out.println(i * ParameterSet.PU_BIN_SIZE + ParameterSet.PU_BIN_SIZE / 2 + "\t" + seq.getAtContentPuBin(i) + "\t" + 
                     pu.getDeltaOfBin(i) + "\t" + transcribedNucleotidesInBinsList[i]);
             
         }
     }
     
     
}
    
       
            
        
        
        
