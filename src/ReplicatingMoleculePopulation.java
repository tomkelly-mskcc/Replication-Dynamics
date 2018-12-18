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
import java.util.Arrays;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 *A population of replicating molecules with the same preRC probability distribution and the same simulation parameters
 * @author tkelly
 */
public class ReplicatingMoleculePopulation {
    private final List<ReplicatingMolecule2> replicatingMolecules; //The list of replicating molelcules

    /**
     * Constructor gets probability distribution for current sequence, sets up a list of replicating molecules,
     * each of which has a pre-initiated molecule with a fixed number of initiators distributed according to
     * the probability distribution and the same set of replication parameters.  Replication of each pre-initiated 
     * molecule is completed by the constructor.  After replication,
     * average right fork frequencies in bins of current molecule can be obtained as well as other averages of replication state
     * variables for the population of molecules
     * @param probabilityDistribution an instance of ProbabilityDistribution2
     * @param params  the parameters controlling replication
     * @param rif1List a list of ref1 intervals for suppression of firing
     */
    public ReplicatingMoleculePopulation(ProbabilityDistribution2 probabilityDistribution, ParameterSet params, List<SeqInterval> rif1List) {
        
        //Make list of replicating molecules and replicate each one
        replicatingMolecules = IntStream.range(0, params.getNumberCells()).
                mapToObj(e-> new ReplicatingMolecule2(new PreInitiatedMolecule(probabilityDistribution, new IRandomNumber()), params, rif1List)).
                collect(Collectors.toList());
    }
    
    /**
     * Returns the average right fork frequency in bins across the molecule averaged over the population of molecules
     * Averaged over all the molecules in the population
     * ORIGINALLY TESTED BY USING PRINT RIGHT FORK FREQUENCIES BELOW TO PRINT ALL RIGHT FORK FREQUENCIES FOR ALL MOLECULES AND THEN AVERAGING BY HAND IN EXCEL
     * Tested again in various ways
     * @return list of AVERAGE right fork frequencies in bins across the molecule
     */
    public List<Double> getPopulationAveRtForkFreqInBins() {
        List<Double> rtForksFrequencyDistribution = new ArrayList<>();
        IntStream.range(0, ParameterSet.NUMBER_BINS).forEach(binIndex -> {
            double ave = replicatingMolecules.stream().mapToDouble(mol -> mol.getMoleculeRtForkFreqInBins().get(binIndex)).average().getAsDouble();
            rtForksFrequencyDistribution.add(ave);
        });
        return rtForksFrequencyDistribution;
    }

    /**
     * Returns the termination frequency in bins over population of replicating molecules
     * @return population average termination frequency in bins
     */
    public List<Double> getPopulationAveTerminationFrequencyInBins () {
        List<Integer> tempTerminationFrequencyDistribution = IntStream.range(0, ParameterSet.NUMBER_BINS).
                mapToObj(e-> 0).collect(Collectors.toList());
        for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
            for (ReplicatingMolecule2 mol : replicatingMolecules) {
                tempTerminationFrequencyDistribution.set(i, tempTerminationFrequencyDistribution.get(i) +
                        + mol.getMoleculeTerminationsInBins().get(i));
            }
        } 
        return tempTerminationFrequencyDistribution.stream().map(i-> (double) i / replicatingMolecules.size()).collect(Collectors.toList());
    }
    
    /**
     * Returns the initiation frequency in bins averaged over the whole population of replicating molecules
     * @return population average initiation frequency in bins
     */
    public List<Double> getPopulationAverageInitiationFrequencyInBins() {
        List<Integer> tempInitiationFrequencyDistribution = IntStream.range(0, ParameterSet.NUMBER_BINS).
                mapToObj(e-> 0).collect(Collectors.toList());
        for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
            for (ReplicatingMolecule2 mol : replicatingMolecules) {
                tempInitiationFrequencyDistribution.set(i, tempInitiationFrequencyDistribution.get(i) +
                        + mol.getMoleculeInitiationsInBins().get(i));
            }
        } 
        return tempInitiationFrequencyDistribution.stream().map(i-> (double) i / replicatingMolecules.size()).collect(Collectors.toList());
    }
    
   
    /**
     * Returns a list of the median time of replication of each bin across the genome
     * Each replicating molecule keeps track time of replication of bin of the size used in pu experiment
     * to a resolution approximating the cycle time (minpercycle)
     * @param binSize  size of bins
     * @param numberBins number of bins
     * @return list of median replication times of each bin
     */
    public List<Double> getPopulationMedianTimeOfReplicationInBins(int binSize, int numberBins) {
       List<Double> populationTimeOfReplicationInBins = new ArrayList<>();
       IntStream.range(0, numberBins).forEach(binIndex -> {
            double[] timeArray = replicatingMolecules.stream().mapToDouble(mol -> mol.getTimeOfMolecluleReplicatonInBins().get(binIndex)).toArray();
            Arrays.sort(timeArray);
            populationTimeOfReplicationInBins.add(timeArray[timeArray.length / 2]);
        });
        return populationTimeOfReplicationInBins;
   }
    
    /**
     * Gets the average fraction replicated of each bin across the genome for a given time
     * Average fraction replicated is the fraction of molecules that have replicated the bin at the given time
     * A bin is replicated if any part of it is replicated
     * @param time  elapsed time of replication
     * @param binSize size of the bin   
     * @param numberBins number of bins
     * @param minPerCycle
     * @return 
     */
    public List<Double> getPopulationAveFxReplicatedInBins(double time, int binSize, int numberBins, double minPerCycle) { 
        //Calculate index of StateRecord corresponding to  time
        int recordIndex = (int) ((time + minPerCycle) / ParameterSet.TIME_RECORDING_INTERVAL);
        
        //find fraction replicated of each bin at given time
        
        //declare and initialize array to hold total fraction replicated over all molecules for each bin
        int[] fxReplicatedArray = new int[numberBins];
        int numberMolecules = replicatingMolecules.size();
        
        //iterate over molecules
        for (int mol = 0; mol < numberMolecules; mol++) {
            int tempRecordIndex = recordIndex;
            //get the recorded states
            List<ReplicatingMolecule2.StateRecord2> stateList = replicatingMolecules.get(mol).getRecordedStates();
            //time is greater than the replication time of this molecule, set tempRecordIndex to last record
            if (tempRecordIndex > stateList.size() - 1) {tempRecordIndex = stateList.size() - 1;}
            
            
            //get the replicatedBinList corresponding of the recod with the index tempRecordIndex
            List<Integer> replicatedBins = stateList.get(tempRecordIndex).getReplicatedBinList();
            //if a bin is in a replicatedBinSegment increment the correpoinding elment of the fraction replicated array
            //the fraction replicated of a bin in a molecule is either 0 or 1
            for (int segmentIndex = 0; segmentIndex < replicatedBins.size(); segmentIndex += 2) {
                for (int binindex = replicatedBins.get(segmentIndex); binindex <= replicatedBins.get(segmentIndex + 1);binindex++) {
                    fxReplicatedArray[binindex]++;
                }
            }
        }
        //compute the average fx replicated by dividing the total by the number of molecules
        List<Double> aveFxReplicatedInBins = new ArrayList<>();
        Arrays.stream(fxReplicatedArray).mapToDouble(n -> (double)n / numberMolecules).forEach(fx ->aveFxReplicatedInBins.add(fx));
        return aveFxReplicatedInBins;      
   }
    
    /**
     * Returns an array containing the times of replication of a given bin for all the replicating molecules
     * @param binIndex the index of the bin
     * @return the array of replication times of the bin
     */
    double[] getPopulationTimesOfReplicationOfBin (int binIndex) {
        return replicatingMolecules.stream().mapToDouble(mol->mol.getTimeOfReplicationOfBin(binIndex)).toArray();
    }
           
   /**
     * Returns a list of elapsed times for the population of molecules
     * @return list of elapsed times
     */
    public List<Double> getElapsedTimeList() {
        return replicatingMolecules.stream().map(e-> e.getElapsedTime()).collect(Collectors.toList());
    }
   
    /**
     * properties of replication states
     */
    public enum StateProperty 
    {TIME,FX_REPLICATED,INITIATIONS, FORKS, PASSIVES, POTENTIALS,CLOSURES;}
    
   /**
    * gets the average value of a state variable for the population of replicating molecules for each time point
    * @param getStateVariable a function that takes a StateRecord and returns a double valued state variable
    * @return the average of the state variable returned by a function
    */ 
   List<Double> getAverageStateVariableList(ToDoubleFunction<ReplicatingMolecule2.StateRecord2> getStateVariable)  {
       //Initialize the average statevariable list - list inludes timepoints during replciation in order
        List<Double> averageStateVariableList = new ArrayList<>();
        //get maximum number of states for any molecule
        int maxNumberStates = replicatingMolecules.stream().map(e-> e.getRecordedStates().size()).max((a,b)-> a-b).get(); //This could be done once                
        //Iterate over states
        for (int i = 0; i < maxNumberStates; i++) {
            double total = 0;
            for (int j = 0; j < replicatingMolecules.size(); j++) {
                int maxIndexStatesJthMolecule = replicatingMolecules.get(j).getRecordedStates().size() - 1;
                total += getStateVariable.applyAsDouble(replicatingMolecules.get(j).getRecordedStates().get(i > maxIndexStatesJthMolecule ? maxIndexStatesJthMolecule : i));
            }
            averageStateVariableList.add( total / replicatingMolecules.size());
        }
        return averageStateVariableList;
   }
       
    /**
     * Prints the averages for state variables for population of molecules at each time point sampled during replication
     * checked against more detailed display by PrintProperties
     */
    void printAverageStateVariables() {
        List<Double> averageTimeList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getTime);
        List<Double> averageInitiationsList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberInitiations);
        List<Double> averageForksList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberForks);
        List<Double> averagePotentialsList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberPotentials);
        List<Double> averageTerminationsList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberTerminations);
        List<Double> averagePassivesList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberPassives);
        List<Double> averageNucsList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNucleotidesReplicated);
        List<Double> averageClosuresList = getAverageStateVariableList(ReplicatingMolecule2.StateRecord2::getNumberClosures);
        List<Double> averageFxReplicatedList = averageNucsList.stream().map(e-> e / ParameterSet.SEQUENCE_LENGTH).collect(Collectors.toList());
        System.out.println("Time" + "\t" + "Initiations" + "\t" + "Forks" + "\t" + "Potentials" + "\t" + "Passives" + "\t" + "FxReplicated"
                + "\t" + "Terminations" + "\t" + "Closures");
        for (int i = 0; i < averageTimeList.size(); i++) {
            System.out.println(averageTimeList.get(i) + "\t" + averageInitiationsList.get(i) + "\t" + averageForksList.get(i) + "\t" +
                    averagePotentialsList.get(i) + "\t" + averagePassivesList.get(i) + "\t" + averageFxReplicatedList.get(i) + "\t" +
                    averageTerminationsList.get(i) + "\t" + averageClosuresList.get(i));
        }
        
       
        
    }
}