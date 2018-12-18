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
import java.util.List;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import static java.util.stream.Collectors.summingDouble;
import java.util.stream.IntStream;

/**
 * A population of replicating molecules with same set of parameters
 * For Combing experiments
 * @author tkelly
 */
public class ReplicatingMoleculePopulationForCombingExpt {
        private final List<ReplicatingMoleculeForCombingExpt> replicatingMolecules; //List of repllicating molecules in simulation
        private final List<CombedMolecule> sortedCombedMoleculeList; //list of observed combed molecules sorted by fraction relicated
    
        /**
         * Creates a list of instances of ReplicatingMolecule
         * @param probabilityDistribution the probability distribution for pre-RCs
         * @param params  the set of parameters for replication
         * @param sortedMolList list of observed combed molecules sorted by fraction replicated
         */
        public ReplicatingMoleculePopulationForCombingExpt(ProbabilityDistribution2 probabilityDistribution, 
            ParameterSet params, List<CombedMolecule> sortedMolList) {
            
            this.sortedCombedMoleculeList = sortedMolList;
            
            //Make list of replicating molecules and replicate each one
            replicatingMolecules = IntStream.range(0, params.getNumberCells()).
                    mapToObj(e-> new ReplicatingMoleculeForCombingExpt(new PreInitiatedMolecule(probabilityDistribution, 
                            new IRandomNumber()), params, sortedMolList)).
                    collect(Collectors.toList());
        
    }
    
    /**
    * gets the average values of a state variable for the population of replicating molecules for each time point
    * for printing state variables
    * @param getStateVariable a function that takes a StateRecord and returns a double valued state variable
    * @return the average of the state variable returned by a function
    */ 
   List<Double> getAverageStateVariableList(ToDoubleFunction<ReplicatingMoleculeForCombingExpt.StateRecord> getStateVariable)  {
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
    public void printAverageStateVariables() {
        List<Double> averageTimeList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getTime);
        List<Double> averageInitiationsList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberInitiations);
        List<Double> averageForksList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberForks);
        List<Double> averagePotentialsList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberPotentials);
        List<Double> averageTerminationsList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberTerminations);
        List<Double> averagePassivesList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberPassives);
        List<Double> averageNucsList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNucleotidesReplicated);
        List<Double> averageClosuresList = getAverageStateVariableList(ReplicatingMoleculeForCombingExpt.StateRecord::getNumberClosures);
        List<Double> averageFxReplicatedList = averageNucsList.stream().map(e-> e / replicatingMolecules.get(0).getSeqLength()).collect(Collectors.toList());
        System.out.println("Time" + "\t" + "Initiations" + "\t" + "Forks" + "\t" + "Potentials" + "\t" + "Passives" + "\t" + "FxReplicated"
                + "\t" + "Terminations" + "\t" + "Closures");
        for (int i = 0; i < averageTimeList.size(); i++) {
            System.out.println(averageTimeList.get(i) + "\t" + averageInitiationsList.get(i) + "\t" + averageForksList.get(i) + "\t" +
                    averagePotentialsList.get(i) + "\t" + averagePassivesList.get(i) + "\t" + averageFxReplicatedList.get(i) + "\t" +
                    averageTerminationsList.get(i) + "\t" + averageClosuresList.get(i));
        }    
    }
    
    /**
     * Returns a list of elapsed times for the population of molecules
     * @return list of elapsed times
     */
    public List<Double> getElapsedTimeList() {
        return replicatingMolecules.stream().map(e-> e.getElapsedTime()).collect(Collectors.toList());
    }
    
    
    /**
     * get the averages over the cell population of variables corresponding to the states of the combed molecules (specified by fraction Replicated of the combed molecule)
     * Note that the state recored from the simulation also include zero and 2.0 fraction replicated so we need to use only those from index 1 to index (number records -2)
     * @return a list of the population averages
     */
    public List<PopulationAveragesCombingSimulation> getAverageStateVariablesForCombedMolecules()  {
        //make new averages record
        List<PopulationAveragesCombingSimulation> averagesRecordList = new ArrayList<>();  //need average record for each state (combed molecule state at fx replicated) ca 160 of them
        int numberMols = replicatingMolecules.size(); //the number of replicating molecules in the cell population - for averaging
        int numberStates =  sortedCombedMoleculeList.size(); // equals the number of combed molecules

        //Iterate over states - i.e.  i is index of state of fxReplicated of each combed molecule
        for (int i = 0; i < numberStates; i++) { 
            
            //initialize totals
            int totalInitiations = 0;
            int totalTerminations = 0;
            int totalPassives = 0;
            int totalForks = 0;
            int totalAdjustedForks = 0;
            int totalClosures = 0;
            double totalTime = 0;
            double totalFxReplicated = 0;
            double totalAdjustedFxReplicated = 0;
            List<Integer> aggregatedIodList = new ArrayList<>();
            List<Integer> adjustedAggregatedIodList = new ArrayList<>();
            
            //now get averages over cell population for each state correspoinding to the fx replicated of a combed molelcule.  j is index over cells in population
            for (int j = 0; j < replicatingMolecules.size(); j++) {
                List<ReplicatingMoleculeForCombingExpt.StateRecord> recordedStates = replicatingMolecules.get(j).getRecordedStates();
                if (recordedStates.size() != numberStates + 2) {
                    throw new IndexOutOfBoundsException("recorded states for combed molecule not equal to number of combed moleucles plus two");
                }
                ReplicatingMoleculeForCombingExpt.StateRecord combedMoleculeRecordedState = recordedStates.get(i + 1);
                totalInitiations += combedMoleculeRecordedState.getNumberInitiations();
                totalTerminations += combedMoleculeRecordedState.getNumberTerminations();
                totalPassives += combedMoleculeRecordedState.getNumberPassives();
                totalForks += combedMoleculeRecordedState.getNumberForks();
                totalClosures += combedMoleculeRecordedState.getNumberClosures();
                totalTime += combedMoleculeRecordedState.getTime();
                totalFxReplicated += combedMoleculeRecordedState.getFractionReplicated();
                
                aggregatedIodList.addAll(getInterCentroidDistances(combedMoleculeRecordedState.getSegmentList()));
            }
            double weight = sortedCombedMoleculeList.get(i).getMoleculeProps().getMolLength() / ParameterSet.SEQUENCE_LENGTH;
            averagesRecordList.add(new PopulationAveragesCombingSimulation( (double) totalInitiations /  numberMols, (double) totalTerminations / numberMols, (double) totalPassives / numberMols, 
                    (double) totalForks / numberMols, (double) totalClosures / numberMols, totalTime /numberMols, totalFxReplicated / numberMols, weight, aggregatedIodList, sortedCombedMoleculeList.get(i),
                    (double) totalAdjustedForks / numberMols, totalAdjustedFxReplicated / numberMols, adjustedAggregatedIodList));
        }
        return averagesRecordList;
   }
    
   
    /**
     * returns list of inter-centroid distances from segment list of in state record of replicating molecule 
     * @param segmentList the list of replicated segments from a state record
     */
    List<Integer> getInterCentroidDistances(List<Integer> segmentlist) {
        //initialize segment list straight or ajusted
        List<Integer> interCentroidList = new ArrayList<>();
        
        //test invarient - even n umber of elements in segment list
        if (segmentlist.size() % 2 != 0) {throw new IllegalArgumentException("segment list does not have even number element");}
        
        //initilaize intermediate variables
        int startIndex = 0;
        int numberSegmentBoundaries = segmentlist.size();
        
        //set start to remove starting terminal segment
        if (segmentlist.get(0) == 0) {numberSegmentBoundaries -= 2; startIndex = 2;}
        
        //set end to remove ending terminal segment
        if (segmentlist.get(segmentlist.size() - 1) == ParameterSet.SEQUENCE_LENGTH - 1) {numberSegmentBoundaries -= 2;}
        
        //return empty list if remaining elements less than two segments
        if (numberSegmentBoundaries < 4) {return interCentroidList;}
        
        //for remianing segments get intercentroids
        for (int i = startIndex; i < startIndex + numberSegmentBoundaries - 3; i += 2) {
            interCentroidList.add((segmentlist.get(i + 3) +  segmentlist.get(i + 2)) / 2 - (segmentlist.get(i + 1) + segmentlist.get(i)) / 2);
        }
        return interCentroidList;
    }    
        
    
    /**
     * The complementary CDFs observed and predicted by simulation
     */
    public class ComplementaryCdfs {
        private final List<Double> observedComplementaryCdf;
        private final  List<Double> predictedComplementaryCdf;

        /**
         * Convenience class to hold complementary CDFs of observed and predicted combed molecules with same fraction replicated
         * @param observedComplementaryCdf
         * @param predictedComplementaryCdf 
         */
        public ComplementaryCdfs(List<Double> observedComplementaryCdf, List<Double> predictedComplementaryCdf) {
            this.observedComplementaryCdf = observedComplementaryCdf;
            this.predictedComplementaryCdf = predictedComplementaryCdf;
        }

        /**
         * Returns complementary CDF for observed combed molecules
         * @return complementary CDF for observed molecules
         */
        List<Double> getObservedComplementaryCdf() {
            return observedComplementaryCdf;
        }

        /**
         * Returns complementary CDF for molecules predicted by simulation with same fraction replicated as observed
         * @return complementary CDF of predicted
         */
        List<Double> getPredictedComplementaryCdf() {
            return predictedComplementaryCdf;
        }
    }
    
    
    public ComplementaryCdfs getPredictedVsObservedComplementaryCdfOfIods(Predicate<PopulationAveragesCombingSimulation> fxFilter) {
        //get list of average state variables and filter it 
        List<PopulationAveragesCombingSimulation> averages = getAverageStateVariablesForCombedMolecules();
        List<PopulationAveragesCombingSimulation> filteredAverages = averages.stream().filter(fxFilter::test).collect(Collectors.toList());
        
        //get totalObsIods weight for all simulated moleucles in the filtered population and generate list of weighted IOds for simulated moleucules
        double totalWeight = filteredAverages.stream().collect(summingDouble(e-> e.getWeight() * e.getAggregateIodList().size()));
        List<PopulationAveragesCombingSimulation.WeightedIod> weightedPredictedIodList = filteredAverages.stream().map(e-> e.getWeightedAggregateIodList()).
                collect(ArrayList::new, ArrayList::addAll, ArrayList::addAll);
        
        //Make list of observed iods and get totalObsIods number of them
        List<Double> observedIodList = filteredAverages.stream().map(ave-> ave.getReferenceMolecule().getMoleculeProps().getMolIodLengthList()).
                collect(ArrayList::new, ArrayList::addAll, ArrayList::addAll);
        
        //make cCDFs
        List<Double> observedCcdf = new ArrayList<>();
        List<Double> predictedCcdf = new ArrayList<>();
        
        for (int i = 0; i < ParameterSet.NUMBER_INTERVALS_FOR_CCDF; i++) {
            int totalObsIods = 0;
            for (double len : observedIodList) {
                if (len >= i * 10) {totalObsIods++;}
            }
            observedCcdf.add((double) totalObsIods / observedIodList.size());
            double totalPredIods = 0.0;
            for (PopulationAveragesCombingSimulation.WeightedIod weightedIod : weightedPredictedIodList) {
                if (weightedIod.iod >= i * 10000) {totalPredIods += weightedIod.weight;}
            }
            predictedCcdf.add(totalPredIods / totalWeight);
        }        
            return new ComplementaryCdfs(observedCcdf, predictedCcdf);
    }
}

