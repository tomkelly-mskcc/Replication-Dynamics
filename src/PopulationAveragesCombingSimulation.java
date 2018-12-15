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

import java.util.List;
import java.util.stream.Collectors;

/**
 * The averages of state variables for a population of simulated replicating molecules - each molecule
 * in the population corresponds to a combed molecules observed by Kaykov and Nurse 
 * @author tkelly
 */
public class PopulationAveragesCombingSimulation {
    private final double aveNumberInitiations;
    private final double aveNumberTerminations;
    private final double  aveNumberPassives;
    private final double aveNumberForks;
    private final double adjustedAveNumberForks;
    private final double forksPerMb;
    private final double adjustedForksPerMb;
    private final double aveNumberClosures;
    private final double time;
    private  final double fxReplicated;
    private final double adjustedFxReplicated;
    private final double weight;
    private  final List<Integer> aggregateIodList;
    private final List<Integer> adjustedAgregateIodList;
    private final List<WeightedIod> weightedAggregateIodList;
    private final List<WeightedIod> adjustedWeightedAggregateIodList;
    private final CombedMolecule referenceMolecule;
    private final double referenceMoleculeFxReplicated;
    private final double referenceMoleculeNumberforks;
    private final double referenceMoleculeForksPerMb;
    private final double referenceMoleculeSize;

    /**
     * Creates the instance containing all the population averages for combing experiments
     * This constructor is used for Combing experiments
     * @param numberInitiations average number of initiations
     * @param numberTerminations average number of 
     * @param numberPassives average number of 
     * @param numberForks average number of 
     * @param numberClosures average number of 
     * @param time elapsed time
     * @param fxReplicated fraction replicated
     * @param weight
     * @param aggregatedIodList
     * @param refMol
     * @param adjustedNumberForks
     * @param adjustedFxReplicated
     * @param adjustedAggregatedIodList 
     */
    public PopulationAveragesCombingSimulation(double numberInitiations, double numberTerminations, double numberPassives, double numberForks, double numberClosures, double time, 
            double fxReplicated, double weight, List<Integer> aggregatedIodList, CombedMolecule refMol, double adjustedNumberForks, double adjustedFxReplicated, List<Integer> adjustedAggregatedIodList) {
        this.aveNumberInitiations = numberInitiations;
        this.aveNumberTerminations = numberTerminations;
        this.aveNumberPassives = numberPassives;
        this.aveNumberForks = numberForks;
        this.adjustedAveNumberForks = adjustedNumberForks;
        this.aveNumberClosures = numberClosures;
        this.time = time;
        this.fxReplicated = fxReplicated;
        this.adjustedFxReplicated = adjustedFxReplicated;
        this.weight = weight;  // observed molecule length
        this.aggregateIodList = aggregatedIodList;  //list of iods fro all moleucles
        this.adjustedAgregateIodList = adjustedAggregatedIodList;
        
        weightedAggregateIodList = aggregatedIodList.stream().map(e-> new WeightedIod(e, this.weight)).collect(Collectors.toList());
        adjustedWeightedAggregateIodList = adjustedAggregatedIodList.stream().map(e-> new WeightedIod(e, this.weight)).collect(Collectors.toList());
        
        referenceMolecule = refMol;
        referenceMoleculeFxReplicated = refMol.getMoleculeProps().getFxReplicated();
        referenceMoleculeNumberforks = refMol.getMoleculeProps().getNumberForks();
        referenceMoleculeSize = refMol.getMoleculeProps().getMolLength();
        referenceMoleculeForksPerMb = refMol.getMoleculeProps().getNumberForks() * 1000 / refMol.getMoleculeProps().getMolLength();
        
        forksPerMb = (double) numberForks * 1000000 / ParameterSet.SEQUENCE_LENGTH;
        adjustedForksPerMb = (double) adjustedNumberForks * 1000000 /ParameterSet.SEQUENCE_LENGTH;
    }
    
    
    public class WeightedIod {
        int iod;
        double weight;

        public WeightedIod(int iod, double weight) {
            this.iod = iod;
            this.weight = weight;
        }
    }
    
    
    
    public double getFxReplicated() {
        return fxReplicated;
    }

    public double getNumberClosures() {
        return aveNumberClosures;
    }

    public List<Integer> getAggregateIodList() {
        return aggregateIodList;
    }

    public double getNumberForks() {
        return aveNumberForks;
    }

    public double getNumberInitiations() {
        return aveNumberInitiations;
    }

    public double getNumberPassives() {
        return aveNumberPassives;
    }

    public double getNumberTerminations() {
        return aveNumberTerminations;
    }

    public double getTime() {
        return time;
    }

    public double getWeight() {
        return weight;
    }

    public List<WeightedIod> getWeightedAggregateIodList() {
        return weightedAggregateIodList;
    }

    public CombedMolecule getReferenceMolecule() {
        return referenceMolecule;
    }

    public double getForksPerMb() {
        return forksPerMb;
    }

    public double getReferenceMoleculeForksPerMb() {
        return referenceMoleculeForksPerMb;
    }

    public double getReferenceMoleculeFxReplicated() {
        return referenceMoleculeFxReplicated;
    }

    public double getReferenceMoleculeNumberforks() {
        return referenceMoleculeNumberforks;
    }

    public double getReferenceMoleculeSize() {
        return referenceMoleculeSize;
    }

    public List<Integer> getAdjustedAgregateIodList() {
        return adjustedAgregateIodList;
    }

    public double getAdjustedAveNumberForks() {
        return adjustedAveNumberForks;
    }

    public double getAdjustedForksPerMb() {
        return adjustedForksPerMb;
    }

    public double getAdjustedFxReplicated() {
        return adjustedFxReplicated;
    }

    public List<WeightedIod> getAdjustedWeightedAggregateIodList() {
        return adjustedWeightedAggregateIodList;
    }
    
    
    
        
        

    
}
