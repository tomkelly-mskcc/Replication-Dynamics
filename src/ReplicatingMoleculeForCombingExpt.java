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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A replicating Molecule 
 * for Combing experiment
 * @author tkelly
 */
public class ReplicatingMoleculeForCombingExpt {
    
    //Constants of replication from parameters
    private final List<PotentialInitiationSite2> potentialInitSitesList;  //the list of potential initiation sites encoding states of each site during DNA replication
    private final int seqLength;  //the length of the DNA sequeccne of the replicating molecule
    
    // Variables related to replication progress  
    // just need enough to calculate setup for each synthetic cycle
    private double elapsedTime = 0;  //elapsed replication time - increments are calculated by synthesis engine
    private double fxReplicated = 0; //fraction replicated
    private final DnaSynthesisCycle2 cycle;  //the synthesis cycle for this replicating molecule
    
    //Variables needed to record states and put results in bins as needed
    private final List<ReplicatingMoleculeForCombingExpt.StateRecord> statesAtRegularIntervals;  // list of states at regular intervals - now approximately 0.01 fraction replicated
    private final List<Double> timeOfReplicatonInBins;  //list of time of replication of bins
    private final  int binSize; //size of bin
    
    
    /**
     * Creates a replicating molecule and replicates it
     * Collects a record of molecular states when fraction replicated is the same as an observed combed molecule
     * @param mol the replicating molecule with potential initiation sites (preRC positions)
     * @param parameters the parameter set for replication
     * @param molList the list of observed combed molecules to simulate
     */
    public ReplicatingMoleculeForCombingExpt(PreInitiatedMolecule mol, ParameterSet parameters, List<CombedMolecule> molList) {
        //initialize iparameters for dynamic simulation
        seqLength = ParameterSet.SEQUENCE_LENGTH;
        binSize = ParameterSet.PU_BIN_SIZE;
        int numberBins = seqLength / binSize;
        double nextRecordingTarget;
        int nextMoleculeTargetIndex;
        
        //get list of positions of initiators
        List<Integer> sitePositionList = mol.get();

        //create list of PotentialInitiationSite2 from positions of potential initiations - sitePositionList is list of 
        //mutable containers for states of replication of each potential initiation site - defines intermediate states 
        //of replication
        potentialInitSitesList = new ArrayList<>();
        for (int i = 0; i < sitePositionList.size(); i++) {
            potentialInitSitesList.add(new PotentialInitiationSite2(sitePositionList.get(i), i));
        }
 
        //initialize endpoints of genome - these are special potential initiation sites becasue they are "active",
        //but they don't have forks - they just define the ends of the molecule
        //replications begins with them pointing at each other
        PotentialInitiationSite2 leftEndPIS;  //left end potential initiation site - special - no forks
        PotentialInitiationSite2 rightEndPIS; //right end potential initiation site - special - no forks
        leftEndPIS = new PotentialInitiationSite2(-1, -1);  
        rightEndPIS = new PotentialInitiationSite2(seqLength, seqLength);  
        leftEndPIS.setRightActivePIS(rightEndPIS);
        leftEndPIS.activateTerminus();
        leftEndPIS.setPositionRightFork(-1);
        rightEndPIS.setLeftActivePIS(leftEndPIS);
        rightEndPIS.activateTerminus();
        rightEndPIS.setPositionLeftFork(seqLength);
        
        //List of recorded states
        statesAtRegularIntervals = new ArrayList<>();
        
        //initialize list of replication time
        timeOfReplicatonInBins = IntStream.range(0, numberBins).mapToObj(i -> 0.0).collect(Collectors.toList());
        
        //get own instance of SynthesisCycle
        cycle = new DnaSynthesisCycle2(potentialInitSitesList, leftEndPIS, rightEndPIS); //list of potential initaition sites plus information to reference end points (left end points to right end)
    
        //Set initial recording target         
        nextRecordingTarget = ParameterSet.TIME_RECORDING_START;  //initial recording targe
        nextMoleculeTargetIndex = 0;
        
        //save starting state of replication at zero time
        statesAtRegularIntervals.add(new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, 0, 0));
        
        
        //Replicate the molecule and save intermittant records of state of replication
        int elongationRate = parameters.getElongationRate();
        while ( ! IsEndOfReplication(cycle)) { 
            int nucleotidesOfForkMovementForCycle = calculateNucleotidesOfForkMovementForCycle(parameters);
            double firingProbabilityPercycle = calculateFiringProbabilityForCycle(parameters);             
            cycle.start(firingProbabilityPercycle, nucleotidesOfForkMovementForCycle);
            elapsedTime += (double) nucleotidesOfForkMovementForCycle / elongationRate;
            fxReplicated = cycle.getFractionReplicated();
            
            //Create a state record based on elapsed time or on fraction replicated depending on whether molList passed
            if (molList  == null) {
                if (elapsedTime >= nextRecordingTarget) {
                    ReplicatingMoleculeForCombingExpt.StateRecord currentStateRecord = new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, elapsedTime,fxReplicated);
                    statesAtRegularIntervals.add(currentStateRecord);
                    while (nextRecordingTarget <= elapsedTime) {nextRecordingTarget += ParameterSet.TIME_RECORDING_INTERVAL;}
                }
            } else {
                if (nextMoleculeTargetIndex < molList.size() && fxReplicated >= molList.get(nextMoleculeTargetIndex).getMoleculeProps().getFxReplicated()) {
                    ReplicatingMoleculeForCombingExpt.StateRecord currentStateRecord = new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, elapsedTime, fxReplicated);
                    statesAtRegularIntervals.add(currentStateRecord);
                    nextMoleculeTargetIndex++;
                    while (nextMoleculeTargetIndex < molList.size() && molList.get(nextMoleculeTargetIndex).getMoleculeProps().getFxReplicated() <= fxReplicated) {
                        currentStateRecord = new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, elapsedTime, fxReplicated);
                        statesAtRegularIntervals.add(currentStateRecord);
                        nextMoleculeTargetIndex++;
                    }                    
                }
            }
            
        }
        //System.out.println(elapsedTime);
        //record state at completion of replication if not already recorded
        while (statesAtRegularIntervals.size() < 161) {statesAtRegularIntervals.add(new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, elapsedTime, 1.0));}
        //if (statesAtRegularIntervals.get(statesAtRegularIntervals.size() - 1).getNucleotidesReplicated() != seqLength) {statesAtRegularIntervals.add(new ReplicatingMoleculeForCombingExpt.StateRecord(cycle, elapsedTime, 1.0,1.0));}
        if (statesAtRegularIntervals.size() != 161) {throw new IndexOutOfBoundsException("size not 161");}
    }
    
    /**
     * Test that replication is complete
     * @return true if complete
     */
    private boolean IsEndOfReplication(DnaSynthesisCycle2 cycle) {
        //test if complete genome is replicated - if not return false
        //if complete genome replicated, set field timeOfCompletion and return true
        //return elapsedTime ;  //can make this more general - use Fx replicated or whatever
        return !(cycle.getNumberNucleotidesReplicated() < seqLength);
    }
    
    /**
     * Calculate number of nucleotides forks that don't terminate will move
     * Now just minPerCycle * elongationRate could be function of replication variables
     * @return default number nucleotides of fork movement per cycle
     */
    private int calculateNucleotidesOfForkMovementForCycle(ParameterSet params) {
        return (int) (params.getMinPerCycle() * params.getElongationRate());
    }
 
    /**
     * Calculate firing probability per potential per cycle
     * If time constant is zero, then calculate product of minPerCycle and maxFiringProbabillityPerMin
     * If time constant is not zero, approach maximum firing probability per cycle exponentially using value of time constant
     * @param parmas the parameters for DNA synthesis
     * @return firing probability per potential per cycle
     */
    private double calculateFiringProbabilityForCycle(ParameterSet params) {  
        //if (params.getTimeConstant() == 0 || elapsedTime > params.getTimeConstant()) {return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin();} //original
        if (elapsedTime * params.getTimeConstant() > params.getMaxFiringProbabilityPerMin()) {return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin();}
        //return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin() * (1 - Math.exp(-(elapsedTime + params.getMinPerCycle()) / params.getTimeConstant()));
        //return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin() * (1 - Math.exp(-elapsedTime / params.getTimeConstant()));
        //return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin() * elapsedTime / params.getTimeConstant(); //original
        return params.getMinPerCycle() * elapsedTime * params.getTimeConstant();
        //return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin() * elapsedTime * elapsedTime / (elapsedTime * elapsedTime + params.getTimeConstant());
    }
    
    /**
     * Returns the list of recorded states
     * @return the list of recorded states
     */
    List<ReplicatingMoleculeForCombingExpt.StateRecord> getRecordedStates () {return statesAtRegularIntervals;}
    
    /**
     * returns elapsed time of replication
     * @return elapsed time
     */
    double getElapsedTime() {return elapsedTime;}
        
    /**
     * returns the length of the replicating sequence
     * @return sequence length in nucleotides
     */
    int getSeqLength() {
        return seqLength;
    }
    
    /**
     * Returns the synthesis cycle engine
     * @return synthesis cycle
     */
    DnaSynthesisCycle2 getCycle() {
        return cycle;
    }
    
    /**
     * Returns the list of potential Initiation sites which are containers for the state of each potential origin of this molecule during replication
     * @return list of potential initiation sites for this molecule
     */
    List<PotentialInitiationSite2> getPotentialInitSitesList() {  //used only for testing
        return potentialInitSitesList;
    }
    
    /**
     * Prints the major variables in the given synthesis cycle - mainly for testing
     * @param cycle the instance of DNASynthesisCycle of this molecule
     */
    void printState(DnaSynthesisCycle2 cycle) {
        System.out.println();
        System.out.println("Number Nucleotides Replicated   " + cycle.getNumberNucleotidesReplicated());
        System.out.println("Number Intiations   " + cycle.getInitiations());
        System.out.println("Number closures   " + cycle.getClosures());
        System.out.println("NumberTerminations   " + cycle.getTerminations());
        System.out.println("Number Passives   " + cycle.getPassives());
        System.out.println("Number Actives   " + cycle.getActives());
        System.out.println("Number Potentials   " + cycle.getPotentials());
        System.out.println("Elapsed Time   " + elapsedTime);
        System.out.format("Fraction Replicated   %.3f%n" ,(double)cycle.getNumberNucleotidesReplicated()/seqLength);
        System.out.println();
        
        if (!cycle.getLeftEnd().getRightActivePIS().IsLeftForkActive()) { //first active left fork is  not active 
            System.out.format("%10d   ", (int) 1);
            }
        potentialInitSitesList.stream().forEach(p -> {
            if (p.IsActive()) {
                if (p.IsLeftForkActive())     {
                    System.out.format("%10d   ",p.getPositionLeftFork());
                }
                if (p.IsRightForkActive()) {
                    System.out.format("%10d   %n",p.getPositionRightFork());
                }
            }
        });
        if (!cycle.getRightEnd().getLeftActivePIS().IsRightForkActive()) {
            System.out.format("%10d   %n",seqLength);
        }
    }
    
    
    
 /**
 * Inner Class for recording intermediate states of replication
 */
    
    public class StateRecord {
        private final int numberActive; //number of active potential initiation sites (pre-RCs)
        private final int numberInitiations; //number of potential initiation sites (pre-RCs) that have initiated
        private final int numberTerminations; //number of potential initiation sites (pre-RCs) that have terminated
        private final int numberPotentials; //number of potential initiation sites (pre-RCs) that are still potential
        private final int  numberPassives; ////number of potential initiation sites (pre-RCs) that have been passively replicated
        private final int numberForks;  //the number of forks
        private final int numberClosures; //the number of closures (convergence of forks)
        private  final double time;  //the time record was recorded
        private final List<Integer> segmentList; //the list of replicated segments in the molecule
        private final List<Integer> replicatedBinList;  //list of values indicating whether bin is more than 50% replicated
        private final int nucleotidesReplicated; //th enumber of nucleotides replicated
        private final int sequenceLength; //the length of the sequence
        private final double fxReplicated; //the fraction replicated
        
        /**
         * Creates a record of the states of the molecule at a given time
         * @param cycle the synthesis engine - an instance of DnaSynthesisCycle2
         * @param totalTime the elapsed time of replication
         * @param fxReplicated the current fraction replicated
         */
        public StateRecord(DnaSynthesisCycle2 cycle, double totalTime, double fxReplicated) {
            sequenceLength = cycle.getRightEnd().getPosition();
            this.time = totalTime;
            numberInitiations = cycle.getInitiations();
            numberClosures = cycle.getClosures();
            numberActive = cycle.getActives();
            numberForks = cycle.getForks();
            numberPassives = cycle.getPassives();
            numberPotentials = cycle.getPotentials();
            numberTerminations = cycle.getTerminations();
            nucleotidesReplicated = cycle.getNumberNucleotidesReplicated();
            this.fxReplicated = fxReplicated;
            segmentList = new ArrayList<>();
            replicatedBinList = new ArrayList<>();
            makeSegmentList();
            
        }        

        /**
         * Makes a list of replicated segments at the elapsed time of this Record 
         * @param siteList  the list of potential Initiations sites at the elapsed time of replication of this record
         * @param cycle  the instance of DNASynthesis Cycle associated with this  replicating molecule
         * @return the list of segments
         */
         private void makeSegmentList() {
            if (cycle.getLeftEnd().getRightActivePIS() == cycle.getRightEnd()) {
                segmentList.add(0);
                replicatedBinList.add(0);
                segmentList.add(sequenceLength - 1);
                replicatedBinList.add((sequenceLength - 1) / binSize);
            } else {
                if ( ! cycle.getLeftEnd().getRightActivePIS().IsLeftForkActive()) {
                     segmentList.add(0);
                     replicatedBinList.add(0);
                 }
                 potentialInitSitesList.stream().filter(site -> site.IsActive()).forEach(activeSite -> {
                     if (activeSite.IsLeftForkActive()) {segmentList.add(activeSite.getPositionLeftFork()); replicatedBinList.add(activeSite.getPositionLeftFork() / binSize);}
                     if (activeSite.IsRightForkActive()) {segmentList.add(activeSite.getPositionRightFork()); replicatedBinList.add(activeSite.getPositionRightFork() / binSize);}
                 });
                 if ( ! cycle.getRightEnd().getLeftActivePIS().IsRightForkActive()) {
                     segmentList.add(sequenceLength - 1);
                     replicatedBinList.add((sequenceLength - 1) / binSize);
                 }
            }
        } 
       
        // return information from recorded state

        /**
         * returns the number initiations
         * @return number initiations
         */   
        public int getNumberInitiations() {return numberInitiations;}

        /**
         * returns the elapsed time
         * @return the elapsed time of replication
         */  
        public double getTime() {return time;}

        /**
         * returns the list of replicated segments
         * @return segment list
         */  
        public List<Integer> getSegmentList () {return segmentList;}

        /**
         * returns the number active potential initiation sites
         * @return number active
         */  
        public int getNumberActive() {return numberActive;}

        /**
         * returns the number terminated potential initiation sites
         * @return number terminated
         */  
        public int getNumberTerminations() {return numberTerminations;}

        /**
         * returns the number potential initiation sites remaining
         * @return number potentials
         */  
        public int getNumberPotentials() {return numberPotentials;}

        /**
         * returns the number potential initiation sites that have been passively replicated
         * @return number passives
         */  
        public int getNumberPassives() {return numberPassives;}

        /**
         * returns the number nucleotides replicated
         * @return number nucleotides replicated
         */  
        public long getNucleotidesReplicated() {return nucleotidesReplicated;}

        /**
         * returns the fraction replicated
         * @return fraction replicated
         */  
        public double getFractionReplicated() {return (double) nucleotidesReplicated/seqLength;}

        /**
         * returns the number forks
         * @return number forks
         */  
        public int getNumberForks() {return numberForks;}

        /**
         * returns the number closed potential initiation sites
         * @return number closed
         */  
        public int getNumberClosures() {return numberClosures;}


         /**
         * Convenience routine to print formated version of the major state record variables
         * @param record 
         */
        void printGenomeState (ReplicatingMolecule2.StateRecord2 record) {
            System.out.format("%30s  %.2f%n", "Time",record.getTime());
            System.out.format("%30s  %d%n", "Number of Initiations",record.getNumberInitiations());
            System.out.format("%30s  %d%n", "Number of Actives",record.getNumberActive());
            System.out.format("%30s  %d%n", "Number of Closures",cycle.getClosures());
            System.out.format("%30s  %d%n", "Number of Terminations",record.getNumberTerminations());
            System.out.format("%30s  %d%n", "Number of Passives",record.getNumberPassives());
            System.out.format("%30s  %d%n", "Number of Potentials",record.getNumberPotentials());
            System.out.format("%30s  %.3f%n", "NucleotidesReplicated", ((double) record.getNucleotidesReplicated()));
            System.out.format("%30s  %.3f%n", "FractionReplicated", ((double) record.getNucleotidesReplicated()/seqLength));
            System.out.format("%30s  %.3f%n", "Number of forks", ((double) record.getNumberForks()));
            System.out.println(); System.out.println();
            //for testing
            //for (int i = 0; i < record.segmentList.size()-1; i=i+2 ) {
                //System.out.println(record.segmentList.get(i) + "   " + record.segmentList.get(i+1));
            //}
            //printState();
        }

    }
}