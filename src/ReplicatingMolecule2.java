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
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A replicating DNA molecule 
 * this a mutable class
 * @author tkelly
 */
public class ReplicatingMolecule2 {
    
    //Constants of replication from parameters
    private final ParameterSet paramSet;
    private final List<PotentialInitiationSite2> potentialInitSitesList;  //the list of potential initiation sites encoding states of each site during DNA replication
    private final int seqLength;  //the length of the DNA sequeccne of the replicating molecule
    
    // Variables related to replication progress  
    // just need enough to calculate setup for each synthetic cycle
    private double elapsedTime = 0;  //elapsed replication time - increments are calculated by synthesis engine
    private final DnaSynthesisCycle2 cycle;  //the synthesis cycle for this replicating molecule
    
    //Variables needed to record states and put results in bins as needed
    private final List<StateRecord2> statesAtRegularIntervals;  // list of states at regular intervals - now approximately 0.01 fraction replicated
    private final List<Double> moleculeRtForkFreqInBins;  //List rightward fork frequencies in bins
    private List<Double> timeOfReplicatonInBins; //List of times of replication of bins
    private final List<Integer> moleculeTerminationsInBins; //List of number of terminations is bins
    private final List<Integer> moleculeInitiationsInBins; //List of number of intiations in bins
    private final  int binSize; //size of Pu bin (not needed it's in the parameters)
        
    /**
    * Initializes variables for replication of molecule, then replicates it
    * * Takes a preInitiated molecule and a set of replication parameters
    * Creates a list of PotentialInitiationSite each of which holds replication state of site and linkage to tother sites
    * Initializes PortentialIniationsites at the ends of the linear DNA molecule (indices 0 and seqLength + 1) - special because no active forks 
    * Creates an instance of SynthesisCycle which updates the list of PotentialIntiationSites (states) for each short synthetic period
    * Creates a list of state records giving selected intermediate states and final state in replication
    * @param mol the pre-initiated molecule - a sorted list of initiation sites
    * @param rif1List list of Rif1 sites for suppression of firing
     * @param parameters the parameters for replication
     */
    public ReplicatingMolecule2(PreInitiatedMolecule mol, ParameterSet parameters, List<SeqInterval> rif1List) {
        
        //initialize iparameters for dynamic simulation
        paramSet = parameters;
        seqLength = ParameterSet.SEQUENCE_LENGTH;
        binSize = ParameterSet.PU_BIN_SIZE;
        int numberBins = seqLength / binSize;
        double nextRecordingTarget;
        //eliminate positions of potential initiation from site positions list based on Rif1 binding
        List<Integer> tempSitePositionList = mol.get();
        List<Integer> sitePositionList = getRif1AdjustedSiteList(tempSitePositionList, rif1List);
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
        timeOfReplicatonInBins = null;
        //get own instance of SynthesisCycle
        cycle = new DnaSynthesisCycle2(potentialInitSitesList, leftEndPIS, rightEndPIS); //list of potential initaition sites plus information to reference end points (left end points to right end)
    
        //Set initial recording target         
        nextRecordingTarget = ParameterSet.TIME_RECORDING_START;  //initial recording targe
        
        //save starting state of replication at zero time
        statesAtRegularIntervals.add(new StateRecord2(cycle, elapsedTime));
                
        //Replicate the molecule and save intermittant records of state of replication
        int elongationRate = parameters.getElongationRate();
        while ( ! IsEndOfReplication(cycle)) { 
            int nucleotidesOfForkMovementForCycle = calculateNucleotidesOfForkMovementForCycle(parameters);
            double firingProbabilityPercycle = calculateFiringProbabilityForCycle(parameters);             
            cycle.start(firingProbabilityPercycle, nucleotidesOfForkMovementForCycle);
            elapsedTime += (double) nucleotidesOfForkMovementForCycle / elongationRate; //this is just the same as the length of the interval min per cycle
            
            //Create a state record based on elapsed time
            if (elapsedTime >= nextRecordingTarget) {
                StateRecord2 currentStateRecord = new StateRecord2(cycle, elapsedTime);
                statesAtRegularIntervals.add(currentStateRecord);
                

                
                //updateTimeOfReplication(currentStateRecord);
                while (nextRecordingTarget <= elapsedTime) {nextRecordingTarget += ParameterSet.TIME_RECORDING_INTERVAL;}
            }
             
        }
        
        //record state at completion of replication if not already recorded
        if (statesAtRegularIntervals.get(statesAtRegularIntervals.size() - 1).getNucleotidesReplicated() != seqLength) {
            statesAtRegularIntervals.add(new StateRecord2(cycle, elapsedTime));
            //updateTimeOfReplication(statesAtRegularIntervals.get(statesAtRegularIntervals.size() - 1));
        }        
        
        moleculeRtForkFreqInBins = getMoleculeRightForkFrequencyDistribution();
        moleculeTerminationsInBins = getMoleculeTerminationDistribution();
        moleculeInitiationsInBins = getMoleculeInitiationDistribution();
    }   
    
    /**
     * creates new sitePosition list lacking any sites in rif1 intervals
     * @param sitePositionList the original list of pre-initiated sites
     * @param rif1Intervals intervals containing Rif1 that significantly alter firing rate
     * @return 
     */
    private List<Integer> getRif1AdjustedSiteList(List<Integer> sitePositionList, List<SeqInterval> rif1Intervals) {
        boolean[] isRif1SiteArray = new boolean[ParameterSet.SEQUENCE_LENGTH];
        for (SeqInterval interval : rif1Intervals) {
            for (int i = interval.start(); i <= interval.end(); i++) {
                isRif1SiteArray[i] = true;
            }
        }
        return sitePositionList.stream().filter(pos-> ! isRif1SiteArray[pos]).collect(Collectors.toList());
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
     * If time constant is not zero, approach maximum firing probability per cycle linearly from zero min to time to reach maximum given by time constant
     * @param parmas the parameters for DNA synthesis
     * @return firing probability per potential per cycle
     */
    private double calculateFiringProbabilityForCycle(ParameterSet params) {  
        if (elapsedTime * params.getTimeConstant() > params.getMaxFiringProbabilityPerMin()) {return params.getMinPerCycle() * params.getMaxFiringProbabilityPerMin();}    
        return params.getMinPerCycle() * elapsedTime * params.getTimeConstant();    
    }
        
    /**
     * Updates a list of the replication time of each bin for this molecule
     * Updating occurs lazily if needed by call to getTimeOfReplicationInBins
     * The time of replication is the time at which a fork first enters the bin
     * Algorithm sets any element of replicated bin list within a replicated sement that is equal to 0.0 (not previously updated) to the current time
     * @param record a state record - called with this set to current state recored
     */
    private void updateTimeOfReplication(StateRecord2 record) {
        //get the list of bins that have replicated as of this record
        List<Integer> replicatedBinList = record.getReplicatedBinList();
        
        if (record.nucleotidesReplicated == 0) {return;}
        //iterate over replicatedBinList
        for (int i = 0; i < replicatedBinList.size() - 1; i+=2) {
            int leftIndex = replicatedBinList.get(i); //Bin icontaining left end of segment
            int rightIndex = replicatedBinList.get(i + 1); //Bin containing right end of segment
            //iterate over segment setting any zero elements to elapsed time
            for (int j = leftIndex; j <= rightIndex; j++) {
                if (timeOfReplicatonInBins.get(j) == 0) {timeOfReplicatonInBins.set(j, record.getTime());}
            }
        }
    } 
    
    /**
     * Returns the list of recorded states
     * @return the list of recorded states
     */
    List<StateRecord2> getRecordedStates () {return statesAtRegularIntervals;}
    
    /**
     * returns elapsed time of replication
     * @return elapsed time
     */
    double getElapsedTime() {return elapsedTime;}

    /**
     * Returns the synthesis cycle engine
     * @return synthesis cycle
     */
    DnaSynthesisCycle2 getCycle() {
        return cycle;
    }
    
    /**
     * Returns a list containing the time of replication of each bin for this molecule
     * @return 
     */
    List<Double> getTimeOfMolecluleReplicatonInBins() {
        if (timeOfReplicatonInBins == null) {
            //initialize elements of time array to 0.0
            timeOfReplicatonInBins = IntStream.range(0, ParameterSet.NUMBER_BINS).mapToObj(i -> 0.0).collect(Collectors.toList());
            //iterate over all state records after the one prior to replication
            for (int i = 1; i < statesAtRegularIntervals.size(); i++) {
                updateTimeOfReplication(statesAtRegularIntervals.get(i));
            }
        }        
        return timeOfReplicatonInBins;
    }

    /**
     * Given a list of double values in  bins the method updates the list adding the fraction of nucleotides
     * in each bin covered by a DNA sement from segmentStart to segmentEnd
     * nucleotides in sequence are indexed [0, sequenceLength -1]
     * first bin is interval in the sequence [0, binSize -1]
     * @param segmentStart index of starting nucleotide in segment
     * @param segmentEnd  index of ending nucleotide in segment
     * @param binSize size of bin in nucleotides
     * @param binList the list of doubles, one for each bin
     */
    static void nucleotideSegmentToBins (int segmentStart, int segmentEnd, int binSize, List<Double> binList) {
        int binStart = segmentStart / binSize;
        int binEnd = segmentEnd / binSize;
        if (binStart == binEnd) {binList.set(binStart, binList.get(binStart) + (double)(segmentEnd - segmentStart + 1) / binSize);}
        else {
            binList.set(binStart, binList.get(binStart) + (double)((binStart + 1) * binSize - segmentStart) / binSize);
            for (int i = binStart + 1; i < binEnd; i++) {
                binList.set(i, 1.0);
            }
            binList.set(binEnd, binList.get(binEnd) + (double)(segmentEnd - binEnd * binSize + 1) / binSize);
        }
    }
    
    /**
     * Returns a list of doubles giving the frequency of rightward forks in bins across the genome
     * The values are 0 or 1 except at the boundaries where forks terminate within a bin for the most part
     * this is the predicted right fork frequencies to be compared with Carr data by At model comparator
     * This test and shown to work fine - see below
     * @return list of right fork frequencies in bins
     */
    private List<Double> getMoleculeRightForkFrequencyDistribution () {
        List<Double> rtForkBinList = IntStream.range(0, seqLength / binSize).mapToObj(i -> 0.0).collect(Collectors.toList());
        potentialInitSitesList.stream().filter(f-> f.IsTerminated()).forEach(potInitSite -> {
            nucleotideSegmentToBins(potInitSite.getPosition(), potInitSite.getPositionRightFork(), binSize, rtForkBinList);
        });
        return rtForkBinList;
    }
  
    /**
     * Returns a list of number terminations in each bin of molecule
     * @return termination list for molecule
     */
    private List<Integer> getMoleculeTerminationDistribution () {
        List<Integer> terminationBinList = IntStream.range(0, ParameterSet.NUMBER_BINS).mapToObj(e->0).collect(Collectors.toList());
        potentialInitSitesList.stream().filter(f-> f.IsTerminated()).map(init-> 
                init.getPositionRightFork() / ParameterSet.PU_BIN_SIZE).filter(b-> b != ParameterSet.NUMBER_BINS - 1).forEach(bin-> 
                        terminationBinList.set(bin, terminationBinList.get(bin) + 1));
        return terminationBinList;
    }
    
    /**
     * Returns a list of number of initiations in each bin of molecule
     * @return initiation list for molecule
     */
    private List<Integer> getMoleculeInitiationDistribution () {
        List<Integer> initiationBinList = IntStream.range(0, ParameterSet.NUMBER_BINS).mapToObj(e->0).collect(Collectors.toList());
        potentialInitSitesList.stream().filter(f-> f.IsTerminated()).map(init-> 
                init.getPosition() / ParameterSet.PU_BIN_SIZE).forEach(bin-> 
                        initiationBinList.set(bin, initiationBinList.get(bin) + 1));
        return initiationBinList;
    }
    
    /**
     * Returns the list of potential Initiation sites which are containers for the state of each potential origin of this molecule during replication
     * @return list of potential initiation sites for this molecule
     */
    List<PotentialInitiationSite2> getPotentialInitSitesList() {  //used only for testing
        return potentialInitSitesList;
    }
    
   /**
     * Returns the right fork frequency in each PU bin for this molecule
     * @return list containing right fork frequencies in each bin
     */
    List<Double> getMoleculeRtForkFreqInBins() {
        return moleculeRtForkFreqInBins;
    }

    /**
     * Returns a list of number of termination sites in each bin of molecule
     * @return termination list for molecule
     */
    List<Integer> getMoleculeTerminationsInBins() {
        return moleculeTerminationsInBins;
    }

    /**
     * Returns a list of number of initiations in each bin of molecule
     * @return list of initiations
     */
    List<Integer> getMoleculeInitiationsInBins() {
        return moleculeInitiationsInBins;
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
    
    public class StateRecord2 {
        private int numberActive = 0; //number of active potential initiation sites (preRCs)
        private int numberInitiations = 0; //number of initiations that have occured
        private int numberTerminations = 0; //number of terminations that have occured
        private int numberPotentials = 0; //number of potential initiations sites that are still potential
        private int  numberPassives = 0; //number of potential intiations sites that have been passively replicated
        private int numberForks = 0; //number forks
        private int numberClosures = 0; //number closures
        private  final double time; //elapsed time of replication
        private final List<Integer> segmentList; //the list of replicated segemnts
        private final List<Integer> replicatedBinList;  //list of values indicating whether bin is more than 50% replicated
        private int nucleotidesReplicated = 0; //number nucleotides replicated
        private final int sequenceLength; //length of sequence
        
        public StateRecord2(DnaSynthesisCycle2 cycle, double totalTime) {
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
            segmentList = new ArrayList<>();
            replicatedBinList = new ArrayList<>();
            makeSegmentList();
        }
        

        /**
         * Makes a list of replicated segments at the elapsed time of this Record 
         * Units for endpoints of segments are nucleotides and bins
         * @param siteList  the list of potential Initiations sites at the elapsed time of replication of this record
         * @param cycle  the instance of DNASynthesis Cycle associated with this  replicating molecule
         * @return 
         */
       private void makeSegmentList() {
           //if no replicated segments just return - nothing to make
            if (cycle.getNumberNucleotidesReplicated() == 0) {
                return;
            }
            
            //if this is the last record - completely replicated the only segment is the whole genome
            if (cycle.getNumberNucleotidesReplicated() == ParameterSet.SEQUENCE_LENGTH) {
                segmentList.add(0);
                replicatedBinList.add(0);
                segmentList.add(ParameterSet.SEQUENCE_LENGTH - 1);
                replicatedBinList.add(ParameterSet.NUMBER_BINS - 1);
                return;
            } 
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

       
       
    // return information from recorded state - maybe just return the state!!
    
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
     * Returns the list of replicated bins
     * @return list of replicated bins
     */
    public List<Integer> getReplicatedBinList() {return replicatedBinList;}
    }
    
     /**
     * Convenience routine to print formated version of the major state record variables
     * @param record 
     */
    final void printGenomeState (ReplicatingMolecule2.StateRecord2 record) {
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
    
    /**
     * Returns the value of a double valued state variable at a given time (to first approximation)
     * If the time is greater than the replication time for this molecule the method
     * returns the value at the end of replication
     * State variable must be type double
     * @param time time in minutes accuracy is the length of the recording interval
     * @param getStateVariable a Function of the state record returning a state value
     * @return the state value at given time
     */
    double getDoubleStateVariableAtTime (double time, Function<ReplicatingMolecule2.StateRecord2, Double> getStateVariable) {
        //calculate index of state record
        int recordIndex = (int) ((time + paramSet.getMinPerCycle()) / ParameterSet.TIME_RECORDING_INTERVAL);
        //if recordIndex exceeds bound of record states list - set index to highest value
        if (recordIndex > statesAtRegularIntervals.size() -1) {
            recordIndex = statesAtRegularIntervals.size() -1;
        }
        ReplicatingMolecule2.StateRecord2 state = statesAtRegularIntervals.get(recordIndex);
        return getStateVariable.apply(state);
    }
    
    
    /**
     * Returns the value of an integer valued state variable at a given time (to first approximation)
     * If the time is greater than the replication time for this molecule the method
     * returns the value at the end of replication
     * State variable must be type integer
     * @param time time in minutes accuracy is the length of the recording interval
     * @param getStateVariable a Function of the state record returning a state value
     * @return the state value at given time
     */
    double getIntegerStateVariableAtTime (double time, Function<ReplicatingMolecule2.StateRecord2, Integer> getStateVariable) {
        //calculate index of state record
        int recordIndex = (int) ((time + paramSet.getMinPerCycle()) / ParameterSet.TIME_RECORDING_INTERVAL);
        //if recordIndex exceeds bound of record states list - set index to highest value
        if (recordIndex > statesAtRegularIntervals.size() -1) {
            recordIndex = statesAtRegularIntervals.size() -1;
        }
        ReplicatingMolecule2.StateRecord2 state = statesAtRegularIntervals.get(recordIndex);
        return getStateVariable.apply(state);
    }
    
    
    /**
    * Returns the time of replication of a particular bin
    * @param binIndex the index of the bin
    * @return 
    */
    double getTimeOfReplicationOfBin (int binIndex) {
        return getTimeOfMolecluleReplicatonInBins().get(binIndex);
    }
}


  

    
    