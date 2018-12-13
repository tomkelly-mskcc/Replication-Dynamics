//Copyright Thomas Kelly, Sloan Kettering Institute, New York, NY 10065

package replicationdynamics;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import static java.util.stream.Collectors.*;

/**
 * main routine for running simulations of combed molecules. 
 * predicts various replication variables for a set of 160 molecules with same fraction replicated
 * as the 160 molecules observed by Kaykov and Nurse
 * Note replication parameters are in instance of ParameterSet and
 * recipes for scanning parameters are in an instance of DefaultParameterRecipies
 * Note: locations of all essential files in Class ParameterSet
 * @author tkelly
 */
public class CombingExperimentSimulation {
    
    public enum Tasks   {   PRINT_COMBED_MOLECULE_PROPERTIES,
                            PRINT_SEQUENCE,
                            PRINT_PARAMETER_SET,
                            PRINT_RMS_DEVIATION_OBS_VS_PREDICTED_FORK_NUMBERS,
                            PRINT_OBSERVED_AND_PREDICTED_FORK_NUMBERS,
                            PRINT_OBSERVED_AND_PREDICTED_CCDFS_ALL,
                            PRINT_OBSERVED_AND_PREDICTED_CCDFS_LESS_THAN_TEN_PERCENT,
                            PRINT_REPLICATION_DYNAMICS,
                            PRINT_ELAPSED_TIME_HISTOGRAM,
                            PRINT_AVERAGE_OBSERVED_NUMBER_FORKS_VS_FRACTION_REPLICATED
                        }
    /**************************
    //SET TASK FOR MAIN ROUTINE
    **************************/
    static Tasks task = Tasks.PRINT_REPLICATION_DYNAMICS;
    
    
    private final static double FX_REPLICATED_BIN_SIZE = .05; //constant for binning molecules by fraction replicated
    
    private static List<CombedMolecule> sortedCombedMolList; //list of observed CombedMolecules sorted according to fraction replicated
    private static List<CombedMolecule> combedMolList;  //unsorted list of observed CombedMolecules
    private static SequenceDNA seq;  //The DNA seqeunce - an instance of SequenceDNA
    private static ReplicatingMoleculePopulationForCombingExpt population;  //the poopulation of replicating molecules with same set of parameters
    private static CombingDataKN combingData; //the observed combing data
    private static ProbabilityDistribution2 probDistribution;  //the probability distribution (uniform for combing expts)

   
    public static void main(String[] args) throws IOException {
                
        //get list of combed molecules from Kaykov and Nurse data - a singleton and sorted list in order of FxReplicated
        combingData = new CombingDataKN(ParameterSet.COMBING_DATA_FOLDER);
        combedMolList = combingData.getMoleculeArray();
        sortedCombedMolList = combedMolList.stream().
                sorted(Comparator.comparingDouble(m-> m.getMoleculeProps().getFxReplicated())).
                filter(e-> e.getMoleculeProps().getFxReplicated() < 1.0).collect(toList());
        
        
        // set parameter set list
        ParameterSetListMaker maker = new ParameterSetListMaker(DefaultParameterRecipes.getDEFAULT_PARAMETER_RECIPES());
        List<ParameterSet> psList = maker.getParameterSetList();
        
        switch  (task) {
            case PRINT_COMBED_MOLECULE_PROPERTIES:
                //Print the properites of Kaykov and Nurse Molecules
                printCombedMoleculeProperties(combedMolList);
                break;
            case PRINT_SEQUENCE:
                Path p1 = Paths.get(ParameterSet.FA_SEQ_FILE);
                seq = new SequenceDNA(p1, new SeqInterval(ParameterSet.SEQUENCE_START, ParameterSet.SEQUENCE_END));
                System.out.println(seq.intervalToString(new SeqInterval(0, seq.getLength() - 1)));
                break;
           }        
        
        //Iterate over all parameter sets
        for (ParameterSet paramSet : psList) {

            //Create uniform probability distribution
            probDistribution = new ProbabilityDistribution2(ParameterSet.SEQUENCE_LENGTH, paramSet);

            //Replicate molecules - create new replicating molecule population with same probability distribution and parameters
            population = new ReplicatingMoleculePopulationForCombingExpt(probDistribution,
                    paramSet, sortedCombedMolList);

            switch  (task) {
                case PRINT_PARAMETER_SET:
                    System.out.println(paramSet.toString());
                    break;
                case PRINT_RMS_DEVIATION_OBS_VS_PREDICTED_FORK_NUMBERS:
                    printRmsDeviationObservedAndPredictedForkNumbers(e-> e.getReferenceMoleculeFxReplicated() <1.0, paramSet);
                    break;
                case PRINT_OBSERVED_AND_PREDICTED_FORK_NUMBERS:
                    System.out.println(paramSet.toString());
                    printObservedVsSimulatedForkNumbers(e->e.getReferenceMolecule().getMoleculeProps().getFxReplicated() < 1.0); 
                    break;
                case PRINT_OBSERVED_AND_PREDICTED_CCDFS_ALL:
                    System.out.println(paramSet.toString());
                    printObservedAndPredictedCcdfs(e->e.getReferenceMolecule().getMoleculeProps().getFxReplicated() > 0
                            && e.getReferenceMolecule().getMoleculeProps().getFxReplicated() <= 1.0);
                    break;
                case PRINT_OBSERVED_AND_PREDICTED_CCDFS_LESS_THAN_TEN_PERCENT:
                    System.out.println(paramSet.toString());
                    printObservedAndPredictedCcdfs(e->e.getReferenceMolecule().getMoleculeProps().getFxReplicated() > 0
                            && e.getReferenceMolecule().getMoleculeProps().getFxReplicated() <= 0.1);
                    break;
                case PRINT_REPLICATION_DYNAMICS:
                    System.out.println(paramSet.toString());
                    population.printAverageStateVariables();
                    break;
                case PRINT_ELAPSED_TIME_HISTOGRAM:
                    printElapsedTimeHistogram();
                    break;
                case PRINT_AVERAGE_OBSERVED_NUMBER_FORKS_VS_FRACTION_REPLICATED:
                    printObservedNumberForksPerChromosomeHistogram();
            }
        }
    }
    
    /**
     * Prints the properties of observed combed molecules in Kaykov and Nurse dataset
     * @param molList the list of combed molecules
     */
    public static void printCombedMoleculeProperties(List<CombedMolecule> molList) {
        System.out.println("COMBED MOLECULE PROPERTIES");
        System.out.println("number\texperiment\tfxReplicated\tlength\tnumberIODs\tnumberForks");
        for (CombedMolecule mol : molList) {
            int molNumber = mol.getMolNumber();
            CombedMoleculeProperties props = mol.getMoleculeProps();
            double fxReplicated = props.getFxReplicated();
            double length = props.getMolLength();
            String expt = mol.getExptName();
            System.out.println(molNumber + "\t" + expt + "\t" + fxReplicated + "\t" + length
                + "\t" + props.getNumberIODs() + "\t" + props.getNumberForks());        
        }
    }
    
    /**
     * Prints a histograms of total replication times (length of S phase)
     */
    public static void printElapsedTimeHistogram() {
        double HISTOGRAM_BIN_SIZE = 1.0;
        List<Double> timeList = population.getElapsedTimeList();
        Map<Integer, Long> timeMap = timeList.stream().collect(groupingBy(t-> (int)(t / HISTOGRAM_BIN_SIZE), counting()));
        System.out.println("DISTRIBUTION OF ELAPSED TIME");
        System.out.println("midpoint\tnumber"); 
        for (int key : timeMap.keySet()) {
            System.out.println((key + .5) * HISTOGRAM_BIN_SIZE + "\t" + timeMap.get(key));
        }
    }
    
    /**
     * Prints a histogram of the average number of forks as function of fraction replicated
     * Bin size for fraction replicated given in FX_REPLICATED_BIN_SIZE
     */
    public static void printObservedNumberForksPerChromosomeHistogram() {
        System.out.println("AVERAGE NUMBER OBSERVED FORKS PER CHROMOSOME AS FUNCTION OF FRACTION REPLICATED");
        System.out.println("midpoint\taveNumber");
        Map<Integer, List<CombedMolecule>> molMap = sortedCombedMolList.stream().collect(groupingBy(mol-> (int) (mol.getMoleculeProps().getFxReplicated() / FX_REPLICATED_BIN_SIZE)));
        for(int key : molMap.keySet()) {
            double averageNumberforks = molMap.get(key).stream().collect(averagingDouble(mol -> mol.getMoleculeProps().getNumberForks() * ParameterSet.SEQUENCE_LENGTH / mol.getMoleculeProps().getMolLength() / 1000));
            System.out.println((key + .5) * FX_REPLICATED_BIN_SIZE + "\t" + averageNumberforks);
        }    
    }
    
    /**
     * Prints observed and predicted complementary CDFs of inter-centroid distances
     * @param predictedFxFilter predicate giving filter for fraction replicated
     */
    public static void printObservedAndPredictedCcdfs(
            Predicate<PopulationAveragesCombingSimulation> predictedFxFilter ) {
        
        //get observed and predicted cCDFs
        ReplicatingMoleculePopulationForCombingExpt.ComplementaryCdfs cdfs = population.getPredictedVsObservedComplementaryCdfOfIods(predictedFxFilter);

        
        System.out.println("OBSERVED vs PREDICTED CCDF OF ICDS");
        System.out.println("bin midpoint\tobserved\tpredicted");
        
        for (int i = 0; i < ParameterSet.NUMBER_INTERVALS_FOR_CCDF; i++) {
            System.out.println(i * ParameterSet.DISTANCE_INTERVAL_FOR_CCDF + "\t" + cdfs.getObservedComplementaryCdf().get(i) + "\t" + cdfs.getPredictedComplementaryCdf().get(i));
        }
    }
    
    
    /**
     * Prints RMS deviation of Observed and Predicted Fork numbers for Kaykov and Nurse dataset
     * @param predictedFxFilter predicate giving filter for fraction replicated
     * @param pset the parameter set
     */
    public static void printRmsDeviationObservedAndPredictedForkNumbers (Predicate<PopulationAveragesCombingSimulation> predictedFxFilter, ParameterSet pset) {
        //get averages list and filter it
        List<PopulationAveragesCombingSimulation> averages = population.getAverageStateVariablesForCombedMolecules();
        List<PopulationAveragesCombingSimulation> filteredAverages = averages.stream().filter(predictedFxFilter).collect(toList());
        
        //calculate RMS deviation of predicted vs observed number forks
        double totalSqDeviation = 0.0;
        for (PopulationAveragesCombingSimulation ave : filteredAverages) {
            totalSqDeviation += (ave.getReferenceMoleculeForksPerMb()- ave.getForksPerMb()) * (ave.getReferenceMoleculeForksPerMb()- ave.getForksPerMb());
        }
        System.out.println("RMS deviation"+ "\t" + "RateIncrease" + "\t" + "MaxRate" );
        System.out.println(Math.sqrt(totalSqDeviation / filteredAverages.size()) + "\t" + pset.getTimeConstant() + "\t" + pset.getMaxFiringProbabilityPerMin());
    }
     
    /**
     * Prints observed and predicted fork numbers for Kaykov and Nurse dataset
     * @param fxFilter predicate giving filter on fraction replicated
     */
    public static void printObservedVsSimulatedForkNumbers (Predicate<PopulationAveragesCombingSimulation> fxFilter) {
            
        System.out.println("SIMULATION VS OBSERVED");
        System.out.println("sim Fx\tobs fx\tsim forks\tobs forks\tobs len\tsim fx\tsim forks per Mb\tobs fx\tobs forks per Mb\tsim forks er mMb\tobs forks per Mb");
        population.getAverageStateVariablesForCombedMolecules().stream().filter(fxFilter::test).forEach( e->System.out.println(e.getFxReplicated() + "\t" + e.getReferenceMoleculeFxReplicated()+
                    "\t" + e.getNumberForks() + "\t" + e.getReferenceMoleculeNumberforks() + "\t" + e.getReferenceMoleculeSize() +
                    "\t"  + e.getFxReplicated() + "\t" + e.getForksPerMb() + "\t" + e.getReferenceMoleculeFxReplicated() + "\t" + 
                    + e.getReferenceMoleculeForksPerMb() + "\t" + e.getForksPerMb() + "\t" + e.getReferenceMoleculeForksPerMb()));
    }
    
    
}

        