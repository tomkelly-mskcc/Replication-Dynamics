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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;

/**
 * A molecule with a fixed number of initiation sites chosen randomly according to probability distribution
 * Probability distribution gives probability of initiation in an initiator site indexed by first nucleotide
 * Index of initiation site is midpoint of corresponding initiator site
 * No overlaps of initiation sites within length of initiator site of each other are allowed
 * List of initiation sites is sorted in increasing sequence index
 * No overlap of initiator sites allowed
 * @author tkelly1
 */
public class PreInitiatedMolecule implements Supplier<List<Integer>>, Serializable {
    
    private final List<Integer> initiatorSitesList;  //unmodifiable sorted list of intiation sites (midpoints of initiator sites)
    private final int sequencelength; //length of sequence

    /**
     * Gets CDF of probability distribution and chooses initiation sites from inverse CDF of random number with uniform probability from 0 to 1
     * @param probDist  an instance of ProbabilityDistribution
     * @param rnd a supplier of a random number
     */
    PreInitiatedMolecule(ProbabilityDistribution2 probDist, Supplier<Double> rnd) {  //Using Dependency Injection for Testing
        int initiatorSiteLength = probDist.getInitiatorSiteLength();
        ParameterSet parameterSet = probDist.getParameterSet();
        List<Double> cdf = probDist.getCdf();
        int halfInitatorLength = initiatorSiteLength / 2;
        sequencelength = ParameterSet.SEQUENCE_LENGTH; 

        //get initiations sites by choosing randomly from probability distribution and preventing overlaps
        int initiatorsLoaded = 0;
        boolean[] disallowedSites = new boolean[cdf.size() + initiatorSiteLength -1];  //to prevent overlaps selected initator sites are set to true - otherwise false
        List<Integer> tempInitiatorSitesList = new ArrayList<>();
        while (initiatorsLoaded < parameterSet.getNumberInitiators()) {  //get number of initiations sites equal to number of initiators
           double r = rnd.get();  //get random number form uniform distribution [0,1]
           int index = Collections.binarySearch(cdf, r);  //inverse cdf to get index of initiator site
           if (index < 0) {index = -index -1;}  //quirk of binary serach engine
           if (!disallowedSites[index] && !disallowedSites[index + initiatorSiteLength - 1]) {  //don't accept initiation site if disallowed - previous on e within iniator site length'
               tempInitiatorSitesList.add(index + halfInitatorLength);
               //System.out.println((index + halfInitatorLength) / 300); //print bin number of preinitation site
               initiatorsLoaded++;
               for (int i = 0; i < initiatorSiteLength; i++) {
                   disallowedSites[index + i] = true;
               }
           }
        }
        Collections.sort(tempInitiatorSitesList);  //sort list of initiation sites
        initiatorSitesList = Collections.unmodifiableList(tempInitiatorSitesList);
        //initiatorSitesList.stream().forEach(e-> System.out.println(e));
    }


    /**
     * returns sorted list of sites of initiation (midpoints of initiator sites
     * @return sorted list of sites of initiation
     */
    @Override
    public List<Integer> get() {  //List of midpoints of initiator sites
        return initiatorSitesList;
    }
    
    /**
     * gets length of sequence of molecule
     * @return the length
     */
    public int getSequencelength() {
        return sequencelength;
    }
}
