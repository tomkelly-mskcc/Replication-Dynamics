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

/**
 *The properties of a molecule in a combing experiment
 * @author tkelly
 */
public class CombedMoleculeProperties {
    // Identity 
    private final CombedMolecule molecule; //the molecule whose prperties are in this instance
    
    //Basic properties
    private double molLength;
    private double fxReplicated;
    private final List<CombedMoleculeSegment> molSegArray;
    private double startPosition;  //poition of start of molecule in genome - not used to this point set to 0
    private char firstType; //type of the first segment G or N
    private char lastType; //type of the last segment G or N
    private int numberForks; 
    
    //Derivative length arrays
    private List<CombedMoleculeReplicatedSegment> molIODArray;
    private List<CombedMoleculeReplicatedSegment> molGArray;
    private List<CombedMoleculeReplicatedSegment> molNArray;

    
    //IOD Props
    private int numberIODs =0;  //number NGNGNs
    private double averageIODLength;  //defined as intercentroid distance in NGNGN
    
    //G Props
    private int numberGs = 0;  //doesn't include terminal Gs
    private int plusNumberGs;  //includes all Gs including terminal ones
    private double averageGLength;    // Doesn't include terminal Gs
    
    //M props
    private int numberNs = 0; //doesn't include terminal Ns
    private int plusNumberNs;  //includes all Ns including terminal ones
    private double averageNLength;    // Doesn't include terminal Gs
    
    
    
    /**
     * Sets the properties of this molecule
     * @param mol the molecule
     */
    public CombedMoleculeProperties(CombedMolecule mol)  {
        molSegArray = mol.getMolSegArray();  //list of segments from 0 (fist) to length -1 (last)
        molecule = mol;
        SetProps();
    }
    
    
    /**
     * Does the heavy lifting for initializing the properties of the molecule 
     */
    private void SetProps() {
        
        double totalIODLength = 0;
        double totalGLength = 0;
        double totalNLength = 0;
        int numForks = 0;

  
        // declaration of derivative length arrays NOTE these arrays do not contain terminal lengths
        molIODArray =  new ArrayList<>();
        molGArray = new ArrayList<>();
        molNArray = new ArrayList<>();
        
        //Fill derivative length arrays from molecular segment array; terminal lengths not included
        startPosition = molecule.startPosition();
        for (int i = 0; i < molSegArray.size(); i++) {   //this is not as efficient as it could be
            if (molSegArray.get(i).getSegmentType() == 'G' && i != 0 && i+3 < molSegArray.size()) {
                double IODlen = molSegArray.get(i).getSegLength()/2 + molSegArray.get(i+1).getSegLength() + molSegArray.get(i+2).getSegLength()/2;
                totalIODLength += IODlen;
                numberIODs++;
                molIODArray.add(new CombedMoleculeReplicatedSegment(i, IODlen, (molSegArray.get(i).getLeftBound()+ molSegArray.get(i).getRightBound())/2, molecule));
           }
           if (molSegArray.get(i).getSegmentType() == 'G' && i != 0 && i+1 < molSegArray.size())  {
                totalGLength += molSegArray.get(i).getSegLength();
                numberGs++;
               molGArray.add(new CombedMoleculeReplicatedSegment(i, molSegArray.get(i).getSegLength(), molSegArray.get(i).getLeftBound(), molecule));
           }
           if (molSegArray.get(i).getSegmentType() == 'N' && i != 0 && i+1 < molSegArray.size()) {
               totalNLength += molSegArray.get(i).getSegLength();
               numberNs++;
               molNArray.add(new CombedMoleculeReplicatedSegment(i, molSegArray.get(i).getSegLength(), molSegArray.get(i).getLeftBound(), molecule));
           }
        }
        averageIODLength = totalIODLength/numberIODs;
        averageGLength = totalGLength/numberGs;
        averageNLength = totalNLength/numberNs;
        
        //Compute number of forks
        numberForks = 2 * numberGs;
        if (molSegArray.get(0).getSegmentType() == 'G') {numberForks++;}
        if (molSegArray.get(molSegArray.size()-1).getSegmentType() == 'G') {numberForks++;}
        
        //compute total Gs and Ns including terminal ones
        plusNumberGs = numberForks - numberGs;
        plusNumberNs = numberNs + 2 -numberForks +2*numberGs;
        
        //Compute total molecular length and Fx Replicated
        molLength = totalGLength + totalNLength + molSegArray.get(0).getSegLength() + molSegArray.get(molSegArray.size() - 1).getSegLength();
        fxReplicated = (totalGLength + (molSegArray.get(0).segmentType == 'G' ? molSegArray.get(0).getSegLength() : 0) +  
                (molSegArray.get(molSegArray.size() - 1).getSegmentType() == 'G' ? molSegArray.get(molSegArray.size() - 1).getSegLength() : 0))/molLength;
        
        //Define terminal types
        firstType = molSegArray.get(0).getSegmentType();
        lastType = molSegArray.get(molSegArray.size() - 1).getSegmentType();
    }
    
    /**
     * Returns fraction replicated of combed molecule
     * @return fraction replicated
     */
    double getFxReplicated() {
        return fxReplicated;
    }
   /**
    * Returns length of combed molecule
    * @return length of molecule
    */
    double getMolLength() {
        return molLength;
    }
    
    /**
     * Returns the position of start of combed molecule
     * @return start position
     */
    double getStartPosition() {
        return startPosition;
    }
    
    /**
     * Returns the number of IODs in combed molecule
     * @return number of IODs
     */
    int getNumberIODs() {
        return numberIODs;
    }
    
    /**
     * Returns the number of Gs i.e. number of replicated segments
     * @return number of Gs in molecule
     */
    int getnumberGs() {
        return numberGs;
    }
    
    /**
     * Returns the number of forks in the combed molecule
     * @return number of forks in molecule
     */
    int getNumberForks() {
        return numberForks;
    }
    
    /**
     * Returns the number of Ns i.e. the number of segments between replicated segments
     * @return 
     */
    int getnumberNs() {
        return numberNs;
    }
    
    public List<CombedMoleculeReplicatedSegment>  getMolIODArray() {
        return molIODArray;
    }
    
    /**
     * Returns a list of the replicated segments
     * @return replicated segments list
     */
    List<CombedMoleculeReplicatedSegment> getMolGArray() {
        return molGArray;
    }
    
    /**
     * Returns a list of N segments i.e. segments between replicated segments
     * @return N list
     */
    List<CombedMoleculeReplicatedSegment> getMolNArray() {
        return molNArray;
    }
    
    /**
     * Returns the type (N or G) of first segment in combed molecule
     * @return type of first segment in combed molecule
     */
    char getFirstType() {
        return firstType;
    }
    
    /**
     * Returns the type (N or G) of last segment in combed molecule
     * @return type of last segment in combed molecule
     */
    char getLastType() {
        return lastType;
    }
    
    /**
     * Returns a list of the segments in the combed molecule
     * @return 
     */
    public List<CombedMoleculeSegment> getmolSegArray(){
        return molSegArray;
    }
    
    public List<Double> getMolIodLengthList() {
        return molIODArray.stream().map(CombedMoleculeReplicatedSegment::length).collect(Collectors.toList());
    }
}
