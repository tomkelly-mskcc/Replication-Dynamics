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

/**
 * A molecule in a combing experiment containing replicated segments
 * Experimental data is in csv files - for each molecule, first row is numerical data, second row is type of data point.  
 * Routine expects  first data type is A, which gives starting coordinate.
 * @author tkelly
 */
public class CombedMolecule {
    private final String exptName;  //The name of the combing experiment - from Kaykov and Nurse
    private final int exptMolNumber;  //The number of the molecule
    private final CombedMoleculeProperties molProps;  //The properties of the molecule
    private double startPosition;  //The starting position of the coordinate 0 (usually) or some absolute position from hybridization to known marker
    private final List<CombedMoleculeSegment> molSegArray;  //list of segments from 0 (first) to length-1 (last) - segments have a type 
    
    public CombedMolecule(String exptName, int exptMolNumber, String[] lin1, String[] lin2) {
        this.exptName = exptName;
        this.exptMolNumber = 0;
        this.molSegArray = new ArrayList<>();
        
        //ParseDataLines(lin1, lin2);  // converts strings in file to doubles (numerical data) and char (types)
        //FillSegmentArray();
        if (lin2[0].charAt(0) == 'A') {
            startPosition = Double.parseDouble(lin1[0]) == 1 ? 0 : Double.parseDouble(lin1[0]);
        } else {
            System.err.println("No Molecular start position (A)%n");
        }
        double leftBound = startPosition;
        for (int i=1; i < lin1.length; i++) {
            double length = Double.parseDouble(lin1[i]);
            molSegArray.add(new CombedMoleculeSegment(i-1, leftBound, leftBound + length, lin2[i].charAt(0))); 
            leftBound += length;
        }
        molProps = new CombedMoleculeProperties(this);
    }
    
    /**
     * Returns the name of the combing experiment
     * @return name of experiment
     */
    String getExptName() {
        return exptName;
    }
    
     /**
     * Returns the number of the combed molecule
     * @return molecule number
     */
    int getMolNumber() {
        return exptMolNumber;
    }
    
    /**
     * Returns the starting position of molecule in genome
     * @return position of the first nucleotide in molecule
     */
    double startPosition () {
        return startPosition;
    }
    
    /**
     * Returns an instance of Combed MoleculeProperties giving the properties of the molecule
     * @return the properties of the molecule
     */
    CombedMoleculeProperties getMoleculeProps() {
        return molProps;
    }
   
    /**
     * Returns a list of the segments in the molecule
     * @return an array of segments in a combed molecule
     */
    List<CombedMoleculeSegment> getMolSegArray() {
        return molSegArray;
    }
}
