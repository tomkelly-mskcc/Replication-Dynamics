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

/**
 * A segment of DNA in a combing experiment - used for lists of various types (N, G, or IOD)
 * segments are given in nucleotides, but resolution of the combing experiment is low
 * consumes smaller memory than CombedMoleculeSegment 
 * @author tkelly
 */

public class CombedMoleculeReplicatedSegment {
    private final double length;  //length of segment
    private final double leftBound;  //left boundary of segment
    private final CombedMolecule molecule; //reference to the molecule containing the segment
    private final int molSegmentNumber; //an index of the segment number in the molecule - really the index of the starts in list of starts,ends....
    
    
    /**
     * Creates instance of RepSegment containing information about the segment
     * @param molSegmentNumber  index of the segment in the molecule (odd numbers)
     * @param length  length of the segment in nucleotides
     * @param leftBound  position of left boundary of segment in molecule
     * @param molecule t5he molecule containing this segment
     */
    public CombedMoleculeReplicatedSegment (int molSegmentNumber, double length, double leftBound, CombedMolecule molecule) {  
        this.length = length;
        this.leftBound = leftBound;
        this.molecule = molecule;
        this. molSegmentNumber = molSegmentNumber;
    }
        
    /**
     * Returns the length of the segment
     * @return segment length
     */
    double length() {
        return length;
    }
    
    /**
     * Returns reference to the combed molecule containing the segment
     * @return the molecule containing the segment
     */
    CombedMolecule getMolecule (){
        return molecule;
    }
    
    /**
     * Returns the number of the segment in the combed molecule
     * @return segment number
     */
    int getMolSegmentNumber() {
        return molSegmentNumber;
    }
    
    /**
     * Returns the position of the left bound of the segment in the combed molecule
     * @return position of left bound in combed molecule
     */
    double getLeftBound () {
        return leftBound;
    }
    
}    
    
    

