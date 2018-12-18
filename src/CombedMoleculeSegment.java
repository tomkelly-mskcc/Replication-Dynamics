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
 * A DNA segment in a combed molecule - used in description of a combed molecule
 * @author tkelly
 */
public class CombedMoleculeSegment {
    int segmentID;  //index of the segment (odd numbers)
    double leftBound;  //left boundary of segment in kb (offset from start)
    double rightBound; //right bound of segment in kb (offset from start)
    double segLength; //segment length in kb
    double segMidpoint; //segment midpoint in kb (offset from start)
    char segmentType; //(N or G or IOD)
    
    public CombedMoleculeSegment (int segmentID, double leftBound, double rightBound, char segmentType) {
        this .segmentID = segmentID;
        this.leftBound = leftBound;
        this.rightBound = rightBound;
        this.segLength = this.rightBound - this.leftBound;
        this.segMidpoint = (this.leftBound + this.rightBound)/2;
        this.segmentType = segmentType;
    }
    
    /**
     * Returns the type of the segment (N or G)
     * @return 
     */
    char getSegmentType() {
        return segmentType;
    }
    
    /**
     * Returns length of segment in nucleotides
     * @return segment length
     */
    double getSegLength() {
        return segLength;
    }
    
    /**
     * Returns index of segment in combed molecule
     * @return index of segment
     */
    int getSegID() {
        return segmentID;
    }
    
    /**
     * Returns left boundary of segment in combed molecule
     * @return left boundary of segment
     */
    double getLeftBound() {
        return leftBound;
    }
    
    /**
     * Returns right boundary of segment in combed molecule
     * @return right boundary of segment
     */
    double getRightBound() {
        return rightBound;
    }

    /**
     * Returns midpoint of segment in combed molecule
     * @return segment midpoint
     */
    double getSegMidpoint() {
        return segMidpoint;
    }
}
