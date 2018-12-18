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
 * Immutable class defining a sequence bin consisting of a bin size
 * and an offset (starting nucleotide of first bin)
 * @author tkelly1
 */
public final class SeqBin {

    /**
     * Reusable standard SecBin with size 300 and offset zero
     */
    public static final SeqBin BIN_300_0 = new SeqBin(300, 0);
    /**
     * Reusable standard SecBin with size 300 and offset 150
     */
    public static final SeqBin BIN_300_150 = new SeqBin(300, 150);
    /**
     * Reusable standard SecBin with size 600 and offset zero
     */
    public static final SeqBin BIN_600_0 = new SeqBin(600, 0);
    /**
     * Reusable standard SecBin with size 600 and offset 300
     */
    public static final SeqBin BIN_600_300 = new SeqBin(600, 300);
    
    public static final SeqBin BIN_1200_0 = new SeqBin(1200,0);
    
    private final int binSize;
    private final int offset; //bins start at 1 + offset
    
    /**
     * Creates a SeqBin given size and offset
     * @param binSize size of bin in nucleotides
     * @param offset coordinate of first nucleotide in first bin
     */
    SeqBin(int binSize, int offset) {
        this.binSize = binSize;
        this.offset = offset;
    }
    
    /** returns size of bin
     * @return number of nucleotides in bin
     */
    int binSize() {return this.binSize;}
    
    /**
     * returns offset of bin
     * @return offset
     */
    int offset() {return this.offset;}
    
    /**
     * conversion utility returns a SeqInterval corresponding to a 
     * bin with given index
     * @param binIndex the index of the bin
     * @return interval corresponding to bin with given index
     */
    SeqInterval BinToInterval(int binIndex){
        return  new SeqInterval((binIndex) * binSize + offset, (binIndex + 1) * (binSize) + offset -1);
    }
    
    /**
     * conversion utility returns the center index of a SeqBin with given bin index
     * @param binIndex  index of bin with this SeqBin
     * @return center of corresponding interval
     */
    int BinIndexToCenterOfInterval (int binIndex) {
        return binIndex * binSize + offset + (binSize - 1) / 2;
    }
}
