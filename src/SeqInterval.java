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
 * An interval in a sequence of nucleotides
 * Immutable Class
 * @author tkelly
 */
public class SeqInterval {
    private final int start;
    private final int end;
    
    /**
     * Initializes a closed interval [min, max].
     * @param  start the smaller endpoint
     * @param  end the larger endpoint
     * 
     */
    public SeqInterval(int start, int end) {
        if (start <= end) {
            this.start = start;
            this.end = end;
        }
        else throw new IllegalArgumentException("Illegal interval");
    }
    
     /**
     * Returns the first endpoint of this interval.
     * @return the first endpoint of this interval
     */
    int start() { 
        return start;
    }

    /**
     * Returns the last endpoint of this interval.
     * @return the max endpoint of this interval
     */
    int end() { 
        return end;
    }
    
    /**
     * returns center of interval rounded down
     * @return interval center
     */
    int center() {
        return (this.end - this.start)/2 + this.start;
    }

    /**
     * Returns true if this interval intersects the specified interval.
     *
     * @param  that the other interval
     * @return {@code true} if this interval intersects the argument interval;
     *         {@code false} otherwise
     */
    boolean intersects(SeqInterval that) {
        if (this.end < that.start) return false;
        return that.end >= this.start;
    }

    /**
     * Returns true if this interval contains the specified value.
     * @param x the value
     * @return {@code true} if this interval contains the value {@code x};
     *         {@code false} otherwise
     */
    boolean contains(int x) {
        return (start <= x) && (x <= end);
    }

   
    /** Returns the size of this interval i.e. the number of locations in the interval inclusive of end.
   * @return the length of this interval inclusive of end
    */
    int length() {
        return end - start + 1;
    }
    
    
    /**
     * Returns a string representation of this interval.
     * @return a string representation of this interval in the form min tab max
     */
    @Override
    public String toString() {
        return "[" + start + ":" + end + "]" ;
    }

    /**
     * Compares this interval to the specified object.
     * @param  other the other interval
     * @return {@code true} if this interval equals the other interval;
     *         {@code false} otherwise
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) return true;
        if (other == null) return false;
        if (other.getClass() != this.getClass()) return false;
        SeqInterval that = (SeqInterval) other;
        return this.start == that.start && this.end == that.end;
    }

    /**
     * Returns an integer hash code for this interval.
     * @return an integer hash code for this interval
     */
    @Override
    public int hashCode() {
        int hash1 = ((Integer) start).hashCode();
        int hash2 = ((Integer) end).hashCode();
        return 31*hash1 + hash2;
    }
}

    

