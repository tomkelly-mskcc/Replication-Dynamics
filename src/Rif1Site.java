/*
 * The MIT License
 *
 * Copyright 2018 tkelly1.
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
 * A Rif1 site annotated as LE by Hayano et al, 2012
 * @author tkelly
 */
public class Rif1Site {
    
    SeqInterval rif1Interval;
    double rif1RedictedProbability;

    public Rif1Site(SeqInterval rif1Interval, double rif1RedictedProbability) {
        this.rif1Interval = rif1Interval;
        this.rif1RedictedProbability = rif1RedictedProbability;
    }

    public SeqInterval getRif1Interval() {
        return rif1Interval;
    }

    public double getRif1RedictedProbability() {
        return rif1RedictedProbability;
    }

    @Override
    public String toString() {
        return String.format("%d\t%d\t%d\t%f", rif1Interval.start(), rif1Interval.end(), rif1Interval.end() - rif1Interval.start() + 1, rif1RedictedProbability);
    }
    
    
}
    

