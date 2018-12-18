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
 * A table of exponential functions exp(k * AtNumber / siteLength of AT number)
 * i.e exp(k * fractionAT) in putative initiator binding site
 * avoids calculating same exponential multiple times (there are only initiatorSitelLength of them
 * @author tkelly
 */
public class AtExponentialTable {
    private final  List<Double> atExponential;

    /**
     * creates an exponential table
     * @param siteLength the length of the initiator binding site in model
     * @param kValue coefficient of fraction AT in exponent
     */
    public AtExponentialTable(int siteLength, double kValue) {
        if (Double.isInfinite(Math.exp(kValue))) {throw new ArithmeticException("kvalue is too big - result infinite");}
        atExponential = new ArrayList<>();
         for (int i = 0; i < siteLength + 1; i++) {
            double atFraction = (double) i / siteLength;
            atExponential.add(Math.exp(kValue * atFraction));
        }
    }

    /**
     * returns value in exponential table by index
     * @param index number of ATs in InitiatorSiteLength
     * @return 
     */
    double getAtExponential(int index) {
        return atExponential.get(index);
    }
}
    

