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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;


/**
 * Rif1 sites described by Hayano et al 2012
 * @author tkelly1
 */
public class Rif1Sites {
    private final List<SeqInterval> rif1IntervalList;
    
    /**
     * Reads coordinates of Rif1 sites from file into a List of SeqInterval
     * @param p path to file containing Rif1 sites
     */   
    public Rif1Sites(Path p) {
        String rawLine;
        BufferedReader inputStream;
        String[] intervalBoundaries;
        SeqInterval interval;
        rif1IntervalList = new ArrayList<>();
        try {
            inputStream = Files.newBufferedReader(p);
            while((rawLine = inputStream.readLine()) != null) {
                intervalBoundaries = rawLine.split(",");
                int start = Integer.parseInt(intervalBoundaries[0]);
                int end = Integer.parseInt(intervalBoundaries[1]);
                rif1IntervalList.add(new SeqInterval(start, end));
            }
        } catch(IOException x) {
                System.err.format("Problem opening Rif1 data file: %s%n", x);                    
        }
    }

    
    /**
     * Returns the list of rif1 sites
     * @return list of SeqInterval describing Rif1 Sites in chromosome
     */
    public List<SeqInterval> getRif1IntervalList() {
        return rif1IntervalList;
    }
}
