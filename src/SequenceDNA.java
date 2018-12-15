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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

/**
 * Immutable class 
 * A DNA sequence
 * provides utility methods for obtaining 
 * properties of sequence
 * @author tkelly
 */
public class SequenceDNA {
    
    private final String sequence;  //The seqeunce of genome/chromosome as a string
    private final int sequenceLength;  //The number of nuycleotides in teh seqeunce
    private final String seqDescription;  //A brief description of the seqeunce
    private Map<Integer,List<Integer>> atFunctionMap = null; //A map of AT contents of initiator sites across the genome
    
    private final int MAX_NUMBER_NUCLEOTIDES = 14000000; //maximum number of nucleotides in seqeunce to be analyzed - not a strong constraint
    private final String NOT_VALID_CHARACTERS = "[^ACGTN]";  //For checking a seqeunce to ensure no non-canoinical letters
    
    /**
     * Creates instance of SequenceDNA from FASTA sequence file
     * @param p path to file
     */
    SequenceDNA(Path p, SeqInterval interval) {
        char[] tempSequence = new char[MAX_NUMBER_NUCLEOTIDES];
        String rawLine;
        BufferedReader inputStream;
        String description = null;
        
        //Get Path for seq file, open a buffered reader, anad read sequence file
        int offset = 0;
        try {
            inputStream = Files.newBufferedReader(p);
            rawLine = inputStream.readLine();
            if (rawLine.startsWith(">")) {
                description = rawLine;
            } else {
                rawLine.getChars(0, rawLine.length(), tempSequence, offset);
                offset += rawLine.length();
            }
            while((rawLine = inputStream.readLine()) != null) {
                rawLine.getChars(0, rawLine.length(), tempSequence, offset);
                offset += rawLine.length();
            }
        } catch(IOException x) {
                System.err.format("Problem opening data file: %s%n", x);                    
        }
        
        //save sequence as a string and check for valid characters
        if (offset < interval.end() - interval.start() + 1)  {throw new IllegalArgumentException("specified seqeunce interval is larger than length of sequence");}
        sequence = (String.copyValueOf(Arrays.copyOfRange(tempSequence, interval.start(), interval.end() + 1))).toUpperCase();
        checkSequence();
        
        //Set getDescription and sequence length
        seqDescription = description == null ? "  interval: " + interval.toString() : description + "  interval: " + interval.toString();
        sequenceLength = sequence.length();
    }
    
    /**
     * Creates instance of SequenceDNA from String containing nucleotides
     * @param description getDescription of sequence 
     * @param nucString sequence as a string
     */
    SequenceDNA(String description, String nucString) {
        this.sequence = nucString;
        checkSequence();
        this.sequenceLength = nucString.length();
        this.seqDescription = description;
    }
    
    
    /**
     * Returns atFunction corresponding to a given initiator site length
     * The function is a list with indices from 0 to sequence length - site length
     * the  index refers to the position of the first nucleotide in the interval of site length nucleotides
     * the algorithm is a bit too fancy - should just call atNumber on intervals of site length at each index 
     * but checked it against the alternative (commented out below)
     * @param initiatorSiteLength
     * @return 
     */
    public List<Integer> atFunction(int initiatorSiteLength) {
        if (atFunctionMap == null) {atFunctionMap = new HashMap<>();}
        if (atFunctionMap.keySet().contains(initiatorSiteLength)) {return atFunctionMap.get(initiatorSiteLength);}
        List<Integer> tempAtFunction = new ArrayList<>();
        //List<Integer> tempAtFunction2 = new ArrayList<>();
        tempAtFunction.add(atNumber(0, initiatorSiteLength - 1)); //number ATs in first initiator site
        //tempAtFunction2.add(atNumber(0, initiatorSiteLength - 1));
        for (int i = 1; i < sequenceLength - initiatorSiteLength + 1; i++) {  //iterates across sequence by one nuc at a time - ++numberAt if new AT at on right end - --numberAT if old AT at left end
            int numberAT = tempAtFunction.get(i - 1);  // number of ATs from i -1 to i + prcSiteLength -1 (total = prcSiteLenth)
            if (isCharAT(i - 1)) {numberAT--;}  // character dropped
            if (isCharAT(i + initiatorSiteLength - 1)) {numberAT++;}  //character added
            tempAtFunction.add(numberAT);
            //tempAtFunction2.add(atNumber(i, i + initiatorSiteLength - 1));
        }
        atFunctionMap.put(initiatorSiteLength, Collections.unmodifiableList(tempAtFunction));
        return atFunctionMap.get(initiatorSiteLength);
    }
    
    /**
     * returns the number of nucleotides in sequence
     * @return number of nucleotides in sequence
     */
    int getLength() {return sequenceLength;}
    
    /**
     * returns the Description of the sequence
     * @return the getDescription of the sequence
     */
    String getDescription() {
        return this.seqDescription;
    }

    /**
     * return the complete sequence as a string
     * @return the sequence
     */
    String getSequence() {
        return sequence;
    }
    
    /**
     * returns true if character at index in sequence is A or T
     * @param index index of character in this sequence
     * @return true if A or T, else false
     */
    boolean isCharAT (int index) {
        return (sequence.charAt(index) == 'A' || sequence.charAt(index) == 'T');
    }
    
    /**
     * returns String containing subsequence corresponding to interval
     * NOte: index of sequence is [0, sequenceLength -1]
     * @param interval  the sequence interval
     * @return the subsequence
     */
    String subSequenceAsString(SeqInterval interval) {
        return sequence.substring(interval.start(), interval.end() + 1);
    }
    
    /**
     * returns new instance of SequenceDNA from subsequence of existing SequenceDNA
     * @param interval the interval of this sequence
     * @return subsequence in interval
     */
    SequenceDNA subSequence(SeqInterval interval) {
        String desc = ">Subsequence of "+seqDescription+" "+interval.toString();
        String seqString = sequence.substring(interval.start(), interval.end() + 1);
        return new SequenceDNA(desc, seqString);
    }
    
    /**
     * returns number of seqBins in this sequence 
     * @param bin instance of seqBin
     * @return number of complete seqBins
     */
    int numberBins (SeqBin bin) {
        return (sequenceLength - bin.offset())/bin.binSize();
    }
   
    /**
     * returns the AT content of an interval of this sequence
     * @param interval the sequence interval defined by instance of SeqInterval
     * @return AT content as fraction
     */
    double atContent(SeqInterval interval) {
        int numberAT = 0;
        for(int i = interval.start(); i < interval.end() + 1; i++) {
            if(sequence.charAt(i) == 'A' || sequence.charAt(i) == 'T') {
                numberAT++;
            }
        }
        return (double)numberAT/(interval.length());
    }
    
    
    /**
     * returns number of ATs in an interval from start to end inclusive of this sequence
     * @param interval  the sequence interval defined start and end indices
     * @return number ATs in interval
     */
    int atNumber(int startIndex, int endIndex) {
        int numberAT = 0;
        for (int i = startIndex; i < endIndex + 1; i++) {
            if(sequence.charAt(i) == 'A' || sequence.charAt(i) == 'T') {
                numberAT++;
            }
        }
        return numberAT;
    }
    
    /**
     * Returns AT content of a sement of sequence given by start and end indices
     * overloads version that takes a seq Interval instead of start and end indices
     * @param startIndex  starting index
     * @param endIndex  ending index
     * @return AT content of segment
     */
    double atContent(int startIndex, int endIndex) {
        int numberAT = 0;
        for(int i = startIndex; i < endIndex + 1; i++) {
            if(sequence.charAt(i) == 'A' || sequence.charAt(i) == 'T') {
                numberAT++;
            }
        }
        return (double)numberAT/(endIndex - startIndex +1);
    }
    
    
    /**
     * Returns the AT content of a bin in a Pu-seq experiment by index of the bin
     * @param binIndex the index - range is 0 to length of seq -1
     * @return AT content of bin
     */
    public double getAtContentPuBin(int binIndex) {
        int start = binIndex * ParameterSet.PU_BIN_SIZE;
        return atContent(start, start + ParameterSet.PU_BIN_SIZE - 1);
    }
    
    /**
     * Returns a list of the AT contents of the bins in the seqeunce
     * @return list of AT contents of bins in sequence
     */
    List<Double> getAtContentInPuBins() {
        return IntStream.range(0, ParameterSet.NUMBER_BINS).mapToDouble(index -> getAtContentPuBin(index)).
                collect(ArrayList::new , ArrayList::add, ArrayList::addAll);
    }
    
    /**
     * returns bin index of a nucleotide number (range of nucleotide numbers is [0,sequenceLength -1]
     * @param nucNumber  the number of a nucleotide in this sequence [0, sequenceLength -1]
     * @param sBin a seqBin giving offset and size of bins
     * @return index of the bin containing the nucleotide at nucNumber
     */
    int nucNumbertoBinIndex(int nucNumber, SeqBin sBin) {
        if (nucNumber < sBin.offset() || nucNumber > (sequenceLength - sBin.offset()) / sBin.binSize()  )  {throw new IllegalArgumentException("Nucleotide number out of range of sequence");}
        return (nucNumber - sBin.offset()) / sBin.binSize();
    }
    
    /**
     * Check sequence for valid characters
     * throws IllegalArgumentException if invalid character found  
     */
    private void checkSequence() {
        Pattern pattern = Pattern.compile(NOT_VALID_CHARACTERS);
        Matcher matcher = pattern.matcher(sequence);
        if (matcher.find()) {
            throw new IllegalArgumentException("Illegal character in sequence");
        }  
    }
     
    /**
     * Returns a string of the contents of a sequence interval (SeqInterval)
     * @param interval an instance of  SeqInterval
     * @return 
     */
    String intervalToString (SeqInterval interval) {
        StringBuilder sb = new StringBuilder();
        sb.append((String.format("%8d\t\t", 1)));
        for (int i = interval.start(); i <= interval.end(); i++) {
            int delta =  i - interval.start();
            if (delta % 10 == 0  && delta  != 0) {sb.append(' ');}
            if (delta % 100 == 0 && delta  != 0 ) {
                //sb.append('\n');
                sb.append(String.format("\n%8d\t", i+1));
                sb.append("\t");
                }
            sb.append(sequence.charAt(i));
        }
        return sb.toString();
    }
    
}
    
    