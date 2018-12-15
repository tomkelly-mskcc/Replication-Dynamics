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
import java.util.Iterator;
import java.util.List;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
* The annotation file for an underlying DNA sequence 
* @author tkelly
*/
public class SequenceDNAAnnotationGenbank {
    
    private final List<String> sequenceAnnotation;  //A string containing the contents of sequence annotation file
    private static final String INTERVAL_MARKER = "(\\d+)..(\\d+)";  //For parsing annotation file
    private final static String TYPE_MARKER = "/gene=\"(.*)\"";  //For parsing annotation file
    
    /**
     * Constructor reads a Genbank file into a List of String
     * makes a list of Transcript
     * @param p path to the Genbank file
     */
    public SequenceDNAAnnotationGenbank(Path p) {

        sequenceAnnotation = new ArrayList<>();
        String rawLine;
        BufferedReader inputStream;
        

        //Get  buffered reader and read lines of GENBANK file into cdsFileArray
        try {
            inputStream = Files.newBufferedReader(p);
            while((rawLine = inputStream.readLine()) != null) {
                sequenceAnnotation.add(rawLine);
            }
            } catch(IOException x) {
                System.err.format("Problem opening data file: %s%n", x);                    
            }
    }   
    
    /**
    * parses the lines of a Genbank file and makes a list of the Transcripts
    * this is the guts of the class
    * @param isTypeIncludedInList A predicate returning true if an RNA type should be included in the list
    * @return a list of sequence intervals corresponding to  transcripts of the given type
    */
    public List<Transcript> getTranscriptsByType(Predicate<Transcript.TRANSCRIPT_TYPE> isTypeIncludedInList) {
        //Initialize local variables
        List<Transcript> transcriptList = new ArrayList();
        SeqInterval rnaInterval;  
        Transcript.TRANSCRIPT_DIRECTION direction;
        Pattern intervalPattern = Pattern.compile(INTERVAL_MARKER);
        Pattern typePattern = Pattern.compile(TYPE_MARKER);
        int elementStart = 0;
        int elementEnd = 0;
        String description = "default";
        Transcript.TRANSCRIPT_TYPE type;
        SeqInterval fivePimeUtrInterval = null;
        SeqInterval cdsInterval = null;
        SeqInterval threePrimeUtrInterval = null;
        int numberIntrons;
        
        
        //Iterate over lines in annotation file
        for (Iterator<String> iter = sequenceAnnotation.iterator(); iter.hasNext();) {
            String line = iter.next();
            
            //find next instance of "gene" in first 16 characters
            if (line.length() > 15 && line.substring(0, 15).contains("gene")) {
                String geneLine = line;
                line = iter.next();
                
                //skip to RNA line  and concatenate with continuation lines
                while ( ! (line.length() > 15 && line.substring(0,15).contains("RNA"))) {line = iter.next();}
                String rnaLine = line;
                while (rnaLine.endsWith(",")) {
                    rnaLine = rnaLine.concat(iter.next());
                }
                line = iter.next();
                
                //skip to "/gene" line
                while ( ! (line.substring(20,28).contains("/gene"))) {iter.next();}
                String descriptionLine = line;
                line = iter.next();
                
                //if mRNA skip to CDS line and concatentate with continuation lines
                if (rnaLine.contains("mRNA")) {
                    while ( ! (line.length() > 15 && line.substring(0,15).contains("CDS"))) {line = iter.next();}
                    String cdsLine = line;
                    while (cdsLine.endsWith(",")) {
                        cdsLine = cdsLine.concat(iter.next());
                    }
                    cdsInterval = getJoinedExons(cdsLine);
                }
                
                //set type and description
                Matcher matcher2 = typePattern.matcher(descriptionLine);
                while (matcher2.find()) {
                    description = matcher2.group(1);
                }
                if (rnaLine.contains("mRNA")) {
                    type = Transcript.TRANSCRIPT_TYPE.MRNA;
                } else {
                    type = getNonMessengerTranscriptType(description);
                }
                
                //set direction
                if (geneLine.contains("complement")) {direction = Transcript.TRANSCRIPT_DIRECTION.TO_LEFT;}
                else {direction = Transcript.TRANSCRIPT_DIRECTION.TO_RIGHT;}
                
                //get transcript interval and number introns
                numberIntrons = 0;
                Matcher matcher = intervalPattern.matcher(geneLine);
                matcher.find();
                rnaInterval = new SeqInterval(Integer.parseInt(matcher.group(1))  - 1, Integer.parseInt(matcher.group(2)) - 1);
                Matcher matcher3 = intervalPattern.matcher(rnaLine);            
                while (matcher3.find()) {
                   numberIntrons++; 
                }
                
                //if mRNA get 5' cds and 3' subintervals
                if (rnaLine.contains("mRNA")) {
                    fivePimeUtrInterval = getFivePrimeUtr(direction, rnaInterval, cdsInterval);
                    threePrimeUtrInterval = getThreePrimeUtr(direction, rnaInterval, cdsInterval);
                }
                
                //add transcript to list
                if (isTypeIncludedInList.test(type)) {
                    transcriptList.add(new Transcript(rnaInterval, type, direction, description, numberIntrons - 1,
                        fivePimeUtrInterval, cdsInterval, threePrimeUtrInterval));
                }
            }
        }    
        return transcriptList;
    }
    
    
    /**
     *Returns the type of an RNA species (not MRNA) from text in the description 
     * @param description test from line containing "/gene" - annotates type of RNA except MRNA
     * @return 
     */
    private Transcript.TRANSCRIPT_TYPE getNonMessengerTranscriptType(String description) {
        if (description.contains("RRNA")) {return Transcript.TRANSCRIPT_TYPE.RRNA;}
        if (description.contains("SNRNA")){return Transcript.TRANSCRIPT_TYPE.SNRNA;}
        if (description.contains("TRNA")) {return Transcript.TRANSCRIPT_TYPE.TRNA;}
        if (description.contains("SNORNA")) {return Transcript.TRANSCRIPT_TYPE.SNORNA;}
        if (description.contains("NCRNA")) {return Transcript.TRANSCRIPT_TYPE.NCRNA;}
        else {return Transcript.TRANSCRIPT_TYPE.PSEUDOGENE;}
    }
    
    

    /**
     * Used to get block for CDS from set of subintervals - block for genes is already in a single interval
     * Block for genes gives total interval for a transcript in Genbank annotation
     * So this method is just to add CDS interval to Transcript, but this value of Transcript is not used any more
     * At this point the only call to this method is from the gettranscriptsByType Parser
     * @param line line of file (possibly concatenated from multiple lines) containing intervals 
     * @return a SeqInterval of the block from beginning to end of subintervals
     */
    private SeqInterval getJoinedExons(String line) {
        Pattern intervalPattern = Pattern.compile(INTERVAL_MARKER);
        Matcher matcher = intervalPattern.matcher(line);
        int intervalStart = 0;
        int intervalEnd = 0;
        int first;
        int second;
        while (matcher.find()) {
            first = Integer.parseInt(matcher.group(1)) - 1;
            second = Integer.parseInt(matcher.group(2)) - 1;
            if (intervalStart == 0) {
                intervalStart = first;
                intervalEnd = second;
            } else {
                if (first < intervalStart) {intervalStart = first;}
                if (second > intervalEnd) {intervalEnd = second;}
            }    
       }
        return new SeqInterval(intervalStart, intervalEnd);
    }
    
    /**
     * Returns a SeqInterval giving the interval of the five prime UTR of a mRNA
     * used  to add f prime UTR to Transcript, but this value in Transcript not used any more
     * @param direction  direction of transcription left or right
     * @param rnaInterval the interval for the mRNA transcript
     * @param cdsInterval  the interval for the CDS of a mRNA
     * @return 
     */
    private SeqInterval getFivePrimeUtr (Transcript.TRANSCRIPT_DIRECTION direction, SeqInterval rnaInterval, SeqInterval cdsInterval) {
        if (direction == Transcript.TRANSCRIPT_DIRECTION.TO_RIGHT) {
            if (cdsInterval.start() - rnaInterval.start() > 0) {return new SeqInterval(rnaInterval.start(), cdsInterval.start() - 1);}
            else {return new SeqInterval(rnaInterval.start(), rnaInterval.start());}                    
            }
        else {
            if (rnaInterval.end() - cdsInterval.end() > 0) {return new SeqInterval(cdsInterval.end() + 1, rnaInterval.end());}
            else {return new SeqInterval(rnaInterval.end(), rnaInterval.end());}   
        }
    }
    
    
    /**
     * Returns a SeqInterval giving the interval of the three prime UTR of a mRNA
     * used to add 3 prime UTR to Transcript instance, but that value of transcript is not used any more
     * @param direction  direction of transcription left or right
     * @param rnaInterval  the interval for the mRNA transcript
     * @param cdsInterval  the interval for the CDS of a mRNA
     * @return 
     */
    private SeqInterval getThreePrimeUtr (Transcript.TRANSCRIPT_DIRECTION direction, SeqInterval rnaInterval, SeqInterval cdsInterval) {
        if (direction == Transcript.TRANSCRIPT_DIRECTION.TO_RIGHT) {
            if (rnaInterval.end() - cdsInterval.end() > 0) {return new SeqInterval(cdsInterval.end() + 1, rnaInterval.end());}
            else {return new SeqInterval(rnaInterval.end(), rnaInterval.end());}                    
            }
        else {
            if (cdsInterval.start()- rnaInterval.start() > 0) {return new SeqInterval(rnaInterval.start(), cdsInterval.start() - 1);}
            else {return new SeqInterval(rnaInterval.start(), rnaInterval.start());}   
        }
    }
    
  
    
    /**
    * returns a list of transcribed blocks based on a list of Transcript
    * checked against a previous versions with more complex algorithms
    * Also checked a few results against the original Genbank file - seems robust
    * only difference was the older one separated block that were adjacent (no non-transcribed space between) 
    * @param transcriptList the list of Transcript
    * @return the list of transcribed blocks
    */
    public List<SeqInterval>  getTranscriptBlocks(List<Transcript> transcriptList) {
        //make an array of booleans indicating whether a postion in the sequence is in a transcript
        boolean[] isInTranscript = getIsInTranscriptArray(transcriptList);

        //make list of indices corresponding to starts and ends of trasncription blocks
        List<Integer> transcriptTransitions = new ArrayList<>();
        if (isInTranscript[0]) {
            transcriptTransitions.add(0);
        }
        for (int i = 1; i < ParameterSet.SEQUENCE_LENGTH; i++) {
            if (isInTranscript[i] && ! isInTranscript[i - 1]) {transcriptTransitions.add(i);}
            else if ( ! isInTranscript[i]  && isInTranscript[i - 1]) {transcriptTransitions.add(i - 1);}
        }
        if (isInTranscript[ParameterSet.SEQUENCE_LENGTH - 1]) {transcriptTransitions.add(ParameterSet.SEQUENCE_LENGTH - 1);}
        
        //make list of intervals corresponding to transcriptblocks
        List<SeqInterval> transcriptIntervals = new ArrayList<>();
            for (int i = 0; i < transcriptTransitions.size() - 1; i += 2) {
                transcriptIntervals.add(new SeqInterval(transcriptTransitions.get(i), transcriptTransitions.get(i + 1)));
            }
        return transcriptIntervals;
    }
    
    
    boolean[] getIsInTranscriptArray (List<Transcript> transcriptList) {
        //make an array of booleans indicating whether a postion in the sequence is in a transcript
        boolean[] isInTranscript = new boolean[ParameterSet.SEQUENCE_LENGTH];
        for (int i = 0; i < transcriptList.size(); i++) {
            for (int j = transcriptList.get(i).getTranscriptInterval().start(); 
                    j <= transcriptList.get(i).getTranscriptInterval().end(); j++) {
                isInTranscript[j] = true;
            }
        }
        return isInTranscript;
    }
    
    
    /**
     * Returns an array containing the number of nucleotides in each bin that are in transcript blocks
     * total should equal the bin size if completely within trasncript block
     * total should equal 0 if completely outside transcript block
     * total should be between 0 and the bin size for bins that have nucleotides both inside and outside transcript blocks 
     * @param transcriptList
     * @return 
     */
    public int[] getNumberTranscribedNucleotidesInPuBins(List<Transcript> transcriptList) {
        boolean[] isInTranscript = getIsInTranscriptArray(transcriptList);
        int[] numberTrueInBinsArray = new int[ParameterSet.NUMBER_BINS];
        for (int i = 0; i < ParameterSet.NUMBER_BINS; i++) {
             for (int j = 0; j < ParameterSet.PU_BIN_SIZE; j++) {
                if (isInTranscript[i * ParameterSet.PU_BIN_SIZE + j]) {
                    numberTrueInBinsArray[i]++;
                }
             }
        }
        return numberTrueInBinsArray;
    }
}
