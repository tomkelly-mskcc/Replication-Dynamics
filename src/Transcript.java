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
 * A transcript - value class holding properties of a transcript
 * @author tkelly
 */
public class Transcript {
    private final SeqInterval transcriptInterval; //A genome segment corresponding to a transcript
    private final TRANSCRIPT_TYPE transcriptType;  //The type of transcript - see enum below
    private final TRANSCRIPT_DIRECTION transcriptDirection;  //The direction of a transcript - see enum below
    private final String description;  //The description of the transcript
    private final int numberIntrons; //the number of introns in the primary transcript
    private final SeqInterval fivePrimeUtr; //The five prime untranslated region of the transcript
    private final SeqInterval cdsblock; //the coding sequence of the transcript
    private final SeqInterval threePrimeUtr;  //The three prime untranslated region of the transcript
    
    
    /**
     * Constructs the transcript from instance of SequenceDNAAnnotationGenbank
     * @param transcriptInterval  //The sequence interval of the genome encompassing the transcript
     * @param type //The type of transcript
     * @param direction  //The direction of the transcript
     * @param description //The description of the transcript
     * @param numbIntrons //The number of introns in the transcript
     */
    public Transcript(SeqInterval transcriptInterval, TRANSCRIPT_TYPE type, TRANSCRIPT_DIRECTION direction, String description, int numbIntrons) {
        this.transcriptInterval = transcriptInterval;
        transcriptType = type;
        this.description = description;
        transcriptDirection = direction;
        numberIntrons = numbIntrons;
        fivePrimeUtr = null;
        threePrimeUtr = null;
        cdsblock = null;
    }
        
    /**
    * Alternate constructor when more information is available
    * @param transcriptInterval  //The sequence interval of the genome encompassing the transcript
    * @param type //The type of transcript
    * @param direction  //The direction of the transcript
    * @param description //The description of the transcript
    * @param numbIntrons //The number of introns in the transcript
    * @param fiveUTR //The five prime untranslated region
    * @param cds //the coding region
    * @param threeUTR //The three prime untranslated region
    */
    public Transcript(SeqInterval transcriptInterval, TRANSCRIPT_TYPE type, TRANSCRIPT_DIRECTION direction, String description, int numbIntrons,
            SeqInterval fiveUTR, SeqInterval cds, SeqInterval threeUTR) {
        this.transcriptInterval = transcriptInterval;
        transcriptType = type;
        this.description = description;
        transcriptDirection = direction;
        numberIntrons = numbIntrons;
        fivePrimeUtr = fiveUTR;
        cdsblock = cds;
        threePrimeUtr = threeUTR;
    }
    
    
    /**
     * The transcript type
     */
    public static enum TRANSCRIPT_TYPE {
        MRNA,
        SNRNA,
        TRNA,
        RRNA,
        SNORNA,
        PSEUDOGENE,
        NCRNA
    }
    
    /**
     * The transcript direction
     */
    public static enum TRANSCRIPT_DIRECTION {
        TO_RIGHT,
        TO_LEFT
    }

    /**
     * REturns the interval in genome (as SeqInterval) of transcript
     * @return the transcript interval in genome
     */
    SeqInterval getTranscriptInterval() {
        return transcriptInterval;
    }

    /**
     * Returns the type of the transcript
     * @return type of transcript
     */
    TRANSCRIPT_TYPE getTranscriptType() {
        return transcriptType;
    }

    /**
     * Returns the description of the transcript
     * @return description of transcript
     */
    String getDescription() {
        return description;
    }

    /**
     * Returns the number of introns in the transcript
     * @return number of introns  in transcript
     */
    int getNumberIntrons() {
        return numberIntrons;
    }

    /**
     * Returns the 5' UTR of the transcript
     * @return 5' UTR of transcript
     */
    SeqInterval getFivePrimeUtr() {
        return fivePrimeUtr;
    }

    /**
     * Returns the 3' UTR of the transcript
     * @return 3' UTR of transcript
     */
    SeqInterval getThreePrimeUtr() {
        return threePrimeUtr;
    }

    /**
     * Returns the codon block of the transcript
     * @return codon block of transcript
     */
    SeqInterval getCdsblock() {
        return cdsblock;
    }
    
    


    @Override
    public String toString() {
        if (transcriptType == TRANSCRIPT_TYPE.MRNA) {
            return String.format("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
                    transcriptType, transcriptDirection, transcriptInterval.start(), 
                    transcriptInterval.end(),description, numberIntrons, 
                    fivePrimeUtr.start(), fivePrimeUtr.end(), cdsblock.start(), cdsblock.end(),
                    threePrimeUtr.start(), threePrimeUtr.end());
        } else
        return String.format("%s\t%s\t%d\t%d\t%s\t%d\t", transcriptType, 
                transcriptDirection, transcriptInterval.start(), transcriptInterval.end(), 
                description, numberIntrons);
    }
    
    public void printMrnaComponents () {
        if (fivePrimeUtr != null) {System.out.printf("%d\t%d\t",fivePrimeUtr.start(), fivePrimeUtr.end());}
        if (cdsblock != null) {System.out.printf("%d\t%d\t",cdsblock.start(), cdsblock.end());}
        if (threePrimeUtr != null) {System.out.printf("%d\t%d\n",threePrimeUtr.start(), threePrimeUtr.end());}
    }
    
    
}
