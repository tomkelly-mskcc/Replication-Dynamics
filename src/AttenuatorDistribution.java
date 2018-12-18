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
 * Creates a list of attenuators at for each initiator binding site in genome
 * Input is a list of transcribed blocks
 * initiator binding sites are indexed by their first nucleotide
 * the indices run from 0 to seqLen - initSiteLen + 1
 * @author tkelly
 */
public class AttenuatorDistribution {
    
    private final List<Double> attenuatorDistribution;

    /**
     * A distribution of attenuation factors due to transcription over the genome
     * @param blockList list of SeqInterval describing transcriptionBlocks
     * @param parameters parameters for replication
     * @param seqLength length of the sequence
     */
    public AttenuatorDistribution(List<SeqInterval>  blockList, ParameterSet parameters, int seqLength) {
        
        int siteLength = parameters.getInitiatorSiteLength();
        
        //Initialize attenuator array with all elements equal to 1.0
        List<Double> tempAttenuatorDistribution = new ArrayList<>();
        for (int i = 0; i < seqLength - parameters.getInitiatorSiteLength() + 1; i++) {
            tempAttenuatorDistribution.add(1.0);
        }
        
        //set elements in attenuator distribution to attenuation factor if they overlap transcription blocks
        int startIndex;
        int endIndex;
        for (int i = 0; i < blockList.size(); i++) {
            //if the start of a transcription block is between 0 and siteLength - 1, set startIndex at 0
            if (blockList.get(i).start() - siteLength < 0) {
                startIndex = 0;
            //so most of the time startIndex is set at start of transcription block is set at startIndex - sitelength
            } else {
                startIndex = blockList.get(i).start() - siteLength;
            }
            
            // if the end of the transcription block is between seqLen -initSiteLen and seqLen -1 then set the endIndex at seqLen - initSiteLen 
            if (blockList.get(i).end() >= tempAttenuatorDistribution.size()) {
                endIndex = tempAttenuatorDistribution.size() - 1;
            //so most of the time we set endIndex at end of block list
            } else {
                endIndex = blockList.get(i).end();
            }
            //from startIndex to endIndex  we set the attnuator distribution to attenuator factor
            //if left end of block falls anywhere in this block attenuator holds because therere is overlap of the initSite with a transcription block
            //this depends on the fact that we are indexing the initSite with the first nucleotide
            for (int j = startIndex; j <= endIndex; j++) {
                tempAttenuatorDistribution.set(j, parameters.getAttenuationFactor());
            }
        }
        attenuatorDistribution = tempAttenuatorDistribution;
    }
    

    /**
     * Returns the list of attenuator values for positions in sequence
     * @return attenuatorDistribution
     */
    List<Double> getAttenuatorDistribution() {
        return attenuatorDistribution;
    }

    @Override
    public String toString () {
        StringBuilder builder = new StringBuilder();
        attenuatorDistribution.stream().forEach(e-> {builder.append(e).append('\n');});
        return builder.toString();
    }
    
    
}
