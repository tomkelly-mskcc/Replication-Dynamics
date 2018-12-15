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
 * only found three major effects of Rif1 and 3 minor so put into a static method
 * @author tkelly
 */
public class Rif1Sites {
    private static final List<SeqInterval> RIF1_INTERVAL_LIST;
    static{
    RIF1_INTERVAL_LIST = new ArrayList<>();
    RIF1_INTERVAL_LIST.add(new SeqInterval(888876, 891984));//major ranks 10
    RIF1_INTERVAL_LIST.add(new SeqInterval(948016, 949125));//significant local ranks 2
    RIF1_INTERVAL_LIST.add(new SeqInterval(2337536, 2339983));//major effect ranks 7
    RIF1_INTERVAL_LIST.add(new SeqInterval(4128445, 4130584));//minor rank 6
    RIF1_INTERVAL_LIST.add(new SeqInterval(3378647, 3380700));//major ranks 1 documented in Masai paper
    RIF1_INTERVAL_LIST.add(new SeqInterval(1079379, 1080184));//small effect rank 9
    }
 
    /**
     * Returns list of 6 sites where Rif1 makes a difference
     * @return Rif1 site list
     */
    static List<SeqInterval> getRIF1_INTERVAL_LIST() {
        return RIF1_INTERVAL_LIST;
    }

    
    
    
    
    
    
    
    
}
