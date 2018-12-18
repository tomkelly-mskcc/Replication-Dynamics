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

import java.util.HashMap;
import java.util.Map;

/**
 * The default parameter recipes - maybe later implement  a file of recipes
 * Each recipe specifies a range of values for a parameter to be tested during parameter optimization
 * Each recipe has two strings.  the first is the name of the parameter
 * and the second specifies a datatype,a starting value of parameter, the number of parameters to be scanned, and the increment
 * @author tkelly
 */
public class DefaultParameterRecipes {
    

     //Initialize default parameter recipes
    private final static Map<String, String> DEFAULT_PARAMETER_RECIPES = new HashMap<>();
    static {
        DEFAULT_PARAMETER_RECIPES.put("initiatorSiteLength", "int,25,1,5"); //length of window for initiator binding
        DEFAULT_PARAMETER_RECIPES.put("numberPreRCs", "int,363,1,50"); //363 for chromosome 2 1040 for genome
        DEFAULT_PARAMETER_RECIPES.put("numberCells", "int,1000,1,1000"); //size of cell poulation -number or replicated molecules
        DEFAULT_PARAMETER_RECIPES.put("exponentialCoefficient","double,9e-11,1,1"); //DEPRECATED
        DEFAULT_PARAMETER_RECIPES.put("exponentialPowerFactor","double,21,1,1");// Now constant in exponential function for AT content
        DEFAULT_PARAMETER_RECIPES.put("maxFiringProbabilityPerMin", "double,0.3,1,.1"); // maximum firing rate - prev optimimization 0.22  .3
        DEFAULT_PARAMETER_RECIPES.put("elongationRate", "int,2000,1,100"); //fork velocity in nucleotides per min
        DEFAULT_PARAMETER_RECIPES.put("timeConstant", "double,0.022,1,.002"); //This is now the rate of increase of the firing rate per minute 
        DEFAULT_PARAMETER_RECIPES.put("minPerCycle", "double,0.01,1,0.1"); //replication cycle time in min
        DEFAULT_PARAMETER_RECIPES.put("attenuationFactor", "double,0.0,1,.05"); //factor for sites in transcripts
        DEFAULT_PARAMETER_RECIPES.put("unusedParameter", "double,0.05,1,0.1"); //for future use
    }

    public DefaultParameterRecipes() {
    }
    
    /**
     * Returns the map of parameters - name maps to recipe
     * @return defensive copy of default parameters  Strings immutable
     */
    public static Map<String, String> getDEFAULT_PARAMETER_RECIPES() {
        Map<String, String> tempMap = new HashMap<>();
        DEFAULT_PARAMETER_RECIPES.keySet().stream().forEach(e-> tempMap.put(e, DEFAULT_PARAMETER_RECIPES.get(e)));
        return tempMap;
    }
}
