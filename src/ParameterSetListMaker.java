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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Maker of a list of parameter sets based on a map of recipes
 *
 * @author tkelly
 */
public class ParameterSetListMaker {

    //These parameters are defined in the class DefaultParameterRecipes
    private final List<Integer> initiatorSiteLenList;
    private final List<Integer> numberPreRCsList;
    private final List<Integer> numberCellsList;
    private final List<Double> exponentialCoefficientList;
    private final List<Double> exponentialPowerFactorList;
    private final List<Double> maxFiringProbabilityPerMinList;
    private final List<Integer> elongationRateList;
    private final List<Double> timeConstantList;
    private final List<Double> minPerCycleList;
    private final List<Double> attenuationFactorList;
    private final List<Double> unusedParameterList;
    

    //List of parameter sets to be tested - all possible combinations of parameter values
    private final List<ParameterSet> parameterSetList;

    /**
     * A list of Parameter Sets from a set of recipes for iterating over all combinations of parameters in the recipes
     * @param parameterRecipeMap A map of the parameter recipes - key is the name of the parameter, value
     * is a String giving the recipe
     */
    public ParameterSetListMaker(Map<String, String> parameterRecipeMap) {
        //make a copy of the default recipes this is for later when 
        //these recipes may be overidden by values from file
        Map<String, String> recipeMap = new HashMap<>();
        parameterRecipeMap.keySet().stream().forEach(e -> recipeMap.put(e, parameterRecipeMap.get(e)));

        //make list of values for each parameter 
        initiatorSiteLenList = new ArrayList<>();
        parseRecipeInteger(recipeMap.get("initiatorSiteLength"), initiatorSiteLenList);
        numberPreRCsList = new ArrayList<>();
        parseRecipeInteger(recipeMap.get("numberPreRCs"), numberPreRCsList);
        numberCellsList = new ArrayList<>();
        parseRecipeInteger(recipeMap.get("numberCells"), numberCellsList);
        exponentialCoefficientList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("exponentialCoefficient"), exponentialCoefficientList);
        exponentialPowerFactorList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("exponentialPowerFactor"), exponentialPowerFactorList);
        maxFiringProbabilityPerMinList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("maxFiringProbabilityPerMin"), maxFiringProbabilityPerMinList);
        elongationRateList = new ArrayList<>();
        parseRecipeInteger(recipeMap.get("elongationRate"), elongationRateList);
        timeConstantList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("timeConstant"), timeConstantList);
        minPerCycleList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("minPerCycle"), minPerCycleList);
        attenuationFactorList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("attenuationFactor"), attenuationFactorList);
        unusedParameterList = new ArrayList<>();
        parseRecipeDouble(recipeMap.get("unusedParameter"), unusedParameterList);

        // Make list of parameter sets from all combination of parameter values
        parameterSetList = new ArrayList<>();
        initiatorSiteLenList.stream().
                forEach(a -> numberPreRCsList.stream().
                forEach(b -> numberCellsList.stream().
                forEach(c -> exponentialCoefficientList.stream().
                forEach(d -> exponentialPowerFactorList.stream().
                forEach(e -> maxFiringProbabilityPerMinList.stream().
                forEach(f -> elongationRateList.stream().
                forEach(g -> timeConstantList.stream().
                forEach(h -> minPerCycleList.stream().
                forEach(i -> attenuationFactorList.stream().
                forEach(j -> unusedParameterList.stream().
                forEach(k -> {
                    ParameterSet tempParameterSet = new ParameterSet(a, b, c, d, e, f, g, h, i, j, k);
                    parameterSetList.add(tempParameterSet);
                })))))))))));
    }

    /**
     * Parses a recipe containing a start value, number of values and increment
     * for values of a parameter all values must be type double adds each
     * parameter to a list of parameters called by constructor
     *
     * @param key the parameter name (and key to corresponding map of default
     * recipes)
     * @param recipe the recipe to be parsed - must be of form "double",double1,
     * double2,...
     * @param list the parameter list to be added to
     */
    private void parseRecipeDouble(String recipe, List<Double> list) {
        String[] tempArray = recipe.split(",");
        if (tempArray[0].equals("double")) {
            double start = Double.parseDouble(tempArray[1]);
            int number = Integer.parseInt(tempArray[2]);
            double increment = Double.parseDouble(tempArray[3]);
            for (int i = 0; i < number; i++) {
                list.add(start + i * increment);
            }
        } else {
            throw new IllegalArgumentException("parameter must be an double");
        }
    }

    /**
     * Parses a recipe containing a start value, number of values and increment
     * for values of a parameter all values must be type integer adds each
     * parameter to a list of parameters called by constructor
     *
     * @param key the parameter name (and key to corresponding map of default
     * recipes)
     * @param recipe the recipe to be parsed - must be of form
     * "integer",integer1, integer2,...
     * @param list the parameter list to be added to
     */
    private void parseRecipeInteger(String recipe, List<Integer> list) {
        String[] tempArray = recipe.split(",");
        if (tempArray[0].equals("int")) {
            int start = Integer.parseInt(tempArray[1]);
            int number = Integer.parseInt(tempArray[2]);
            int increment = Integer.parseInt(tempArray[3]);
            for (int i = 0; i < number; i++) {
                list.add(start + i * increment);
            }
        } else {
            throw new IllegalArgumentException("parameter must be an int");
        }
    }


    /**
     * returns deep copy (defensive) of list of parameter sets
     *
     * @return list of parameter sets with all combinations of parameters
     */
    public List<ParameterSet> getParameterSetList() {
        List<ParameterSet> tempList = new ArrayList<>();
        parameterSetList.stream().forEach(e -> tempList.add(e));
        return tempList;
    }
}
