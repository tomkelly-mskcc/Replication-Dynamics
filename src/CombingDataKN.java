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
import java.nio.file.DirectoryIteratorException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author tkelly
 * 12/24/2015
 * 
 * Reads single molecule data in Nurse format into an array of CombedMolecule Objects
 * each describing the many properties of a combed molecule
 * For example, each CombedMolecule object references an array of CombedMoleculeSegment
 * containing the individual segments within the molecule and their types
 * (G = replicated; N=not replicated). Each Molecule object also 
 * references  a CombedMoleculeProperties object that contains arrays of IODs,
 * Gs, Ns, and many other properties of the given molecule. See class 
 * CombedMoleculeProperties for details.
 * Class also contains a number of static convenience routines for printing lists (probably not optimal)
 * KN is Kaykov-Nurse
 */

public class CombingDataKN {
    
    private final  List<CombedMolecule> molList;  //the list of molecules in a polymerase usage experiment from files reformated by Kelly from Kaykov and Nurse et al
    
    /**
     * Creates an instance of DNACombingDataKN with list of CombedMolecules
     * @param dataFolder the folder containing the combing data files 
     */
    
    public CombingDataKN(String dataFolder) {
        
        //Make an array of all the measured molecules
        molList = new ArrayList<>();

        BufferedReader inputStream = null;
        String rawLine;
        List<Path> inputFileNames = new ArrayList<>();
        
        //Read file paths into array inputFileNames
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(Paths.get(dataFolder))) {
            for (Path file : stream) {
                inputFileNames.add(file);
            }
        } catch (IOException | DirectoryIteratorException x) {System.err.format("Problem getting data files: %s%n", x);}

        //Iterate over paths 
        for (Path p : inputFileNames) {
            //read contents of file
            try {
                inputStream = Files.newBufferedReader(p);
                for (int i = 0;(rawLine = inputStream.readLine()) != null; i++) {
                    String[] dataLine1 = rawLine.split(",");
                    String[] dataLine2 = inputStream.readLine().split(",");
                    String s = p.getFileName().toString().replaceAll(".csv","");
                    molList.add(new CombedMolecule(s, i+1, dataLine1, dataLine2));
                    int k =1;
                } 
            } catch (IOException x) {
                System.err.format("Problem reading data files: %s%n", x);
            }
        }
        
        //done - close inputStream
        try {
            if (inputStream != null) {inputStream.close();}
        } catch (IOException ex2) {
            System.err.format("IOException: %s%n", ex2);
        }
    }
    
    /**
     * Returns a list of CombedMolecule from Kaykov and Nurse data
     * Each instance of CombedMolecule describes properties of a combed molecule
     * @return combed molecule list
     */
    public List<CombedMolecule> getMoleculeArray() {
        return molList;
    }
}
