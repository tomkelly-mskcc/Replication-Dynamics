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
 * The main DNA synthesis engine
 * Each cycle of synthesis is invoked by start() routine
 * For each short interval specified as cycle time in ParameterSet, the engine
 * first initiates replication at sites with probability given in argument 
 * and then elongates all forks by an amount given in argument 
 * The engine updates the states of the potential initiation sites in the potentialInitSiteList
 * 
 * could probably make this code a bit cleaner, but it has been thoroughly tested and it works and is fast
 * it is very complex and beautifying it would only risk introducing errors
 * @author tkelly
 */
public class DnaSynthesisCycle2 {
    
    //Potential Initiation sites
    private final List<PotentialInitiationSite2> potentialInitSiteList;  //Listwhose potentialintiation sites are to be updated
    private final PotentialInitiationSite2 leftEnd; //potential initiation site  at left end of the molecule(no forks)
    private final PotentialInitiationSite2 rightEnd;  //potential initiation site at the right end of the molecule (no forks)
    private PotentialInitiationSite2 lastActivePIS;  // the last potential initiation site visited
    
    private int leftEndFork; //fork moving toward left end
    private int rightEndFork; //fork moving toward right end
    private int forkOnLeft; //fork immediately to the left of this site
    private int forkOnRight;  //fork immediately to the right of this site
    
    // summary numbers for states
    private int terminations = 0;
    private int potentials;
    private int passives = 0;
    private int actives = 0;
    private int initiations = 0;
    private int closures = 0;  //these are events when two forks meet
    private int numberNucleotidesReplicated = 0;
    private double fractionReplicated;
    private int forks = 0;
    private final int seqLength;

    /**
     * Create a DNA synthesis cycle
     * @param potentialSiteList  the list of postentialInitiationSite (which defines the state of a site)
     * @param leftEnd potential initiation site at left end _ not in list - no forks
     * @param rightEnd potential initiation site at right end - not in list - no forks
     */
    public DnaSynthesisCycle2 (List<PotentialInitiationSite2> potentialSiteList, PotentialInitiationSite2 leftEnd, PotentialInitiationSite2 rightEnd) {
        this.leftEnd = leftEnd;  //check
        this.rightEnd = rightEnd;
        seqLength = rightEnd.getPosition();
        this.potentialInitSiteList = potentialSiteList;
        potentials = potentialSiteList.size();
    }
    
    /**
     * Implements the two steps of a single cycle of DNA synthesis
     * updates the state of the genome by calls to methods in PotentialInitiationSite.
     * @param firingProbabilityForCycle
     * @param nucleotidesOfForkMovementForCycle 
     */
    public void start(double firingProbabilityForCycle, int nucleotidesOfForkMovementForCycle) {
        
        /**
         *  INITIATION ENGINE
         */
        lastActivePIS = leftEnd;
        for (PotentialInitiationSite2 site : potentialInitSiteList) { //INITIATION ENGINE  FIX Maybe ISSUE WITH ORDER
            //System.out.println("");
            if (site.IsActive()) {
                lastActivePIS = site;
            }
            else if (site.IsPotential()) {
                if (Math.random() < firingProbabilityForCycle) {                    
                    site.activate(lastActivePIS);
                    //System.out.println(site.getPosition());
                    numberNucleotidesReplicated++;  //result of placement of fork - initiation site is replicated by one of them
                    lastActivePIS = site;
                    initiations++;
                    actives++;
                    potentials--;
                    forks += 2;
                }
            }
        }        
        /**
        *  SYNTHESIS ENGINE
        * 
        * Left end is special case - leftmost active site has left-moving fork until it terminates.
        * The test for termination must take into account that there is no active fork at the left end of the chromosome.
        * Similarly, the number of nucleotides synthesized is increased by the activity of only a single fork.
        * Once a left-moving fork from a leftmost active site terminates, the leftmost active site has only a rightward-moving fork.
        * therefore we need to keep track of the position of a left moving fork until it terminates in order to:
        * 1) apply the correct test for termination, 2)correctly increment the number of nucleotides synthesized, and 
        * 3) mainly to determine whether Potential sites in the region to the left of moving fork are passively replicated.
        */
        if (initiations == 0) {return;}  //first initiation hasn't happened
        
        //DEAL WITH LEFT FORK OF LEFTMOST ACTIVE SITE
        PotentialInitiationSite2 leftSite = leftEnd.getRightActivePIS();  //get leftmost active site
        
        //IF LEFTFORK OF LEFTMOSTSITE IS ACTIVE (if not we do nothing)
        if (leftSite.IsLeftForkActive()) { 
            //TWO CASES will the left fork terminate at the left end or will it not terminate
            //if so, 
            // if right fork was already inactivated, terminate the iste
            
            //CASE I -LEFT FORK TERMINATES - properly update number of nucleotides replicated, set the position of left fork to 0
            //make a note of the position of the left end fork to check for passive repliction
            //set left fork to inactive
            //update number of forks and closures
            //if this inactivates the leftmostPIS - we catch it later
            if ((leftSite.getPositionLeftFork() - nucleotidesOfForkMovementForCycle) <= 0) {  //test for termination of leftmost fork on this cycle
                numberNucleotidesReplicated += leftSite.getPositionLeftFork();  
                leftSite.setPositionLeftFork((int) 0);
                leftEndFork = 0; // note position of left end fork so use to check for passive replication
                leftSite.inactivateLeftFork();
                forks--;
                closures++;
                
                //TERMINATE AND UPDATE LINKAGES
                if (!leftSite.IsRightForkActive()) {  
                    leftSite.terminate();  //Set flags and update linkages
                    terminations++;
                    actives--;
                }
                
               
            //CASE 2 - LEFT FORK DOES NOT TERMINATE
            //if left fork from leftmost site not going to terminate at aleft end then extend fork to left
            //note the position of left most fork to check for passive replication
            //properly increment number of nucleotides synthesized
            } else {
                leftSite.extendLeftFork(nucleotidesOfForkMovementForCycle);
                leftEndFork = leftSite.getPositionLeftFork();  //note position of left fork to determine whether potential sites have been passively replicated
                numberNucleotidesReplicated += nucleotidesOfForkMovementForCycle;
            }
              //CHECK FOR PASSIVE REPLICATION IN LEFT INTERVAL  Not necessary will be picked up in next loop
                /*for (int i = leftSite.getIndex() -1 ; i >= 0; i--) {
                    if (potentialInitSiteList.get(i).getPosition() >= leftEndFork) {
                        potentialInitSiteList.get(i).passivelyReplicate();
                    }
                    break;
                }*/
                
        }
        
        //MAIN LOOP ITERATE ACROSS PIS 
        //UPDATE POSITIONS OF RIGHT FORK AND LEFT FORK OF NEXT ACTIVE
        for (PotentialInitiationSite2 site : potentialInitSiteList) { 
            //if site is not active skip but check thatit doesn't have active forks
            if (!site.IsActive()) {  //don't need this as only updating actives and left and righ ends probably innocuous
                if (site.IsLeftForkActive() || site.IsRightForkActive()) {
                    System.out.println("inactive site with active forks");
                }
            }
            //PROCESS ALL ACTIVE SITES
            if (site.IsActive()) {
                site.linkageCheck(seqLength);  //FIX maybe eliminate in final
                
                //CASE I - RIGHTFORK IS ACTIVE AND IS THE RIGHTMOST ACTIVE SITE
                //TWO SUBCASES - RIGHT FORK WILL TERMINATE OR RIGHT FORK WILL NOT TERMINATE
                if (site.IsRightForkActive()) {  //right fork is active
                    if (site.getRightActivePIS().equals(rightEnd)) { //Site is the rightmost -
                       
                        //CASE1A -RIGHT FORK OF RIGHTMOST ACTIVE SITE WILL TERMINATE AT RIGHT END    
                        if ((site.getPositionRightFork() + nucleotidesOfForkMovementForCycle) >= seqLength - 1) {  //if it terminates
                            site.inactivateRightFork();
                            forks--;
                            numberNucleotidesReplicated += (seqLength - 1 - site.getPositionRightFork());
                            site.setPositionRightFork(seqLength - 1);
                            rightEndFork = seqLength - 1;  //for passive identification
                            closures++;

                            //terminate rightmost site if both forks inactivated - this is necessary because if it needs to be terminated
                            //it would not have been terminated during processing of previous interval (it's not an endpoint of previousinterval
                            if (!site.IsLeftForkActive()) {  //if not, this site will not have been terminated 
                                site.terminate();
                                terminations++;
                                actives--;
                            }
                            
                        //CASE IB - RIGHT FORK OF RIGHTMOST ACTIVE SITE WILL NOT TERMINATE TERMINATE
                        //Extend right fork, save right fork psoition for dealing with passives, update nucs
                        } else {
                            //rightmost fork doesn't terminate
                            site.extendRightFork(nucleotidesOfForkMovementForCycle);
                            rightEndFork = site.getPositionRightFork();  // save for passive replication test on potentils
                            numberNucleotidesReplicated += nucleotidesOfForkMovementForCycle;
                        }
                        
                    //CASE II - RIGHT FORK IS ACTIVE AND SITE IS NOT THE RIGHTMOST
                    //TWO SUBCASE RIGHT AND LEFT FORKS (FROM RIGHT ACTIVE SITE) WILL MEET OR RIGHT AND LEFT FORKS WILL NOT MEET
                    } else { //Site is active and not the rightmost
                        
                        //CASE IIA - RIGHT FORK WILL MEET LEFT FORK FROM RIGHT ACTIVE SITE
                        //inactivate both forks
                        //update positions of both forks
                        //save positions of forks for dealing passives
                        //terminate either or both sites depending on the state of their other forks
                        //update linkages of if either site terminates 
                        //update nucs and update flags
                        if (IsTerminatedOnRight(site,nucleotidesOfForkMovementForCycle)) {  //if forks will meet
                            site.inactivateRightFork(); //inactivate right fork
                            site.getRightActivePIS().inactivateLeftFork(); //inactivate left fork of right active
                            forks -= 2;
                            int spaceBetweenForks = site.getRightActivePIS().getPositionLeftFork() - site.getPositionRightFork() - 1;  //changed long to int
                            if (spaceBetweenForks == 1) {site.setPositionRightFork(site.getPositionRightFork() + 1);}  //only one nucleotide beteen forks just update right fork position
                            if (spaceBetweenForks >=2) {  //more than one nucleotide between forks - change positions of both
                                site.setPositionRightFork(site.getPositionRightFork() + spaceBetweenForks / 2);
                                site.getRightActivePIS().setPositionLeftFork(site.getPositionRightFork() + 1);
                            }
                            forkOnLeft = site.getPositionRightFork();  //for dealing with passives
                            forkOnRight = site.getRightActivePIS().getPositionLeftFork(); //for dealing with passives
                            //FIX - PROBABLY DON'T NEED THIS
                            if (forkOnRight - forkOnLeft != 1) {  //for debugging detects positions of forks not updated correctly after they met
                                System.out.println("forks should have met" + (forkOnRight - forkOnLeft));
                            }
                            numberNucleotidesReplicated += spaceBetweenForks;
                            closures++;
                            if (!site.IsLeftForkActive()) {  //The rule is; any time a fork at a site is terminated check whether the other fork is also terminated
                                site.terminate(); //SHOULD UPDATE LINKAGES OK - now left active pis points to right active pis
                                terminations++;
                                actives--;
                            }
                            if (!site.getRightActivePIS().IsRightForkActive()) {  
                                if (!site.getRightActivePIS().IsActive()) {  //for debugging
                                    System.out.println("redundant termination");
                                }
                                site.getRightActivePIS().terminate();  //now active site to left of site points to activesite to right of right active site
                                terminations++;
                                actives--;
                            }
                            
                        //CASE 2B -RIGHT FORK WILL NOT MEET LEFT FORK FROM RIGHT ACTIVE SITE
                        //Extend both forks
                        //save position of forks for dealing with passives
                        //update nucs
                        } else {
                            site.extendRightFork(nucleotidesOfForkMovementForCycle);
                            forkOnLeft = site.getPositionRightFork();
                            site.getRightActivePIS().extendLeftFork(nucleotidesOfForkMovementForCycle);
                            forkOnRight = site.getRightActivePIS().getPositionLeftFork();
                            numberNucleotidesReplicated += 2 * nucleotidesOfForkMovementForCycle;
                        }
                    }
                }
            
            
            //DETERMINE WHETHER POTENTIAL SITE WAS ACTIVELY REPLICATED AND UPDaTE ACCORDINGLY
            } else if (site.IsPotential()) {
                
                //CASE I - IT IS PASSIVELY REPLICATED
                if(IsPassivelyReplicated(site)) {
                    site.passivelyReplicate();
                    passives++;
                    potentials--;
                }
                
                
                else {
                    if (initiations >0) {
                        int i = 1;
                    }
                }
            }
        }
        fractionReplicated = (double) numberNucleotidesReplicated/seqLength;
    }
    
    /**
     * Returns true if forks will close on the right of the site
     * @param site current potential initiation site
     * @param nucleotidesToMoveForks number of nucleotides for fork movement on this cycle
     * @return true if closure on right
     */
    private boolean IsTerminatedOnRight(PotentialInitiationSite2 site, int nucleotidesToMoveForks) {
        return (site.getRightActivePIS().getPositionLeftFork() - nucleotidesToMoveForks) - (site.getPositionRightFork() + nucleotidesToMoveForks) <= 1;
    }
    
    /**
     * Returns true if site is passively replicated
     * @param site current potential initiation site
     * @return true if passively replicated site
     */
    private boolean IsPassivelyReplicated(PotentialInitiationSite2 site) {
        if (site.getPosition() < leftEnd.getRightActivePIS().getPosition()) {  //if in left interval before the first active
            if (site.getPosition() < leftEndFork) {return false;}
            //return true;
        }
        else if (site.getPosition() > rightEnd.getLeftActivePIS().getPosition()) { //if in right interval after the last active
            if (site.getPosition() > rightEndFork) {return false;}
            //else {return true;}
        }
        else if (site.getPosition() > forkOnLeft && site.getPosition() < forkOnRight) {  //if in other interval
            return false;
        }
        return true;
    }

    /**
     * Returns the number of active potential initiation sites
     * @return number of active sites
     */
    int getActives() {
        return actives;
    }

    /**
     * Returns the number of closures that have occurred
     * @return number of closures
     */
    int getClosures() {
        return closures;
    }

    /**
     * Returns the number of active forks
     * @return number of forks
     */
    int getForks() {
        return forks;
    }

    /**
     * Returns the fraction replicated
     * @return fraction replicated
     */
    double getFractionReplicated() {
        return fractionReplicated;
    }

    /**
     * Returns number of initiation events
     * @return number initiations
     */
    int getInitiations() {
        return initiations;
    }

    /**
     * Returns number of passive replication events
     * @return number passive replications
     */
    int getPassives() {
        return passives;
    }

    /**
     * Returns the number of potential initiation sites
     * @return number of potential sites
     */
    int getPotentials() {
        return potentials;
    }

    /**
     * Returns the number of termination events - a potential initiation site closes both of its forks 
     * @return number terminations
     */
    int getTerminations() {
        return terminations;
    }

    /**
     * Returns the total number of nucleotides replicated
     * @return number of nucleotides replicated
     */
    int getNumberNucleotidesReplicated() {
        return numberNucleotidesReplicated;
    }

    /**
     * Returns the last potential initiation site visited
     * @return last potential site
     */
    PotentialInitiationSite2 getLastActivePIS() {
        return lastActivePIS;
    }

    /**
     * Returns left end site
     * @return left end site
     */
    PotentialInitiationSite2 getLeftEnd() {
        return leftEnd;
    }

    /**
     * Returns right end site
     * @return right end site
     */
    PotentialInitiationSite2 getRightEnd() {
        return rightEnd;
    }
    
    /**
     * Returns the list of potential initiation sites for this replicating molecule
     * @return list of sites
     */
    List<PotentialInitiationSite2> getPotentialInitSiteList() {
        return potentialInitSiteList;
    }

    
    
    
    
    
    
    
    
    
    boolean potentialFlag = false;
    List<Integer> unrepList;
    public void consistencyCheck() {
        
        List<Integer> unrepBlockList = makeUnreplicatedBlockList();
        unrepList = unrepBlockList;
        potentialInitSiteList.stream().forEach(s -> {
            if(s.IsActive()) {
                if ((!s.IsLeftForkActive() && !s.IsRightForkActive()) || s.IsPotential() || s.IsTerminated() || s.IsPassivelyReplicated()) {System.out.println("error active site");}
            }
            else if (s.IsPotential() || s.IsPassivelyReplicated()) {  //must be in unreplicated block if potential or replicated block if passive
                for (int l = 0; l < unrepBlockList.size()-1; l += 2) {
                    if (s.getPosition() >=  unrepBlockList.get(l) && s.getPosition() <=  unrepBlockList.get(l+1)) {
                        potentialFlag = true;  //its in replicated block
                        break;
                    }
                }
                if (s.IsPotential()) {
                    if (!potentialFlag || s.IsActive() || s.IsTerminated() || s.IsPassivelyReplicated()) {
                        System.out.println("error potential site");
                    }
                }
            }
            else if (s.IsPassivelyReplicated()) {
                    if (potentialFlag || s.IsActive() || s.IsTerminated() || s.IsPotential()) {System.out.println("error passive site");}
                }
            potentialFlag = false;
        });
    }
    
    private int lBound = 0;
    private int rBound = 0;
    public List<Integer> makeUnreplicatedBlockList() {
        List<Integer> unreplicatedBlockList = new ArrayList<>();
        PotentialInitiationSite2 leftMostActiveSite = leftEnd.getRightActivePIS();
        if (!leftMostActiveSite.equals(rightEnd)) {
            if(leftMostActiveSite.IsLeftForkActive()) {
                unreplicatedBlockList.add((int) 1);
                unreplicatedBlockList.add(leftMostActiveSite.getPositionLeftFork()-1);
                System.out.format("%d   %d%n", 1, leftMostActiveSite.getPositionLeftFork()-1);
            }
        }
               
        potentialInitSiteList.stream().forEach(s -> {
            if (s.IsRightForkActive() && !s.IsLeftForkActive()) {lBound = s.getPositionRightFork()+1;}
            
            else if (s.IsLeftForkActive()  && !s.equals(leftMostActiveSite) && !s.IsRightForkActive()) {
                rBound = s.getPositionLeftFork() - 1;
                unreplicatedBlockList.add(lBound);
                unreplicatedBlockList.add(rBound);
                System.out.format("%d   %d%n", lBound, rBound);
            }
            else if (s.IsLeftForkActive() && s.IsRightForkActive()) {
                if (!s.equals(leftMostActiveSite)) {
                    rBound = s.getPositionLeftFork() - 1;
                    unreplicatedBlockList.add(lBound);
                    unreplicatedBlockList.add(rBound);
                    System.out.format("%d   %d%n", lBound, rBound);
                }
                lBound = s.getPositionRightFork() + 1;
            }
        });
        PotentialInitiationSite2 rightMostActiveSite = rightEnd.getLeftActivePIS();
        
            if (rightMostActiveSite.IsRightForkActive()) {
                unreplicatedBlockList.add(rightMostActiveSite.getPositionRightFork()+1);
                unreplicatedBlockList.add(seqLength);
                System.out.format("%d   %d%n", rightMostActiveSite.getPositionRightFork()+1, seqLength);
            }
        return unreplicatedBlockList;
    }
    
    
    
}

