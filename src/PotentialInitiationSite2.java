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

/**Mutable
 *A potential initiation site where a preRC may be assembled and activated -
 * the state of these sites is updated continuously as replication proceeds.
 * Each site has a state defined by: potential, active, terminated, passively replicated.
 * During replication the state may be altered.
 * In the beginning all potential initiation sites are in potential state (obviously).
 * In the end the state must be terminated or passively replicated.
 * Note that there are linkages among the active states, so that each active state
 * knows the references to the active states to the left and right.
 * 
 * @author tkelly
 */

public class PotentialInitiationSite2 {
    
    private final int position;  //position in genome/chromosome
    private final int index; //index of site in molecule
    private boolean active = false; 
    private boolean passivelyReplicated = false;
    private boolean leftForkActive = false; //is leftward moving fork from this site active
    private boolean rightForkActive = false; //is rightward moving fork from this site active
    private boolean terminated = false; //both forks no longer active
    private boolean potential = true;
    private int positionRightFork; //the position of the rightward moving fork
    private int positionLeftFork; //the position of the leftward moving fork
    private PotentialInitiationSite2 leftActivePIS = null;  //Active site to the left
    private PotentialInitiationSite2 rightActivePIS = null;  //Active site to the right

    
    /**
     * Constructor initializes the position and index of the potential initiation site from argument
     * @param position
     * @param index
     */
    public PotentialInitiationSite2 (int position, int index) {
        this.position = position;
        this.positionLeftFork = position;
        this.positionRightFork = position;
        this.index = index;
    }
    
    /**
     * 
     * SETTERS - CHANGE STATE
     */
    /**
     * Invoked when a PotentialIS changes state from Potential to Active.
     * Updates the state flags.  Also updates linkages among active PotentialISs.
     * 
     * @param lastActivePIS
     */
    void activate (PotentialInitiationSite2 lastActivePIS) {
        active = true;
        potential = false;
        leftForkActive = true;
        rightForkActive = true;
        updateLinkagesToAdd(lastActivePIS);
    }
    
    /**
     * Activates without setting linkages
     */
    void activateTerminus () {  
        active = true;
        potential = false;
        leftForkActive = true;
        rightForkActive = true;
    }
    
    /**
     *  inactivate the leftward moving fork
     */
    void inactivateLeftFork() {
        leftForkActive = false;
    }
    
    /**
     * inactivate the rightward moving fork
     */
    void inactivateRightFork() {
        rightForkActive = false;
    }
    
    /**
     * terminate the site - active to terminated
     */
    void terminate() { 
        active = false;
        terminated = true;
        updateLinkagesToTerminate();
    }
    
    /**
     * changes state to passively replicated from potential
     */
    void passivelyReplicate() {
        passivelyReplicated = true;
        potential = false;
    }
    
    /**
     * extends leftward moving fork
     * @param nucsToMove the number of nucleotides to move the fork
     */
    void extendLeftFork(int nucsToMove) {positionLeftFork -= nucsToMove;}
    
    /**
     * extends rightward moving fork
     * @param nucsToMove the number of nucleotides to move the fork
     */
    void extendRightFork(int nucsToMove) {positionRightFork += nucsToMove;}
    
    /**
     * initialize the right active PIS as right end PIS
     * @param rightPIS the right end PIS
     */
    void setRightActivePIS (PotentialInitiationSite2 rightPIS) {rightActivePIS = rightPIS;}
    
    /**
     * initialize the left active PIS as left end PIS
     * @param leftPIS left end PIS
     */
    void setLeftActivePIS (PotentialInitiationSite2 leftPIS) {leftActivePIS = leftPIS;}
    
    /**
     * initialize position of rightward moving fork
     * @param pos the position of the site
     */
    void setPositionRightFork(int pos) {positionRightFork = pos;}
    
    /**
     * initialize position of leftward moving fork
     * @param pos the position of the site
     */
    void setPositionLeftFork(int pos) {positionLeftFork = pos;}
        
    /**
     * set rightward moving fork to active
     * @param rightForkFlag right fork active flag
     */
    void setRightForkActive(boolean rightForkFlag) {rightForkActive = rightForkFlag;}
    
    /**
     * set leftward moving fork to active
     * @param leftForkFlag left fork active flag
     */
    void setLeftForkActive(boolean leftForkFlag) {leftForkActive = leftForkFlag;}

    /**
     * set site to active
     * @param active active flag
     */
    void setActive(boolean active) {
        this.active = active;
    }

    /**
     * set ste to potential
     * @param potential potential flag
     */
    void setPotential(boolean potential) {
        this.potential = potential;
    }

    /**
     * set site to passively replicated
     * @param passivelyReplicated passive flag
     */
    void setPassivelyReplicated(boolean passivelyReplicated) {
        this.passivelyReplicated = passivelyReplicated;
    }

    /**
     * set site to terminated
     * @param terminated terminated flag
     */
    void setTerminated(boolean terminated) {
        this.terminated =terminated;
        }
    
    /**
     * creates a string with key state variables
     * @return 
     */
    @Override
    public String toString() {
        String type;
        if (active && !potential && !terminated && !passivelyReplicated) {type = "active";}
        else if (!active && potential && !terminated && !passivelyReplicated) {type = "potential";}
        else if (!active && !potential && terminated && !passivelyReplicated){type = "terminated";}
        else if (!active && !potential && !terminated && passivelyReplicated){type = "passive";}
        else {type = "unknown";}
        return String.format("%d\t%s\t%s\t%d\t%s\t%d", position, type, leftForkActive, positionLeftFork, rightForkActive, positionRightFork);
    }
    
    
   
    
    /**
     * When a PotentialIS changes state from potential to active,
     * this method updates the linkages among the active states
    * @param lastActivePIS
     */
   void updateLinkagesToAdd (PotentialInitiationSite2 lastActivePIS) {
       leftActivePIS = lastActivePIS;
       rightActivePIS = lastActivePIS.getRightActivePIS();
       lastActivePIS.getRightActivePIS().setLeftActivePIS(this);
       lastActivePIS.setRightActivePIS(this);
   }
    
  /**
   * updates linkages when site terminates
   */
   void updateLinkagesToTerminate() {  //set linkages to reflect that this is no longer active - link left and right activePIS
       leftActivePIS.setRightActivePIS(rightActivePIS);
       rightActivePIS.setLeftActivePIS(leftActivePIS);
   }
   
    
    /**
     * GETTERS - RETRIEVE STATE
     * 
     * @return 
     */
        
   /**
    * get position
    * @return position in molecule
    */ 
   int getPosition () {return position;}
        
   /**
    * ret active PIS to left
    * @return left active PIS
    */ 
    PotentialInitiationSite2 getLeftActivePIS () {return leftActivePIS;}
    
    /**
     * get active PIS to left
     * @return left active PIS
     */
    PotentialInitiationSite2 getRightActivePIS () {return rightActivePIS;}
    
    /**
     * get position of the rightward moving fork
     * @return right fork position
     */
    int getPositionRightFork() {return positionRightFork;}
    
    /**
     * get position of leftward moving fork
     * @return left fork position
     */
    int getPositionLeftFork() {return positionLeftFork;}
        
    /**
     * is site active
     * @return true if active
     */
    boolean IsActive () {return active;} //return true if initiation site is active (has initiated and one or both forks remain active)
        
    /**
     * is site potential
     * @return true if potential
     */
    boolean IsPotential () {return potential;}  //return true if site has not yet initiated synthesis
        
    /**
     * is site terminate
     * @return true if terminated
     */
    boolean IsTerminated() {return terminated;}  //return true if two forks have met on both sides - site no longer active
    
    /**
     * is site passively replicated
     * @return true if passively replicated
     */
    boolean IsPassivelyReplicated() {return passivelyReplicated;}
    
    /**
     * is left fork active
     * @return true if left fork is active
     */
    boolean IsLeftForkActive() {return leftForkActive;}
    
    /**
     * is right fork active
     * @return true if right fork active
     */
    boolean IsRightForkActive() {return rightForkActive;}

    /**
     * Returns index of site in molecule
     * @return index
     */
    int getIndex() {
        return index;
    }
    
    
    
    
    
   
    
    /**
     * FOR TESTING PURPOSES - insert after test for active site in synthesis engine
     * Ensures that active sites are linked up appropriately 
     */
    public enum SiteType {RIGHT_ACTIVE,LEFT_ACTIVE,BOTH_ACTIVE, NEITHER_ACTIVE, EITHER_ACTIVE}
    
    SiteType getSiteType() {
        if (leftForkActive && rightForkActive) {return SiteType.BOTH_ACTIVE;}
        else if (leftForkActive && ! rightForkActive) {return SiteType.LEFT_ACTIVE;}
        else if (rightForkActive && !leftForkActive) {return SiteType.RIGHT_ACTIVE;}
        else if (rightForkActive || leftForkActive) {return SiteType.EITHER_ACTIVE;}
        else return SiteType.NEITHER_ACTIVE;
    }
    
    void linkageCheck(int seqLength) {
        switch (getSiteType()) {
            case BOTH_ACTIVE:
                if ((rightActivePIS.getSiteType() == SiteType.LEFT_ACTIVE || rightActivePIS.getSiteType() == SiteType.BOTH_ACTIVE) &&
                        (leftActivePIS.getSiteType() == SiteType.RIGHT_ACTIVE || leftActivePIS.getSiteType() == SiteType.BOTH_ACTIVE)) {
                    return;
                }
                break;
            case LEFT_ACTIVE:
                if ((leftActivePIS.getSiteType() == SiteType.RIGHT_ACTIVE || leftActivePIS.getSiteType() == SiteType.BOTH_ACTIVE)  &&
                        (rightActivePIS.getSiteType() == SiteType.RIGHT_ACTIVE || rightActivePIS.getPosition() == seqLength)) {
                    return;
                }   break;
            case RIGHT_ACTIVE:
                if ((rightActivePIS.getSiteType() == SiteType.LEFT_ACTIVE || rightActivePIS.getSiteType() == SiteType.BOTH_ACTIVE) &&
                        (leftActivePIS.getSiteType() == SiteType.LEFT_ACTIVE || leftActivePIS.getPosition() == -1)) {
                    return;
                }   break;
            case NEITHER_ACTIVE:
                break;
            case EITHER_ACTIVE:
                break;
            default:  
                break;
        }
        System.out.println("Linkage Error" + "\t" + position);
    }
}

    
    

