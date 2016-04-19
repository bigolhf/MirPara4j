/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

/**
 *
 * @author simonray
 */
public class miRParaPredParams {
    
    private int                     window      = 500;
    private int                     step        = 250;
    private int                     start       = 1;
    private int                     distance    = 60;
    private double                  cutoff      = 0.8;
    private String                  model       = "overall";
    private int                     level       = 1;
    private int                     minMiRNAlen = 1;
    private int                     maxMiRNAlen = 40;

    
    
    
    
    /**
     * @return the window
     */
    public int getWindow() {
        return window;
    }

    /**
     * @param window the window to set
     */
    public void setWindow(int window) {
        this.window = window;
    }

    /**
     * @return the step
     */
    public int getStep() {
        return step;
    }

    /**
     * @param step the step to set
     */
    public void setStep(int step) {
        this.step = step;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the distance
     */
    public int getDistance() {
        return distance;
    }

    /**
     * @param distance the distance to set
     */
    public void setMinHairpinLength(int distance) {
        this.distance = distance;
    }

    /**
     * @return the cutoff
     */
    public double getCutoff() {
        return cutoff;
    }

    /**
     * @param cutoff the cutoff to set
     */
    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
    }

    /**
     * @return the model
     */
    public String getModel() {
        return model;
    }

    /**
     * @param model the model to set
     */
    public void setModel(String model) {
        this.model = model;
    }

    /**
     * @return the level
     */
    public int getLevel() {
        return level;
    }

    /**
     * @param level the level to set
     */
    public void setLevel(int level) {
        this.level = level;
    }

    /**
     * @return the minMiRNAlen
     */
    public int getMinMiRNAlen() {
        return minMiRNAlen;
    }

    /**
     * @param minMiRNAlen the minMiRNAlen to set
     */
    public void setMinMiRNAlen(int minMiRNAlen) {
        this.minMiRNAlen = minMiRNAlen;
    }

    /**
     * @return the maxMiRNAlen
     */
    public int getMaxMiRNAlen() {
        return maxMiRNAlen;
    }

    /**
     * @param maxMiRNAlen the maxMiRNAlen to set
     */
    public void setMaxMiRNAlen(int maxMiRNAlen) {
        this.maxMiRNAlen = maxMiRNAlen;
    }
    
}
