/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.util.HashMap;

/**
 *
 * @author weibo
 */
public class MatMiRNA extends SimRNA{

    public MatMiRNA(){
        super();
    }
    public MatMiRNA(String id,String seq){
        super(id,seq);
    }

    private int maxInternalLoopSize=0;
    private int internalLoop_num=0;
    private int unpairedBase_num=0;
    private double unpairedBase_rate=0;

    private int strand;
    private char firstBase;
    private double stability=0;
    private char dangleBaseOne;
    private char dangleBaseTwo;

    private int miStart;
    private int miEnd;

    private String note="";

    private HashMap featureSet;

    public void setFeatureSet(){
        featureSet=new HashMap();
        featureSet.put("miRNA_id", this.getId());
        featureSet.put("miRNA_sequence", this.getSeq());
        featureSet.put("miRNA_structure", this.getStr());
        featureSet.put("miRNA_energy", this.getEnergy());
        featureSet.put("miRNA_size", this.getLength());
        featureSet.put("miRNA_GC_content", this.getGC_content() );
        featureSet.put("miRNA_A_content", this.getA_content());
        featureSet.put("miRNA_U_content", this.getU_content());
        featureSet.put("miRNA_G_content", this.getG_content());
        featureSet.put("miRNA_C_content", this.getC_content());
        featureSet.put("miRNA_pair_number", this.getPair_num());
        featureSet.put("miRNA_G-U_number", this.getGU_num());
        featureSet.put("miRNA_unpair_number", this.getUnpairedBase_num());
        featureSet.put("miRNA_unpair_rate", this.getUnpairedBase_rate());
        featureSet.put("miRNA_internalLoop_number", this.getInternalLoop_num());
        featureSet.put("miRNA_internalLoop_size", this.getInternalLoopSize());
        featureSet.put("miRNA_start", this.getMiStart());
        featureSet.put("miRNA_end", this.getMiEnd());
        featureSet.put("miRNA_stability", this.getStability());
        featureSet.put("miRNA_firstBase", this.getFirstBase());
        featureSet.put("overhang_base1", this.getDangleBaseOne());
        featureSet.put("overhang_base2", this.getDangleBaseTwo());
        featureSet.put("strand", this.getStrand());

        featureSet.put("miStart", this.getStart());
        featureSet.put("miEnd", this.getEnd());
    }

    public HashMap getFeatureSet(){
        return featureSet;
    }

    /**
     * @return the maxInternalLoopSize
     */
    public int getInternalLoopSize() {
        return maxInternalLoopSize;
    }

    /**
     * @param maxInternalLoopSize the maxInternalLoopSize to set
     */
    public void setInternalLoopSize(int maxInternalLoopSize) {
        this.maxInternalLoopSize = maxInternalLoopSize;
    }

    /**
     * @return the internalLoop_num
     */
    public int getInternalLoop_num() {
        return internalLoop_num;
    }

    /**
     * @param internalLoop_num the internalLoop_num to set
     */
    public void setInternalLoop_num(int internalLoop_num) {
        this.internalLoop_num = internalLoop_num;
    }

    /**
     * @return the unpairedBase_num
     */
    public int getUnpairedBase_num() {
        return unpairedBase_num;
    }

    /**
     * @param unpairedBase_num the unpairedBase_num to set
     */
    public void setUnpairedBase_num(int unpairedBase_num) {
        this.unpairedBase_num = unpairedBase_num;
    }

    /**
     * @return the unpairedBase_rate
     */
    public double getUnpairedBase_rate() {
        return unpairedBase_rate;
    }

    /**
     * @param unpairedBase_rate the unpairedBase_rate to set
     */
    public void setUnpairedBase_rate(double unpairedBase_rate) {
        this.unpairedBase_rate = unpairedBase_rate;
    }

    /**
     * @return the strand
     */
    public int getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(int strand) {
        this.strand = strand;
    }

    /**
     * @return the firstBase
     */
    public char getFirstBase() {
        return firstBase;
    }

    /**
     * @param firstBase the firstBase to set
     */
    public void setFirstBase(char firstBase) {
        this.firstBase = firstBase;
    }

    /**
     * @return the stability
     */
    public double getStability() {
        return stability;
    }

    /**
     * @param stability the stability to set
     */
    public void setStability(double stability) {
        this.stability = stability;
    }

    /**
     * @return the dangleBaseOne
     */
    public char getDangleBaseOne() {
        return dangleBaseOne;
    }

    /**
     * @param dangleBaseOne the dangleBaseOne to set
     */
    public void setDangleBaseOne(char dangleBaseOne) {
        this.dangleBaseOne = dangleBaseOne;
    }

    /**
     * @return the dangleBaseTwo
     */
    public char getDangleBaseTwo() {
        return dangleBaseTwo;
    }

    /**
     * @param dangleBaseTwo the dangleBaseTwo to set
     */
    public void setDangleBaseTwo(char dangleBaseTwo) {
        this.dangleBaseTwo = dangleBaseTwo;
    }

    /**
     * @return the miEnd
     */
    public int getMiEnd() {
        return miEnd;
    }

    /**
     * @param miEnd the miEnd to set
     */
    public void setMiEnd(int miEnd) {
        this.miEnd = miEnd;
    }

    /**
     * @return the miStart
     */
    public int getMiStart() {
        return miStart;
    }

    /**
     * @param miStart the miStart to set
     */
    public void setMiStart(int miStart) {
        this.miStart = miStart;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }

}
