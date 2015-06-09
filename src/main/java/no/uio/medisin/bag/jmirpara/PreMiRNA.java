/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author weibo
 */
public class PreMiRNA extends SimRNA{

    public PreMiRNA(){
        super();
    }
    public PreMiRNA(String id,String seq){
        super(id,seq);
    }

    private int maxInternalLoopSize=0;
    private int internalLoop_num=0;
    private int unpairedBase_num=0;
    private double unpairedBase_rate=0;

    private int lowerStemSize=0;
    private int upperStemSize=0;
    private int topStemSize=0;
    private int upperStart=0;
    private int upperEnd=0;
    private int lowerStemInternalLoopSize=0;
    private int topStemInternalLoopSize=0;
    private int lowerStemInternalLoop_num=0;
    private int topStemInternalLoop_num=0;
    private int lowerStemUnpairedBase_num=0;
    private int topStemUnpairedBase_num=0;
    private double lowerStemUnpairedBase_rate=0;
    private double topStemUnpairedBase_rate=0;

    private HashMap featureSet;

    private ArrayList<MatMiRNA> mis=new ArrayList<MatMiRNA>();
    private int index=0;

    public void addProduct(MatMiRNA miRNA){
        mis.add(miRNA);
        index+=1;
    }


    public int SizeOfProduct(){
        return mis.size();
    }

    public MatMiRNA NextMiRNA(){
        index-=1;
        return mis.get(index);
    }

    public void setFeatureSet(){
        featureSet=new HashMap();
        featureSet.put("preRNA_sequence", this.getSeq());
        featureSet.put("preRNA_structure", this.getStr());
        featureSet.put("preRNA_energy", this.getEnergy());
        featureSet.put("preRNA_size", this.getLength());
        featureSet.put("preRNA_GC_content", this.getGC_content() );
        featureSet.put("preRNA_A_content", this.getA_content());
        featureSet.put("preRNA_U_content", this.getU_content());
        featureSet.put("preRNA_G_content", this.getG_content());
        featureSet.put("preRNA_C_content", this.getC_content());
        featureSet.put("preRNA_pair_number", this.getPair_num());
        featureSet.put("preRNA_G-U_number", this.getGU_num());
        featureSet.put("preRNA_unpair_number", this.getUnpairedBase_num());
        featureSet.put("preRNA_unpair_rate", this.getUnpairedBase_rate());
        featureSet.put("preRNA_internalLoop_number", this.getInternalLoop_num());
        featureSet.put("preRNA_internalLoop_size", this.getInternalLoopSize());
        featureSet.put("upperStem_start", this.getUpperStart());
        featureSet.put("upperStem_end", this.getUpperEnd());
        featureSet.put("upperStem_size", this.getUpperStemSize());
        featureSet.put("lowerStem_size", this.getLowerStemSize());
        featureSet.put("lowerStem_unpair_number", this.getLowerStemUnpairedBase_num());
        featureSet.put("lowerStem_unpair_rate", this.getLowerStemUnpairedBase_rate());
        featureSet.put("lowerStem_internalLoop_number", this.getLowerStemInternalLoop_num());
        featureSet.put("lowerStem_internalLoop_size", this.getLowerStemInternalLoopSize());
        featureSet.put("topStem_size", this.getTopStemSize());
        featureSet.put("topStem_unpair_number", this.getTopStemUnpairedBase_num());
        featureSet.put("topStem_unpair_rate", this.getTopStemUnpairedBase_rate());
        featureSet.put("topStem_internalLoop_number", this.getTopStemInternalLoop_num());
        featureSet.put("topStem_internalLoop_size", this.getTopStemInternalLoopSize());
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
     * @return the lowerStemSize
     */
    public int getLowerStemSize() {
        return lowerStemSize;
    }

    /**
     * @param lowerStemSize the lowerStemSize to set
     */
    public void setLowerStemSize(int lowerStemSize) {
        this.lowerStemSize = lowerStemSize;
    }

    /**
     * @return the upperStemSize
     */
    public int getUpperStemSize() {
        return upperStemSize;
    }

    /**
     * @param upperStemSize the upperStemSize to set
     */
    public void setUpperStemSize(int upperStemSize) {
        this.upperStemSize = upperStemSize;
    }

    /**
     * @return the topStemSize
     */
    public int getTopStemSize() {
        return topStemSize;
    }

    /**
     * @param topStemSize the topStemSize to set
     */
    public void setTopStemSize(int topStemSize) {
        this.topStemSize = topStemSize;
    }

    /**
     * @return the lowerStemInternalLoopSize
     */
    public int getLowerStemInternalLoopSize() {
        return lowerStemInternalLoopSize;
    }

    /**
     * @param lowerStemInternalLoopSize the lowerStemInternalLoopSize to set
     */
    public void setLowerStemInternalLoopSize(int lowerStemInternalLoopSize) {
        this.lowerStemInternalLoopSize = lowerStemInternalLoopSize;
    }

    /**
     * @return the topStemInternalLoopSize
     */
    public int getTopStemInternalLoopSize() {
        return topStemInternalLoopSize;
    }

    /**
     * @param topStemInternalLoopSize the topStemInternalLoopSize to set
     */
    public void setTopStemInternalLoopSize(int topStemInternalLoopSize) {
        this.topStemInternalLoopSize = topStemInternalLoopSize;
    }

    /**
     * @return the lowerStemInternalLoop_num
     */
    public int getLowerStemInternalLoop_num() {
        return lowerStemInternalLoop_num;
    }

    /**
     * @param lowerStemInternalLoop_num the lowerStemInternalLoop_num to set
     */
    public void setLowerStemInternalLoop_num(int lowerStemInternalLoop_num) {
        this.lowerStemInternalLoop_num = lowerStemInternalLoop_num;
    }

    /**
     * @return the topStemInternalLoop_num
     */
    public int getTopStemInternalLoop_num() {
        return topStemInternalLoop_num;
    }

    /**
     * @param topStemInternalLoop_num the topStemInternalLoop_num to set
     */
    public void setTopStemInternalLoop_num(int topStemInternalLoop_num) {
        this.topStemInternalLoop_num = topStemInternalLoop_num;
    }

    /**
     * @return the lowerStemUnpairedBase_num
     */
    public int getLowerStemUnpairedBase_num() {
        return lowerStemUnpairedBase_num;
    }

    /**
     * @param lowerStemUnpairedBase_num the lowerStemUnpairedBase_num to set
     */
    public void setLowerStemUnpairedBase_num(int lowerStemUnpairedBase_num) {
        this.lowerStemUnpairedBase_num = lowerStemUnpairedBase_num;
    }

    /**
     * @return the topStemUnpairedBase_num
     */
    public int getTopStemUnpairedBase_num() {
        return topStemUnpairedBase_num;
    }

    /**
     * @param topStemUnpairedBase_num the topStemUnpairedBase_num to set
     */
    public void setTopStemUnpairedBase_num(int topStemUnpairedBase_num) {
        this.topStemUnpairedBase_num = topStemUnpairedBase_num;
    }

    /**
     * @return the lowerStemUnpairedBase_rate
     */
    public double getLowerStemUnpairedBase_rate() {
        return lowerStemUnpairedBase_rate;
    }

    /**
     * @param lowerStemUnpairedBase_rate the lowerStemUnpairedBase_rate to set
     */
    public void setLowerStemUnpairedBase_rate(double lowerStemUnpairedBase_rate) {
        this.lowerStemUnpairedBase_rate = lowerStemUnpairedBase_rate;
    }

    /**
     * @return the topStemUnpairedBase_rate
     */
    public double getTopStemUnpairedBase_rate() {
        return topStemUnpairedBase_rate;
    }

    /**
     * @param topStemUnpairedBase_rate the topStemUnpairedBase_rate to set
     */
    public void setTopStemUnpairedBase_rate(double topStemUnpairedBase_rate) {
        this.topStemUnpairedBase_rate = topStemUnpairedBase_rate;
    }

    public void setUpperStart(int upperStart) {
        this.upperStart = upperStart;
    }

    /**
     * @return the upperStart
     */
    public int getUpperStart() {
        return upperStart;
    }

    /**
     * @return the upperEnd
     */
    public int getUpperEnd() {
        return upperEnd;
    }

    /**
     * @param upperEnd the upperEnd to set
     */
    public void setUpperEnd(int upperEnd) {
        this.upperEnd = upperEnd;
    }



}

