
package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author weibo
 */
public class PriMiRNA extends SimRNA{

    public PriMiRNA() {
        super();
    }

    public PriMiRNA(String id, String seq) {
        super(id, seq);
    }

    private int basalSegEnd = 0;//count from 1
    private int basalSegSize = 0;
    private String terminalLoopSeq;
    private int terminalLoopStart = 0;//count from 1
    private int terminalLoopEnd = 0;//count from 1
    private int terminalLoopSize = 0;
    private int internalLoopSize = 0;
    private int internalLoop_num = 0;
    private int unpairedBase_num = 0;
    private double unpairedBase_rate = 0;
    private String seq5; //5' seq of pri
    private String str5; //5' structure of pri
    private String seq3; //3' seq of pri
    private String str3; //3' structure of pri
    private String midBase; //mid-base of terminal loop
    private String priLine1;
    private String priLine2;
    private String priLine3;
    private String priLine4;
    private String priPlot;
    private int basalBaseNum = 0;
    private int[] strIndex;
    private HashMap featureSet;


    private ArrayList<PreMiRNA> pres=new ArrayList<PreMiRNA>();
    private int index=0;

    public void addProduct(PreMiRNA preRNA){
        pres.add(preRNA);
        index+=1;
    }

    public int sizeOfProduct(){
        return pres.size();
    }

    public PreMiRNA nextPreRNA(){
        index-=1;
        return pres.get(index);
    }

    public void setFeatureSet(){
        featureSet=new HashMap();
        featureSet.put("priRNA_id", this.getId());
        featureSet.put("priRNA_start", this.getStart());
        featureSet.put("priRNA_end", this.getEnd());
        featureSet.put("priRNA_sequence", this.getSeq());
        featureSet.put("priRNA_structure", this.getStr());
        featureSet.put("priRNA_energy", this.getEnergy());
        featureSet.put("priRNA_size", this.getLength());
        featureSet.put("priRNA_plot", MiGenesis.setStrPlot(this));
        featureSet.put("priRNA_GC_content", this.getGC_content() );
        featureSet.put("priRNA_A_content", this.getA_content());
        featureSet.put("priRNA_U_content", this.getU_content());
        featureSet.put("priRNA_G_content", this.getG_content());
        featureSet.put("priRNA_C_content", this.getC_content());
        featureSet.put("priRNA_pair_number", this.getPair_num());
        featureSet.put("priRNA_G-U_number", this.getGU_num());
        featureSet.put("priRNA_unpair_number", this.getUnpairedBase_num());
        featureSet.put("priRNA_unpair_rate", this.getUnpairedBase_rate());
        featureSet.put("priRNA_internalLoop_number", this.getInternalLoop_num());
        featureSet.put("priRNA_internalLoop_size", this.getInternalLoopSize());
        featureSet.put("basalSegment_size", this.getBasalSegSize());
        featureSet.put("basalSegment_end", this.getBasalSegEnd());
        featureSet.put("terminalLoop_size", this.getTerminalLoopSize());
    }

    public HashMap getFeatureSet(){
        return featureSet;
    }

    /**
     * @return the basalSegEnd
     */
    public int getBasalSegEnd() {
        return basalSegEnd;
    }

    /**
     * @param basalSegEnd the basalSegEnd to set
     */
    public void setBasalSegEnd(int basalSegEnd) {
        this.basalSegEnd = basalSegEnd;
    }

    /**
     * @return the basalSegSize
     */
    public int getBasalSegSize() {
        return basalSegSize;
    }

    /**
     * @param basalSegSize the basalSegSize to set
     */
    public void setBasalSegSize(int basalSegSize) {
        this.basalSegSize = basalSegSize;
    }

    /**
     * @return the terminalLoopSeq
     */
    public String getTerminalLoopSeq() {
        return terminalLoopSeq;
    }

    /**
     * @param terminalLoopSeq the terminalLoopSeq to set
     */
    public void setTerminalLoopSeq(String terminalLoopSeq) {
        this.terminalLoopSeq = terminalLoopSeq;
    }

    /**
     * @return the terminalLoopStart
     */
    public int getTerminalLoopStart() {
        return terminalLoopStart;
    }

    /**
     * @param terminalLoopStart the terminalLoopStart to set
     */
    public void setTerminalLoopStart(int terminalLoopStart) {
        this.terminalLoopStart = terminalLoopStart;
    }

    /**
     * @return the terminalLoopEnd
     */
    public int getTerminalLoopEnd() {
        return terminalLoopEnd;
    }

    /**
     * @param terminalLoopEnd the terminalLoopEnd to set
     */
    public void setTerminalLoopEnd(int terminalLoopEnd) {
        this.terminalLoopEnd = terminalLoopEnd;
    }

    /**
     * @return the terminalLoopSize
     */
    public int getTerminalLoopSize() {
        return terminalLoopSize;
    }

    /**
     * @param terminalLoopSize the terminalLoopSize to set
     */
    public void setTerminalLoopSize(int terminalLoopSize) {
        this.terminalLoopSize = terminalLoopSize;
    }

    /**
     * @return the internalLoopSize
     */
    public int getInternalLoopSize() {
        return internalLoopSize;
    }

    /**
     * @param internalLoopSize the internalLoopSize to set
     */
    public void setInternalLoopSize(int internalLoopSize) {
        this.internalLoopSize = internalLoopSize;
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
     * @return the seq5
     */
    public String getSeq5() {
        return seq5;
    }

    /**
     * @param seq5 the seq5 to set
     */
    public void setSeq5(String seq5) {
        this.seq5 = seq5;
    }

    /**
     * @return the str5
     */
    public String getStr5() {
        return str5;
    }

    /**
     * @param str5 the str5 to set
     */
    public void setStr5(String str5) {
        this.str5 = str5;
    }

    /**
     * @return the seq3
     */
    public String getSeq3() {
        return seq3;
    }

    /**
     * @param seq3 the seq3 to set
     */
    public void setSeq3(String seq3) {
        this.seq3 = seq3;
    }

    /**
     * @return the str3
     */
    public String getStr3() {
        return str3;
    }

    /**
     * @param str3 the str3 to set
     */
    public void setStr3(String str3) {
        this.str3 = str3;
    }

    /**
     * @return the midBase
     */
    public String getMidBase() {
        return midBase;
    }

    /**
     * @param midBase the midBase to set
     */
    public void setMidBase(String midBase) {
        this.midBase = midBase;
    }

    /**
     * @return the priLine1
     */
    public String getPriLine1() {
        return priLine1;
    }

    /**
     * @param priLine1 the priLine1 to set
     */
    public void setPriLine1(String priLine1) {
        this.priLine1 = priLine1;
    }

    /**
     * @return the priLine2
     */
    public String getPriLine2() {
        return priLine2;
    }

    /**
     * @param priLine2 the priLine2 to set
     */
    public void setPriLine2(String priLine2) {
        this.priLine2 = priLine2;
    }

    /**
     * @return the priLine3
     */
    public String getPriLine3() {
        return priLine3;
    }

    /**
     * @param priLine3 the priLine3 to set
     */
    public void setPriLine3(String priLine3) {
        this.priLine3 = priLine3;
    }

    /**
     * @return the priLine4
     */
    public String getPriLine4() {
        return priLine4;
    }

    /**
     * @param priLine4 the priLine4 to set
     */
    public void setPriLine4(String priLine4) {
        this.priLine4 = priLine4;
    }

    /**
     * @return the priPlot
     */
    public String getPriPlot() {
        return priPlot;
    }

    /**
     * @param priPlot the priPlot to set
     */
    public void setPriPlot(String priPlot) {
        this.priPlot = priPlot;
    }

    /**
     * @return the basalBaseNum
     */
    public int getBasalBaseNum() {
        return basalBaseNum;
    }

    /**
     * @param basalBaseNum the basalBaseNum to set
     */
    public void setBasalBaseNum(int basalBaseNum) {
        this.basalBaseNum = basalBaseNum;
    }

    

    

    /**
     * @return the strIndex
     */
    public int[] getStrIndex() {
        return strIndex;
    }

    /**
     * @param strIndex the strIndex to set
     */
    public void setStrIndex(int[] strIndex) {
        this.strIndex = strIndex;
    }

    
}
