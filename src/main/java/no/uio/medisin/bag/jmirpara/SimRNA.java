
package no.uio.medisin.bag.jmirpara;

/**
 * simple RNA class,
 * extended by PriMiRNA,PreMiRNA and MiRNA
 * extends SimSeq
 * @author weibo
 */
public class SimRNA extends SimSeq {

    public SimRNA(){
        super();
    }
    public SimRNA(String id,String seq){
        super(id,seq);
    }
    public SimRNA(SimSeq seq){
        this.setId(seq.getId());
        this.setName(seq.getName());
        this.setSeq(seq.getSeq());
        this.setStart(seq.getStart());
        this.setEnd(seq.getEnd());
        this.setLength(seq.getLength());
    }
    private String str="";
    private float energy=0;
    private float GC_content=0;
    private int GU_num=0;
    private int pair_num=0;
    private float A_content=0;
    private float U_content=0;
    private float G_content=0;
    private float C_content=0;

    /**
     * @return the str
     */
    public String getStr() {
        return str;
    }

    /**
     * @param str the str to set
     */
    public void setStr(String str) {
        this.str = str;
    }

    /**
     * @return the energy
     */
    public float getEnergy() {
        return energy;
    }

    /**
     * @param energy the energy to set
     */
    public void setEnergy(float energy) {
        this.energy = energy;
    }

    /**
     * @return the GC_content
     */
    public float getGC_content() {
        return GC_content;
    }

    /**
     * @param GC_content the GC_content to set
     */
    public void setGC_content(float GC_content) {
        this.GC_content = GC_content;
    }

    /**
     * @return the GU_num
     */
    public int getGU_num() {
        return GU_num;
    }

    /**
     * @param GU_num the GU_num to set
     */
    public void setGU_num(int GU_num) {
        this.GU_num = GU_num;
    }

    /**
     * @return the pair_num
     */
    public int getPair_num() {
        return pair_num;
    }

    /**
     * @param pair_num the pair_num to set
     */
    public void setPair_num(int pair_num) {
        this.pair_num = pair_num;
    }

    /**
     * @return the A_content
     */
    public float getA_content() {
        return A_content;
    }

    /**
     * @param A_content the A_content to set
     */
    public void setA_content(float A_content) {
        this.A_content = A_content;
    }

    /**
     * @return the U_content
     */
    public float getU_content() {
        return U_content;
    }

    /**
     * @param U_content the U_content to set
     */
    public void setU_content(float U_content) {
        this.U_content = U_content;
    }

    /**
     * @return the G_content
     */
    public float getG_content() {
        return G_content;
    }

    /**
     * @param G_content the G_content to set
     */
    public void setG_content(float G_content) {
        this.G_content = G_content;
    }

    /**
     * @return the C_content
     */
    public float getC_content() {
        return C_content;
    }

    /**
     * @param C_content the C_content to set
     */
    public void setC_content(float C_content) {
        this.C_content = C_content;
    }


}