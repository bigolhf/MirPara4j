
package no.uio.medisin.bag.jmirpara;

/**
 * simple RNA class
 * stores basic information about an RNASequence
 * 
 * extended by PriMiRNA, PreMiRNA and MiRNA classes
 * extends SimSeq
 * @author weibo
 */
public class SimpleRNASequence extends SimpleSeq {

    
    
    /**
     * Basic Class Constructor
     * 
     */
    public SimpleRNASequence(){
        super();
    }
    
    
    
    
    /**
     * Class Constructor from Sequence specified as Strings
     * 
     * @param id
     * @param seq 
     */
    public SimpleRNASequence(String id, String seq){
        super(id,seq);
    }
    
    
    
    
    /**
     * Class Constructor from SimpleSeq
     * 
     * @param seq 
     */
    public SimpleRNASequence(SimpleSeq seq){
        this.setID(seq.getId());
        this.setName(seq.getName());
        this.setSeq(seq.getSeq());
        this.setAbsStartInQuerySeq(seq.getStart());
        this.setAbsEndInQuerySeq(seq.getEnd());
        this.setLength(seq.getLength());
    }
    
    
    
    private String      structureString="";         // in Vienna Bracket Notation
    private float       energy=0;
    private float       GC_content=0;
    private int         GU_num=0;
    private int         pair_num=0;
    private float       A_content=0;
    private float       U_content=0;
    private float       G_content=0;
    private float       C_content=0;

    
    @Override
    public String toString(){
        String str = "ID:" + this.getId() + "\t"
                    + "Name::" + this.getName() + "\t"
                    + "Start:" + this.getStart() + "\t"
                    + "End:" + this.getEnd() + "\t"
                    + "Len:" + this.getLength() + "\t"
          + "\n";
        return str;
    }
    /**
     * @return str
     */
    public String getStructureStr() {
        return structureString;
    }

    /**
     * @param str 
     */
    public void setStructureString(String str) {
        this.structureString = str;
    }

    /**
     * @return energy
     */
    public float getEnergy() {
        return energy;
    }

    /**
     * @param energy 
     */
    public void setEnergy(float energy) {
        this.energy = energy;
    }

    /**
     * @return GC_content
     */
    public float getGC_content() {
        return GC_content;
    }

    /**
     * @param GC_content 
     */
    public void setGC_content(float GC_content) {
        this.GC_content = GC_content;
    }

    /**
     * @return GU_num
     */
    public int getGU_num() {
        return GU_num;
    }

    /**
     * @param GU_num 
     */
    public void setNumberOfGUPairs(int GU_num) {
        this.GU_num = GU_num;
    }

    /**
     * @return pair_num
     */
    public int getNumberOfPairs() {
        return pair_num;
    }

    /**
     * @param pair_num 
     */
    public void setNumberOfPairedBases(int pair_num) {
        this.pair_num = pair_num;
    }

    /**
     * @return A_content
     */
    public float getA_content() {
        return A_content;
    }

    /**
     * @param A_content 
     */
    public void setA_content(float A_content) {
        this.A_content = A_content;
    }

    /**
     * @return U_content
     */
    public float getU_content() {
        return U_content;
    }

    /**
     * @param U_content 
     */
    public void setU_content(float U_content) {
        this.U_content = U_content;
    }

    /**
     * @return G_content
     */
    public float getG_content() {
        return G_content;
    }

    /**
     * @param G_content
     */
    public void setG_content(float G_content) {
        this.G_content = G_content;
    }

    /**
     * @return C_content
     */
    public float getC_content() {
        return C_content;
    }

    /**
     * @param C_content 
     */
    public void setC_content(float C_content) {
        this.C_content = C_content;
    }


}