/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

/**
 * a simple sequence object that contains some basic functionality
 * and stores basic sequence features
 * 
 * @author weibo
 */
public class SimpleSeq {
    private String id="";
    private String seq="";
    private int length=0;
    private int start=0;
    private int end=0;
    private String name="";


    /**
     * Empty Class Constructor
     */
    public SimpleSeq(){

    }
    
    
    /**
     * Constructor 
     * 
     * @param id
     * @param seq 
     */
    public SimpleSeq(String id,String seq){
        
        this.id=id;
        this.seq=seq;
        this.length=seq.length();
        
    }

    /**
     * Reverse complement RNA/DNA sequence
     * if sequence contains 'U' or 'u' it is assumed the sequence is RNA
     * 
     * @param seqIn
     * @return 
     */
    public static String Complement(String seqIn){
        
        StringBuilder Complement = new StringBuilder();
        char [] strReversed = new StringBuilder(seqIn).reverse().toString().toCharArray();

        for (char nt: strReversed) {
            switch (nt){
                case 'a':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("t");
                    else
                        Complement.append("u");
                    break;
                    
                case 'A':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("T");
                    else
                        Complement.append("U");
                    break;
                    
                case 'c':
                    Complement.append("g");
                    break;
                    
                case 'C':
                    Complement.append("G");
                    break;
                    
                case 'g':
                    Complement.append("c");
                    break;
                    
                case 'G':
                    Complement.append("C");
                    break;
                    
                case 't':
                    Complement.append("a");
                    break;
                    
                case 'T':
                    Complement.append("A");
                    break;
                    
                case 'u':
                    Complement.append("a");
                    break;
                    
                case 'U':
                    Complement.append("A");
                    break;
                    
                default:
                    Complement.append("N");
                    break;
                    
            }

        }
        return Complement.toString();
    }



    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id 
     */
    public void setID(String id) {
        this.id = id;
    }

    /**
     * @return the seq
     */
    public String getSeq() {
        return seq;
    }

    /**
     * @param seq 
     */
    public void setSeq(String seq) {
        this.seq = seq;
    }

    /**
     * @return length
     */
    public int getLength() {
        return length;
    }

    /**
     * @param length 
     */
    public void setLength(int length) {
        this.length = length;
    }

    /**
     * @return start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start
     */
    public void setAbsStartInQuerySeq(int start) {
        this.start = start;
    }

    /**
     * @return end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end
     */
    public void setAbsEndInQuerySeq(int end) {
        this.end = end;
    }

    /**
     * @return name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name 
     */
    public void setName(String name) {
        this.name = name;
    }

    


}