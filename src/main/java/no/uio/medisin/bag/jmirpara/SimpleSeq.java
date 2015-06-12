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
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id 
     */
    public void setId(String id) {
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
    public void setStart(int start) {
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
    public void setEnd(int end) {
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