package no.uio.medisin.bag.jmirpara.train;

import java.util.ArrayList;

/**
 *
 * @author sr
 */
public class MiRBaseEntry {

  private String priAccession;
  private String priID;
  private String priSequence;
  private String orgnism;
  private String taxonomy;

  private ArrayList<MatEntry> mat=new ArrayList<MatEntry>();

  public int getMatStart(int i){
      return mat.get(i).matStart;
  }
  public int getMatEnd(int i){
      return mat.get(i).matEnd;
  }
  public String getMatAcc(int i){
      return mat.get(i).matAccession;
  }
  public String getMatID(int i){
      return mat.get(i).matID;
  }
  public String getEvidence(int i) {
        return mat.get(i).evidence;
  }

  public int getMatNum(){
      return mat.size();
  }

  public void addMat(){
      mat.add(new MatEntry());
  }

  public void setMatStart(int s){
      mat.get(mat.size()-1).matStart=s;
  }
  public void setMatEnd(int e){
      mat.get(mat.size()-1).matEnd=e;
  }
  public void setMatAcc(String a){
      mat.get(mat.size()-1).matAccession=a;
  }
  public void setMatID(String i){
      mat.get(mat.size()-1).matID=i;
  }
  public void setEvidence(String e){
      mat.get(mat.size()-1).evidence=e;
  }



  /**
     * @return the priAccession
     */
    public String getPriAccession() {
        return priAccession;

    }

    /**
     * @param priAccession the priAccession to set
     */
    public void setPriAccession(String priAccession) {
        this.priAccession = priAccession;
    }

    /**
     * @return the priID
     */
    public String getPriID() {
        return priID;
    }

    /**
     * @param priID the priID to set
     */
    public void setPriID(String priID) {
        this.priID = priID;
    }

    /**
     * @return the priSequence
     */
    public String getPriSequence() {
        return priSequence;
    }

    /**
     * @param priSequence the priSequence to set
     */
    public void setPriSequence(String priSequence) {
        this.priSequence = priSequence;
    }


    /**
     * @return the orgnism
     */
    public String getOrgnism() {
        return orgnism;
    }

    /**
     * @param orgnism the orgnism to set
     */
    public void setOrgnism(String orgnism) {
        this.orgnism = orgnism;
    }

    /**
     * @return the taxology
     */
    public String getTaxonomy() {
        return taxonomy;
    }

    /**
     * @param taxology the taxology to set
     */
    public void setTaxonomy(String taxonomy) {
        this.taxonomy = taxonomy;
    }


    private class MatEntry{
      private int matStart=0;
      private int matEnd=0;
      private String matAccession="";
      private String matID="";
      private String evidence="";
      private int matNum=0;
    }


}
