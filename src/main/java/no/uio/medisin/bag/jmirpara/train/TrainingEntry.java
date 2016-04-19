/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara.train;

import java.util.ArrayList;
import no.uio.medisin.bag.core.PriMiRNA;

/**
 *
 * @author weibo
 */
public class TrainingEntry {

    private PriMiRNA pri;
    ArrayList<MiEntry> mis;

    public TrainingEntry(){
        pri=new PriMiRNA();
        mis=new ArrayList<MiEntry>();
    }

    public void setPriRNA(PriMiRNA pri){
        this.pri=pri;
    }
    public PriMiRNA getPriRNA(){
        return pri;
    }
    public void addMiEntry(){
        mis.add(new MiEntry());
    }
    public void addMiEntry(int start,int size,int strand,int label){
        mis.add(new MiEntry(start,size,strand,label));
    }
    public int NumOfMi(){
        return mis.size();
    }
    public void setStart(int s){
        mis.get(mis.size()-1).start=s;
    }
    public int getStart(int i){
        return mis.get(i).start;
    }
    public void setSize(int s){
        mis.get(mis.size()-1).size=s;
    }
    public int getSize(int i){
        return mis.get(i).size;
    }
    public void setStrand(int s){
        mis.get(mis.size()-1).strand=s;
    }
    public int getStrand(int i){
        return mis.get(i).strand;
    }
    public void setLabel(int l){
        mis.get(mis.size()-1).label=l;
    }
    public int getLabel(int i){
        return mis.get(i).label;
    }

//    public void setMiId(String id){
//        mis.get(mis.size()-1).id=id;
//    }
//    public String getMiId(int i){
//        return mis.get(i).id;
//    }
//    public void setMiSeq(String seq){
//        mis.get(mis.size()-1).seq=seq;
//    }
//    public String getMiSeq(int i){
//        return mis.get(i).seq;
//    }


    public void resetMis(){
        mis=new ArrayList<MiEntry>();
    }

    private class MiEntry{
        private int start;
        private int size;
        private int strand;
        private int label;

//        private String id;
//        private String seq;
        public MiEntry(){

        }
        public MiEntry(int start,int size,int strand,int label){
            this.start=start;
            this.size=size;
            this.strand=strand;
            this.label=label;
//            this.seq=pri.getSeq().substring(start, start+size);
        }
    }


}
