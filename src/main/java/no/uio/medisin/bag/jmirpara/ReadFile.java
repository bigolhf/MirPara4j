/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


/**
 * load the input file
 * @author weibo
 */
public class ReadFile {

    private File fname;
    private String format;
    private BufferedReader br;
    private ArrayList<SimSeq> seqs;
    private double progress=0;
    private String title;
    private StringBuilder s;

    public ReadFile(File fname) throws FileNotFoundException, IOException{
        this.fname=fname;
        br = new BufferedReader(new FileReader(fname));
        title = br.readLine();
        if(title.substring(0, 1).equals(">") == false){
            throw new java.io.IOException("bad FASTA format on first line\n");
        }
        title=title.substring(1, title.length());
    }
//    public ReadFile(File fname) throws IOException{
//        this.fname=fname;
//        readFastaFile();
//    }
    public ReadFile(String fname){
        this.fname=new File(fname);
    }
    public ReadFile(File fname,String format){
        this.fname=fname;
        this.format=format;
    }
    
    public boolean hasSeq(){
        if(br==null) return false;
        else return true;
    }
    
    public SimSeq getOneSeq() throws IOException{
        SimSeq seq=null;
        String line = "";
        StringBuilder s = new StringBuilder();

        
        line=br.readLine();
        while(line!=null){
            if(line.startsWith(">")==true){
                
                seq=new SimSeq(title,s.toString());
                seq.setName(title);
                
                title=line.substring(1, line.length());
                s=null;
                break;
            }
            else{
                line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                s.append(line);
                line=br.readLine();
            }
        }
        if(line==null){
            br.close();
            br=null;
        }
        if(s!=null){
            seq=new SimSeq(title,s.toString());
            seq.setName(title);
        }
        System.out.println("Loaded sequence "+title);
        return seq;
    }
    

    /**
    * Read in an alignment in FASTA format.  Stores results in a SimpleSeq List
    *
    * @throws java.io.IOException
    */
    public void readFastaFile() throws java.io.IOException{
        System.out.print("Start loading data "+fname.getName()+" ...");
        long totalSize=fname.length();
        long readedSize=0;
        
        seqs=new ArrayList<SimSeq>();

        br = new BufferedReader(new FileReader(fname));
        String line = "";
        String title = "";
        StringBuilder seq = new StringBuilder();

        title = br.readLine();
        readedSize+=title.length()+1;
        if(title.substring(0, 1).equals(">") == false){
            throw new java.io.IOException("bad FASTA format on first line\n");
        }
        title = title.substring(1, title.length());

        while((line = br.readLine()) != null){
            readedSize+=line.length()+1;
            if(line.startsWith(">")==true){
                line = line.substring(1, line.length());
                addSequence(title, seq.toString());
                title = line;
                seq = new StringBuilder();
            }
            else{
                line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                seq.append(line);
            }
            progress=readedSize*1.0/totalSize;
            if(progress>1) progress=1;
            System.out.print(Output.decimal(progress*100)+"%"+Output.backspace(Output.decimal(progress*100)+"%"));
       }
       System.out.println();
       addSequence(title,seq.toString());
       br.close();
    }

    private void addSequence(String title, String seq) {
        SimSeq entry=new SimSeq(title,seq);
        entry.setName(title);
        getSeqs().add(entry);
        /**********debug***********************/
        System.out.println("Load sequence "+title);
        /**********debug***********************/
    }

    /**
     * store the sequences from file as SimSeq
     * @return ArrayList<SimSeq>: the seqs
     */
    public ArrayList<SimSeq> getSeqs() {
        return seqs;
    }


  }
