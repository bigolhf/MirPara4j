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
public class ReadFastaFile {

    private File fastaFilename;
    private String format;
    private BufferedReader brFA;
    private ArrayList<SimpleSeq> seqs;
    private double progress=0;
    private String title;

    
    
    /**
     * read query sequences in FASTA format
     * 
     * @param fname : query sequence filename
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public ReadFastaFile(File fname) throws FileNotFoundException, IOException{
        this.fastaFilename=fname;
        brFA = new BufferedReader(new FileReader(fname));
        title = brFA.readLine();
        if(title.substring(0, 1).equals(">") == false){
            throw new java.io.IOException("bad FASTA format on first line\n");
        }
        title=title.substring(1, title.length());
    }
    
    
    
    
    /**
     * check whether there are any more entries in the FASTA file
     * @return 
     */
    public boolean hasSeq(){
        
        if(brFA==null) return false;
        return true;
        
    }
    
    public SimpleSeq getOneSeq() throws IOException{
        SimpleSeq seq=null;
        String line = "";
        StringBuilder s = new StringBuilder();

        
        line=brFA.readLine();
        while(line!=null){
            if(line.startsWith(">")==true){
                
                seq=new SimpleSeq(title,s.toString());
                seq.setName(title);
                
                title=line.substring(1, line.length());
                s=null;
                break;
            }
            else{
                line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                s.append(line);
                line=brFA.readLine();
            }
        }
        if(line==null){
            brFA.close();
            brFA=null;
        }
        if(s!=null){
            seq=new SimpleSeq(title,s.toString());
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
        System.out.print("Start loading data "+fastaFilename.getName()+" ...");
        long totalSize=fastaFilename.length();
        long readedSize=0;
        
        seqs=new ArrayList<SimpleSeq>();

        brFA = new BufferedReader(new FileReader(fastaFilename));
        String line = "";
        String title = "";
        StringBuilder seq = new StringBuilder();

        title = brFA.readLine();
        readedSize+=title.length()+1;
        if(title.substring(0, 1).equals(">") == false){
            throw new java.io.IOException("bad FASTA format on first line\n");
        }
        title = title.substring(1, title.length());

        while((line = brFA.readLine()) != null){
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
       brFA.close();
    }

    private void addSequence(String title, String seq) {
        SimpleSeq entry=new SimpleSeq(title,seq);
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
    public ArrayList<SimpleSeq> getSeqs() {
        return seqs;
    }


  }
