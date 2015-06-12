/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author weibo
 */
public class InputData {

    private ArrayList<SimpleRNASequence> seqs;
    private ArrayList<SimpleRNASequence> segs;
    
    public ArrayList load(String fname) throws IOException{
        readFastaFile(fname);
        return seqs;
    }
   

    /**
    * Read in an alignment in FASTA format.  Stores results in a SimpleSeq List
    *
    * @throws java.io.IOException
    */
    private void readFastaFile(String fname) throws java.io.IOException{

        seqs=new ArrayList<SimpleRNASequence>();

        BufferedReader br = new BufferedReader(new FileReader(fname));
        String line = "";
        String title = "";
        String seq = "";

        title = br.readLine().replaceAll("[\\/\\s\\|,]+", "_");
        if(title.substring(0, 1).equals(">") == false){
            throw new java.io.IOException("bad FASTA format on first line\n");
        }
        title = title.substring(1, title.length()).replaceAll("[\\/\\s\\|,]+", "_");

        while((line = br.readLine()) != null){
            if(line.startsWith(">")==true){
                line = line.substring(1, line.length());
                addSequence(title, seq);
                title = line.replaceAll("[\\/\\s\\|,]+", "_");
                seq = "";
            }
            else{
                line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                seq = seq.concat(line);
            }
       }
       addSequence(title,seq);
       br.close();
    }

    private void addSequence(String title, String seq) {
        SimpleRNASequence entry=new SimpleRNASequence(title,seq);
        entry.setName(title);
        seqs.add(entry);
        /**********debug***********************/
        System.out.println("Loads sequence "+title);
        /**********debug***********************/
    }

    /**
     * cut a sequence into a series of fragments
     * with specified window and step size
     */
    public ArrayList<SimpleRNASequence> slide(int window, int step, SimpleRNASequence seq){
        segs=new ArrayList<SimpleRNASequence>();
        int length=seq.getLength();
        int n=(length-window)/step;
        if(n<0) n=0;
        if(length-n*step+window>19) n+=1;
        int start=0,end=0;
        for(int i=0;i<=n;i++){
            start=i*step;
            if(start>=length) break;
            end=start+window;
            if(end>length) end=length;
            String id=seq.getName()+"_"+(start+1)+"-"+end;
            String subseq=seq.getSeq().substring(start,end);

            SimpleRNASequence frag=new SimpleRNASequence(id,subseq);
            frag.setStart(start+1);// count from 1
            frag.setEnd(end); //count from 1
            frag.setName(seq.getName());
            segs.add(frag);
        }
        /**********debug***********************/
        System.out.println(segs.size()+" fragments are generated from "+seq.getId());
        /**********debug***********************/
        return segs;
    }
}
