/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author weibo
 */
public class foldMirbase {
    public static void main(String[] argv) throws IOException{
        PathSet.setLibDir("/home/weibo/NetBeansProjects/JmiPara2/lib");
        ReadFastaFile hps=new ReadFastaFile(new File("/home/weibo/Desktop/jmi/mirbase_data/v13/hairpin.fa"));
        BufferedWriter br=new BufferedWriter(new FileWriter("/home/weibo/Desktop/sl.rna"));
        ArrayList<SimpleSeq> seqs=hps.getSeqs();
        for(SimpleSeq seq:seqs){
            SimpleRNASequence rna=new SimpleRNASequence(seq);
            ParaToolKit.foldRNA(rna);
            br.write(rna.getId()+"\n"+rna.getSeq()+"\n"+rna.getStr()+"\n"+rna.getEnergy()+"\n");
        }
        br.close();
    }
}
