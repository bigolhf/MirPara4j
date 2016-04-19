/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import no.uio.medisin.bag.core.SimpleRNASequence;
import no.uio.medisin.bag.core.SimpleSequenceSet;
import no.uio.medisin.bag.core.SimpleSeq;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author weibo
 */
@Deprecated
public class foldMirbase {
    public static void main(String[] argv) throws IOException{
        PathSet.setLibDir("/home/weibo/NetBeansProjects/JmiPara2/lib");
        SimpleSequenceSet hps=new SimpleSequenceSet(new File("/home/weibo/Desktop/jmi/mirbase_data/v13/hairpin.fa"));
        BufferedWriter br=new BufferedWriter(new FileWriter("/home/weibo/Desktop/sl.rna"));
        ArrayList<SimpleSeq> seqs=hps.getSeqs();
        for(SimpleSeq seq:seqs){
            SimpleRNASequence rna=new SimpleRNASequence(seq);
            ParaToolKit.foldRNA(rna);
            br.write(rna.getId()+"\n"+rna.getSeq()+"\n"+rna.getStructureStr()+"\n"+rna.getEnergy()+"\n");
        }
        br.close();
    }
}
