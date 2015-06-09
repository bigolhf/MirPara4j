/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.BufferedWriter;
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
        ReadFile hps=new ReadFile("/home/weibo/Desktop/jmi/mirbase_data/v13/hairpin.fa");
        BufferedWriter br=new BufferedWriter(new FileWriter("/home/weibo/Desktop/sl.rna"));
        ArrayList<SimSeq> seqs=hps.getSeqs();
        for(SimSeq seq:seqs){
            SimRNA rna=new SimRNA(seq);
            ParaToolKit.foldRNA(rna);
            br.write(rna.getId()+"\n"+rna.getSeq()+"\n"+rna.getStr()+"\n"+rna.getEnergy()+"\n");
        }
        br.close();
    }
}
