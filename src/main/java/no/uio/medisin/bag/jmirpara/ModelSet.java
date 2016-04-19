/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

/**
 * 
 * @author weibo
 */
public class ModelSet {

    public static String[] model(String m){
        if(m.equals("overall"))
            return overall;
        else if(m.equals("animal"))
            return animal;
        else if(m.equals("virus"))
            return virus;
        else if(m.equals("plant"))
            return plant;
        else if(m.equals("range"))
            return range;
        else if(m.equals("nopri"))
            return nopri;
        else return all;

    }

    private static String[] nopri=new String []{   
        "miRNA_start","miRNA_end","upperStem_start","upperStem_end",
        "preRNA_size","lowerStem_size","upperStem_size","topStem_size",
        "preRNA_energy","miRNA_stability",
        "preRNA_GC_content","miRNA_GC_content",
        "preRNA_A_content","preRNA_U_content","preRNA_G_content","preRNA_C_content",
        "miRNA_A_content","miRNA_U_content","miRNA_G_content","miRNA_C_content",
        "preRNA_pair_number","miRNA_pair_number",
        "preRNA_unpair_number","preRNA_unpair_rate",
        "miRNA_unpair_number","miRNA_unpair_rate",
        "lowerStem_unpair_number","lowerStem_unpair_rate",
        "topStem_unpair_number","topStem_unpair_rate",
        "preRNA_G-U_number","miRNA_G-U_number",
        "preRNA_internalLoop_number","preRNA_internalLoop_size",
        "miRNA_internalLoop_number","miRNA_internalLoop_size",
        "lowerStem_internalLoop_number","lowerStem_internalLoop_size",
        "topStem_internalLoop_number","topStem_internalLoop_size",
        "miRNA_firstBase","overhang_base1","overhang_base2"
    };//43

    private static String[] range=new String []{
        "priRNA_size","preRNA_size","miRNA_size",
        "miRNA_start","miRNA_end","upperStem_start","upperStem_end",
        "basalSegment_size","lowerStem_size","upperStem_size","topStem_size","terminalLoop_size",
        "priRNA_energy","preRNA_energy","miRNA_stability",
        "priRNA_GC_content","preRNA_GC_content","miRNA_GC_content",
        "priRNA_A_content","priRNA_U_content","priRNA_G_content","priRNA_C_content",
        "preRNA_A_content","preRNA_U_content","preRNA_G_content","preRNA_C_content",
        "miRNA_A_content","miRNA_U_content","miRNA_G_content","miRNA_C_content",
        "priRNA_pair_number","preRNA_pair_number","miRNA_pair_number",
        "priRNA_unpair_number","priRNA_unpair_rate",
        "preRNA_unpair_number","preRNA_unpair_rate",
        "miRNA_unpair_number","miRNA_unpair_rate",
        "lowerStem_unpair_number","lowerStem_unpair_rate",
        "topStem_unpair_number","topStem_unpair_rate",
        "priRNA_G-U_number","preRNA_G-U_number","miRNA_G-U_number",
        "priRNA_internalLoop_number","priRNA_internalLoop_size",
        "preRNA_internalLoop_number","preRNA_internalLoop_size",
        "miRNA_internalLoop_number","miRNA_internalLoop_size",
        "lowerStem_internalLoop_number","lowerStem_internalLoop_size",
        "topStem_internalLoop_number","topStem_internalLoop_size"};//56

    //my feature combination
//    private static final String[] animal=new String[]{
//        "preRNA_G_content","preRNA_internalLoop_number","miRNA_stability","miRNA_U_content",
//        "miRNA_internalLoop_number","overhang_base1","topStem_size","topStem_internalLoop_size",
//        "upperStem_start","preRNA_size","preRNA_GC_content","preRNA_pair_number",
//        "miRNA_pair_number","miRNA_G-U_number","overhang_base2","lowerStem_unpair_rate",
//        "lowerStem_internalLoop_number","topStem_unpair_number","topStem_unpair_rate",
//        "topStem_internalLoop_number","upperStem_size","preRNA_unpair_rate"
//    };//22
//
//    private static final String[] overall=new String[]{
//        "preRNA_energy","preRNA_A_content","preRNA_U_content","preRNA_C_content",
//        "preRNA_unpair_rate","preRNA_internalLoop_number","miRNA_stability","miRNA_G_content",
//        "miRNA_GC_content","miRNA_unpair_number","miRNA_firstBase","lowerStem_internalLoop_size",
//        "topStem_size","topStem_internalLoop_size","preRNA_size","preRNA_GC_content",
//        "preRNA_pair_number","miRNA_pair_number","miRNA_G-U_number","overhang_base2",
//        "lowerStem_unpair_rate","lowerStem_internalLoop_number","topStem_unpair_number",
//        "topStem_unpair_rate","topStem_internalLoop_number","upperStem_size"
//    };//26
//
//    private static final String[] plant=new String []{
//        "preRNA_A_content","preRNA_G_content","preRNA_U_content","preRNA_C_content",
//        "preRNA_G-U_number","preRNA_unpair_rate","miRNA_G_content","miRNA_U_content",
//        "miRNA_C_content","miRNA_GC_content","miRNA_internalLoop_number","miRNA_firstBase",
//        "overhang_base1","lowerStem_unpair_number","upperStem_start","upperStem_end",
//        "preRNA_size","preRNA_GC_content","preRNA_pair_number","miRNA_pair_number",
//        "miRNA_G-U_number","overhang_base2","lowerStem_unpair_rate","lowerStem_internalLoop_number",
//        "topStem_unpair_number","topStem_unpair_rate","topStem_internalLoop_number",
//        "upperStem_size"
//    };//28
//
//    private static final String[] virus=new String []{
//        "preRNA_size","preRNA_energy","preRNA_A_content","preRNA_G_content",
//        "preRNA_U_content","preRNA_C_content","preRNA_GC_content","preRNA_unpair_number",
//        "preRNA_unpair_rate","preRNA_internalLoop_size","miRNA_stability",
//        "miRNA_A_content","miRNA_G_content","miRNA_U_content","miRNA_GC_content",
//        "miRNA_unpair_rate","lowerStem_size","lowerStem_unpair_rate",
//        "lowerStem_internalLoop_number","topStem_size","topStem_unpair_number",
//        "topStem_unpair_rate","topStem_internalLoop_size","topStem_internalLoop_number",
//        "upperStem_size","upperStem_end"
//    };//26


    private static final String[] all=new String []{
        "priRNA_id","miRNA_id",

        "priRNA_sequence","preRNA_sequence","miRNA_sequence",

        "priRNA_start","priRNA_end","priRNA_size",
        "preRNA_size",
        "miRNA_start","miRNA_end","miRNA_size",
        "basalSegment_size","basalSegment_end",
        "lowerStem_size",
        "upperStem_start","upperStem_end","upperStem_size",
        "topStem_size",
        "terminalLoop_size",

        "priRNA_structure","priRNA_energy","preRNA_structure","preRNA_energy",

        "priRNA_GC_content","preRNA_GC_content","miRNA_GC_content",
        "priRNA_A_content","priRNA_U_content","priRNA_G_content","priRNA_C_content",
        "preRNA_A_content","preRNA_U_content","preRNA_G_content","preRNA_C_content",
        "miRNA_A_content","miRNA_U_content","miRNA_G_content","miRNA_C_content",

        "priRNA_pair_number","preRNA_pair_number","miRNA_pair_number",
        "priRNA_unpair_number","priRNA_unpair_rate",
        "preRNA_unpair_number","preRNA_unpair_rate",
        "miRNA_unpair_number","miRNA_unpair_rate",
        "lowerStem_unpair_number","lowerStem_unpair_rate",
        "topStem_unpair_number","topStem_unpair_rate",
        "priRNA_G-U_number","preRNA_G-U_number","miRNA_G-U_number",
        "priRNA_internalLoop_number","priRNA_internalLoop_size",
        "preRNA_internalLoop_number","preRNA_internalLoop_size",
        "miRNA_internalLoop_number","miRNA_internalLoop_size",
        "lowerStem_internalLoop_number","lowerStem_internalLoop_size",
        "topStem_internalLoop_number","topStem_internalLoop_size",

        "miRNA_stability","strand","miRNA_firstBase",
        "overhang_base1","overhang_base2"
    };//70


    //Wu's feature combination
        private static final String[] overall=new String []{
        "lowerStem_unpair_rate", "preRNA_size", "lowerStem_internalLoop_number", "upperStem_size",
        "miRNA_internalLoop_size", "miRNA_firstBase", "miRNA_A_content", "miRNA_GC_content",
        "preRNA_GC_content", "preRNA_U_content", "preRNA_A_content", "topStem_internalLoop_size",
        "preRNA_unpair_rate", "miRNA_pair_number", "preRNA_pair_number", "topStem_internalLoop_number",
        "topStem_unpair_rate", "miRNA_C_content", "miRNA_start", "miRNA_unpair_rate",
        "miRNA_G_content", "overhang_base2", "preRNA_C_content", "preRNA_G_content",
        "miRNA_U_content"};//25

    private static final String[] animal=new String []{
        "topStem_internalLoop_size", "topStem_internalLoop_number", "topStem_size", "upperStem_size",
        "miRNA_GC_content", "miRNA_internalLoop_size", "miRNA_internalLoop_number", "miRNA_A_content",
        "miRNA_C_content", "miRNA_G_content", "miRNA_U_content", "miRNA_start",
        "miRNA_unpair_rate", "overhang_base1", "preRNA_GC_content", "preRNA_size",
        "preRNA_A_content", "preRNA_C_content", "preRNA_G_content", "preRNA_U_content",
        "preRNA_unpair_rate", "miRNA_stability", "lowerStem_unpair_rate", "topStem_unpair_rate"};//24

    private static final String[] virus=animal;

    private static final String[] plant=new String []{
        "miRNA_firstBase", "topStem_internalLoop_size", "lowerStem_internalLoop_size", "lowerStem_internalLoop_number",
        "upperStem_size", "miRNA_GC_content", "miRNA_G-U_number", "miRNA_internalLoop_size",
        "miRNA_internalLoop_number", "miRNA_A_content", "miRNA_G_content", "miRNA_U_content",
        "miRNA_pair_number", "overhang_base1", "preRNA_GC_content", "preRNA_A_content",
        "preRNA_C_content", "preRNA_G_content", "preRNA_U_content", "preRNA_unpair_rate",
        "miRNA_stability", "lowerStem_unpair_rate", "topStem_unpair_rate", "upperStem_start"};//24

    //new optimum
//    private static final String[] overall=new String []{
//        "preRNA_size","miRNA_G_content","miRNA_firstBase","preRNA_A_content","miRNA_C_content","preRNA_U_content","preRNA_unpair_rate","topStem_unpair_rate","miRNA_unpair_rate","lowerStem_unpair_rate","preRNA_G_content"};//11
//
//    private static final String[] animal=new String []{
//        "preRNA_size","miRNA_firstBase","miRNA_A_content","preRNA_A_content","miRNA_U_content","miRNA_GC_content","preRNA_G_content","topStem_unpair_rate","miRNA_C_content","miRNA_G_content","lowerStem_unpair_rate","miRNA_unpair_rate","preRNA_unpair_rate"};//13
//
//    private static final String[] virus=new String []{
//        "preRNA_size","miRNA_firstBase","miRNA_U_content","topStem_unpair_rate","preRNA_GC_content","preRNA_A_content","miRNA_G_content",
//"miRNA_C_content","miRNA_unpair_rate","preRNA_unpair_rate","preRNA_A_content","miRNA_internalLoop_number"
//    };//12
//
//    private static final String[] plant=new String []{
//        "miRNA_unpair_number","miRNA_firstBase","miRNA_internalLoop_number","miRNA_G_content","miRNA_internalLoop_size","preRNA_G_content","topStem_unpair_rate","miRNA_A_content","preRNA_GC_content",
//"miRNA_U_content","miRNA_unpair_rate","miRNA_C_content","preRNA_A_content","preRNA_U_content","miRNA_GC_content","preRNA_unpair_rate","preRNA_C_content"};//17
}
