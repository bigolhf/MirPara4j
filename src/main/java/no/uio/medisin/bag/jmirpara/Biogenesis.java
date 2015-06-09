
package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;

/**
 *
 * @author weibo
 */
public class Biogenesis {


//    public ArrayList<PriMiRNA> findPriRNA(SimRNA seq, int threshold){
//        ParaToolKit.foldRNA(seq);
//        return ParaToolKit.extractHairpin(seq, threshold);
//    }
//
//
//
//    public void processPriRNA(PriMiRNA priRNA){
//        if(priRNA.getStr().equals("") || priRNA.getEnergy()==0)
//            ParaToolKit.foldRNA(priRNA);
//        parsePri(priRNA);
//    }
//    
//    private void parsePri(PriMiRNA priRNA){
//        ParaToolKit.parsePriStr(priRNA);
//
//        String priLine1=priRNA.getPriLine1();
//        String priLine2=priRNA.getPriLine2();
//        String priLine3=priRNA.getPriLine3();
//        String priLine4=priRNA.getPriLine4();
//
//        priRNA.setGC_content(ParaToolKit.contentGC(priRNA.getSeq()));
//        priRNA.setA_content(ParaToolKit.contentNT(priRNA.getSeq(), "[Aa]"));
//        priRNA.setU_content(ParaToolKit.contentNT(priRNA.getSeq(), "[Uu]"));
//        priRNA.setG_content(ParaToolKit.contentNT(priRNA.getSeq(), "[Gg]"));
//        priRNA.setC_content(ParaToolKit.contentNT(priRNA.getSeq(), "[Cc]"));
//        //GU pair number
//        priRNA.setGU_num(ParaToolKit.pairGUnum(priLine2, priLine3));
//        //pair number
//        priRNA.setPair_num(ParaToolKit.pairNum(priLine2));
//        //priRNA unpaired base do not include basal part and terminalloop
//        priRNA.setUnpairedBase_num(ParaToolKit.unpairedNum(priLine1,priLine4) - priRNA.getBasalBaseNum() -
//                (priRNA.getTerminalLoopSeq().length() - priRNA.getMidBase().length()));
//        priRNA.setUnpairedBase_rate(ParaToolKit.unpairedRate(
//                priRNA.getUnpairedBase_num(), priRNA.getPair_num()));
//        //priRNA internal loop include basal part, do not include terminal loop
//        priRNA.setInternalLoop_num(ParaToolKit.internalLoopNum(priLine2) - 1);
//        priRNA.setInternalLoopSize(ParaToolKit.internalLoopSize(
//                priLine1.replaceAll("\\s+\\w+$", " "), priLine4.replaceAll("\\s+\\w+$", " ")));
//
//        priRNA.setFeatureSet();
//    }
//
//    public int numOfmir(PriMiRNA priRNA){
//        return 5*(priRNA.getTerminalLoopStart()-(20+24)+priRNA.getLength());
//    }
//    
//
//    public void findMiRNA(PriMiRNA priRNA){
//        for(int w=20;w<25;w++){//miRNA size
//            int i;
//           //5' strand
//            for(i=0;i<priRNA.getTerminalLoopStart()-w;i++) //miRNA start point
//                maturateMiRNA(priRNA,i,w);
//
//            //3' strand
//            for(i=priRNA.getTerminalLoopEnd()+1;i<priRNA.getLength()-w;i++)
//                maturateMiRNA(priRNA,i,w);
//                
//        }
//    }
//    
//
//    public void maturateMiRNA(PriMiRNA priRNA, int start, int size){
//        PreMiRNA preRNA = new PreMiRNA();
//        MatMiRNA miRNA = new MatMiRNA();
//
//        int strand=ParaToolKit.strand(priRNA, start);
//        int upperStart=0,upperEnd=0;
//        int[] strIndex=priRNA.getStrIndex();
//        String priLine1=priRNA.getPriLine1();
//        String priLine2=priRNA.getPriLine2();
//        String priLine3=priRNA.getPriLine3();
//        String priLine4=priRNA.getPriLine4();
//        String midBase=priRNA.getMidBase();
//        int basalEnd=priRNA.getBasalSegEnd();//count from 1
//        String preLine1,preLine2,preLine3,preLine4,
//                lowerLine1,lowerLine2,lowerLine3,lowerLine4,
//                topLine1,topLine2,topLine3,topLine4,
//                miLine1,miLine2,miLine3,miLine4;
//        String dangle;
//
//        if (strand == 5) {
//            upperStart = strIndex[start];//count from 0
//            upperEnd = strIndex[start + size - 1];//count from 0
//
//        } else if (strand == 3) {
//            upperStart = strIndex[start + size - 1];
//            upperEnd = strIndex[start];
//        }
//
//        priRNA.setPriPlot(ParaToolKit.setStrPlot(priRNA,upperStart,upperEnd, strand));
//        
//        //---------------set pre-miRNA
//        preLine1 = priLine1.substring(upperStart);
//        preLine2 = priLine2.substring(upperStart);
//        preLine3 = priLine3.substring(upperStart);
//        preLine4 = priLine4.substring(upperStart);
//        preRNA.setUpperStart(upperStart + 1);//count from 1
//        preRNA.setUpperEnd(upperEnd + 1);//count from 1
//        preRNA.setSeq(ParaToolKit.plot2Seq(preLine1, preLine2) + midBase + ParaToolKit.reverse(ParaToolKit.plot2Seq(preLine4, preLine3)));
//        preRNA.setLength(preRNA.getSeq().length());
//        preRNA.setStart(ParaToolKit.plot2Seq(priLine1.substring(0, upperStart), priLine2.substring(0, upperStart)).length());
//        preRNA.setEnergy(ParaToolKit.preMFE(priRNA.getStr(), preRNA.getSeq(), preRNA.getStart()));
//        preRNA.setGC_content(ParaToolKit.contentGC(preRNA.getSeq()));
//        preRNA.setA_content(ParaToolKit.contentNT(preRNA.getSeq(), "[Aa]"));
//        preRNA.setU_content(ParaToolKit.contentNT(preRNA.getSeq(), "[Uu]"));
//        preRNA.setG_content(ParaToolKit.contentNT(preRNA.getSeq(), "[Gg]"));
//        preRNA.setC_content(ParaToolKit.contentNT(preRNA.getSeq(), "[Cc]"));
//        //GU pair number
//        preRNA.setGU_num(ParaToolKit.pairGUnum(preLine2, preLine3));
//        //pair number
//        preRNA.setPair_num(ParaToolKit.pairNum(preLine2));
//        //preRNA unpaired base, do not include terminal loop bases
//        preRNA.setUnpairedBase_num(ParaToolKit.unpairedNum(preLine1, preLine4) - (priRNA.getTerminalLoopSize() - midBase.length()));
//        preRNA.setUnpairedBase_rate(ParaToolKit.unpairedRate(preRNA.getUnpairedBase_num(), preRNA.getPair_num()));
//        //preRNA internal loop, do not include terminal loop
//        preRNA.setInternalLoop_num(ParaToolKit.internalLoopNum(preLine2) - 1);
//        preRNA.setInternalLoopSize(ParaToolKit.internalLoopSize(preLine1.replaceAll("\\s+\\w+$", " "), preLine4.replaceAll("\\s+\\w+$", " ")));
//
//        //set lowerStem
//        if (basalEnd < upperStart) {
//            lowerLine1 = priLine1.substring(basalEnd + 1, upperStart);
//            lowerLine2 = priLine2.substring(basalEnd + 1, upperStart);
//            lowerLine3 = priLine3.substring(basalEnd + 1, upperStart);
//            lowerLine4 = priLine4.substring(basalEnd + 1, upperStart);
//            //lower stem length
//            preRNA.setLowerStemSize(upperStart - basalEnd - 1);
//            //lower stem unpaired base
//            preRNA.setLowerStemUnpairedBase_num(ParaToolKit.unpairedNum(lowerLine1, lowerLine4));
//            preRNA.setLowerStemUnpairedBase_rate(ParaToolKit.unpairedRate(preRNA.getLowerStemUnpairedBase_num(), ParaToolKit.pairNum(lowerLine2)));
//            //lower stem internal loop
//            preRNA.setLowerStemInternalLoop_num(ParaToolKit.internalLoopNum(lowerLine2));
//            preRNA.setLowerStemInternalLoopSize(ParaToolKit.internalLoopSize(lowerLine1, lowerLine4));
//        }
//
//        preRNA.setUpperStemSize(upperEnd - upperStart + 1);
//
//        //set topStem
//        int stemEnd = strIndex[priRNA.getTerminalLoopStart()-1];//count from 1
//        if (upperEnd + 1 < stemEnd) {
//            topLine1 = priLine1.substring(upperEnd + 1, stemEnd);
//            topLine2 = priLine2.substring(upperEnd + 1, stemEnd);
//            topLine3 = priLine3.substring(upperEnd + 1, stemEnd);
//            topLine4 = priLine4.substring(upperEnd + 1, stemEnd);
//            //top stem length
//            preRNA.setTopStemSize(stemEnd - upperEnd - 1);
//            //top stem unpaired base
//            preRNA.setTopStemUnpairedBase_num(ParaToolKit.unpairedNum(topLine1, topLine4));
//            preRNA.setTopStemUnpairedBase_rate(ParaToolKit.unpairedRate(preRNA.getTopStemUnpairedBase_num(), ParaToolKit.pairNum(topLine2)));
//            //top stem internal loop
//            preRNA.setTopStemInternalLoop_num(ParaToolKit.internalLoopNum(topLine2));
//            preRNA.setTopStemInternalLoopSize(ParaToolKit.internalLoopSize(topLine1, topLine4));
//        }
//        preRNA.setFeatureSet();
//
//
//        //-----------------set miRNA
//        miLine1 = priLine1.substring(upperStart, upperEnd + 1);
//        miLine2 = priLine2.substring(upperStart, upperEnd + 1);
//        miLine3 = priLine3.substring(upperStart, upperEnd + 1);
//        miLine4 = priLine4.substring(upperStart, upperEnd + 1);
//        miRNA.setSeq(priRNA.getSeq().substring(start, start + size));
//        miRNA.setLength(miRNA.getSeq().length());
//        miRNA.setGC_content(ParaToolKit.contentGC(miRNA.getSeq()));
//        miRNA.setA_content(ParaToolKit.contentNT(miRNA.getSeq(), "[Aa]"));
//        miRNA.setU_content(ParaToolKit.contentNT(miRNA.getSeq(), "[Uu]"));
//        miRNA.setG_content(ParaToolKit.contentNT(miRNA.getSeq(), "[Gg]"));
//        miRNA.setC_content(ParaToolKit.contentNT(miRNA.getSeq(), "[Cc]"));
//        //GU pair number
//        miRNA.setGU_num(ParaToolKit.pairGUnum(miLine2, miLine3));
//        //pair number
//        miRNA.setPair_num(ParaToolKit.pairNum(miLine2));
//        //miRNA unpaired base
//        miRNA.setUnpairedBase_num(ParaToolKit.unpairedNum(miLine1, miLine4));
//        miRNA.setUnpairedBase_rate(ParaToolKit.unpairedRate(miRNA.getUnpairedBase_num(), miRNA.getPair_num()));
//        //miRNA internal loop
//        miRNA.setInternalLoop_num(ParaToolKit.internalLoopNum(miLine2));
//        miRNA.setInternalLoopSize(ParaToolKit.internalLoopSize(miLine1, miLine4));
//        //first base
//        miRNA.setFirstBase(priRNA.getSeq().charAt(start));
//        miRNA.setMiStart(start + 1);//count from 1
//        miRNA.setMiEnd(start + size);//count from 1
//        miRNA.setStart(start + priRNA.getStart());//count from 1
//        miRNA.setEnd(miRNA.getStart() + size - 1);//count from 1
//        miRNA.setLength(size);
//        miRNA.setStrand(strand);
//        //miRNA dangle
//        if (strand == 5) {
//            dangle = ParaToolKit.dangleSeq(priRNA.getSeq(),miLine4, miLine3, strand);
//        } else {
//            dangle = ParaToolKit.dangleSeq(priRNA.getSeq(),miLine1, miLine2, strand);
//        }
//        miRNA.setDangleBaseOne(dangle.charAt(1));
//        miRNA.setDangleBaseTwo(dangle.charAt(0));
//        //stability
//        miRNA.setStability(ParaToolKit.stability(miLine2, miLine3, strand));
//        //miRNA id
//        miRNA.setName(priRNA.getName());
//        miRNA.setId(miRNA.getName() + "_MIR_" + miRNA.getStart()
//                + "-" + miRNA.getLength());
//
//        miRNA.setFeatureSet();
////
////        //store product in priRNA
////        preRNA.addProduct(miRNA);
////        priRNA.addProduct(preRNA);
//    }



}
