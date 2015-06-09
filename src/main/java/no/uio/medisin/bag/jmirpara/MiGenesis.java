
package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * this class include methods to parse the features of pri-miRNA,pre-miRNA and miRNA.
 * the order of parsing:
 * parsePrimiRNA();
 * maturateMiRNA();
 * gatherFeatures();
 * @author weibo
 */
public class MiGenesis {

    private PriMiRNA priRNA;
    private PreMiRNA preRNA;
    private MatMiRNA miRNA;

    private String priLine1;
    private String priLine2;
    private String priLine3;
    private String priLine4;



    private static String regSL="(\\.*)(\\([\\.\\(]*\\()(\\.+)(\\)[\\.\\)]*\\))(\\.*)";
                                // (1)........(2).........(3).......(4)..........(5)
    private Pattern pattern = Pattern.compile(regSL);
    private Matcher m;

    public MiGenesis(PriMiRNA pri){
        priRNA=pri;
    }

    /**
     * parse features of a pri-miRNA
     */
    public void parsePrimiRNA(){
//        if(priRNA.getStr().equals("") || priRNA.getEnergy()==0)
//            SLScaner.foldRNA(priRNA);

        parsePriStr(priRNA);

        priLine1=priRNA.getPriLine1();
        priLine2=priRNA.getPriLine2();
        priLine3=priRNA.getPriLine3();
        priLine4=priRNA.getPriLine4();
        //nucleotide content
        priRNA.setGC_content(contentGC(priRNA.getSeq()));
//        priRNA.setA_content(contentNT(priRNA.getSeq(), "[Aa]"));
//        priRNA.setU_content(contentNT(priRNA.getSeq(), "[Uu]"));
//        priRNA.setG_content(contentNT(priRNA.getSeq(), "[Gg]"));
//        priRNA.setC_content(contentNT(priRNA.getSeq(), "[Cc]"));
        priRNA.setA_content(contentNT(priRNA.getSeq(), 'A')+contentNT(priRNA.getSeq(), 'a'));
        priRNA.setU_content(contentNT(priRNA.getSeq(), 'U')+contentNT(priRNA.getSeq(), 'u'));
        priRNA.setG_content(contentNT(priRNA.getSeq(), 'G')+contentNT(priRNA.getSeq(), 'g'));
        priRNA.setC_content(contentNT(priRNA.getSeq(), 'C')+contentNT(priRNA.getSeq(), 'c'));
        //GU pair number
        priRNA.setGU_num(pairGUnum(priLine2, priLine3));
        //pair number
        priRNA.setPair_num(pairNum(priLine2));
        //priRNA unpaired base do not include basal part and terminalloop
        priRNA.setUnpairedBase_num(unpairedNum(priLine1,priLine4) - priRNA.getBasalBaseNum() -
                (priRNA.getTerminalLoopSeq().length() - priRNA.getMidBase().length()));
        priRNA.setUnpairedBase_rate(unpairedRate(
                priRNA.getUnpairedBase_num(), priRNA.getPair_num()));
        //priRNA internal loops DO NOT include basal part, do not include terminal loop
        priRNA.setInternalLoop_num(internalLoopNum(priLine2));
        priRNA.setInternalLoopSize(internalLoopSize( priLine1,priLine4));

        priRNA.setFeatureSet();
    }

    /**
     * maturate miRNA from a pri-miRNA.
     * parse the features of pre-miRNA and miRNA
     * @param priRNA
     * @param start--miRNA start, count from 0
     * @param size--miRNA length
     */
    public void maturateMiRNA(int start, int size){

        preRNA = new PreMiRNA();
        miRNA = new MatMiRNA();

        int strand=strand(priRNA, start, size);
        int upperStart=0,upperEnd=0;
        int[] strIndex=priRNA.getStrIndex();
        
        String midBase=priRNA.getMidBase();
        int basalEnd=priRNA.getBasalSegEnd();//count from 1

        priLine1=priRNA.getPriLine1();
        priLine2=priRNA.getPriLine2();
        priLine3=priRNA.getPriLine3();
        priLine4=priRNA.getPriLine4();
        String preLine1;
        String preLine2;
        String preLine3;
        String preLine4;
        String miLine1;
        String miLine2;
        String miLine3;
        String miLine4;
        String lowerLine1;
        String lowerLine2;
        String lowerLine3;
        String lowerLine4;
        String topLine1;
        String topLine2;
        String topLine3;
        String topLine4;
        
        String dangle;

        if (strand == 5) {
            upperStart = strIndex[start];//count from 0
            upperEnd = strIndex[start + size - 1];//count from 0

        } else if (strand == 3) {
            upperStart = strIndex[start + size - 1];
            upperEnd = strIndex[start];
        }
        else return;//do something...

        

        priRNA.setPriPlot(setStrPlot(priRNA,upperStart,upperEnd, strand));

        //---------------set pre-miRNA
        preLine1 = priLine1.substring(upperStart);
        preLine2 = priLine2.substring(upperStart);
        preLine3 = priLine3.substring(upperStart);
        preLine4 = priLine4.substring(upperStart);
        preRNA.setUpperStart(upperStart + 1);//count from 1
        preRNA.setUpperEnd(upperEnd + 1);//count from 1
        preRNA.setSeq(plot2Seq(preLine1, preLine2) + midBase + reverse(plot2Seq(preLine4, preLine3)));
        preRNA.setLength(preRNA.getSeq().length());
        preRNA.setStart(plot2Seq(priLine1.substring(0, upperStart), priLine2.substring(0, upperStart)).length());
        preRNA.setEnergy(preMFE(priRNA.getStr(), preRNA.getSeq(), preRNA.getStart()));
        preRNA.setGC_content(contentGC(preRNA.getSeq()));
//        preRNA.setA_content(contentNT(preRNA.getSeq(), "[Aa]"));
//        preRNA.setU_content(contentNT(preRNA.getSeq(), "[Uu]"));
//        preRNA.setG_content(contentNT(preRNA.getSeq(), "[Gg]"));
//        preRNA.setC_content(contentNT(preRNA.getSeq(), "[Cc]"));
        preRNA.setA_content(contentNT(preRNA.getSeq(),'A')+contentNT(preRNA.getSeq(),'a'));
        preRNA.setU_content(contentNT(preRNA.getSeq(),'U')+contentNT(preRNA.getSeq(),'u'));
        preRNA.setG_content(contentNT(preRNA.getSeq(),'G')+contentNT(preRNA.getSeq(),'g'));
        preRNA.setC_content(contentNT(preRNA.getSeq(),'C')+contentNT(preRNA.getSeq(),'c'));
        //GU pair number
        preRNA.setGU_num(pairGUnum(preLine2, preLine3));
        //pair number
        preRNA.setPair_num(pairNum(preLine2));
        //preRNA unpaired base, do not include terminal loop bases
        preRNA.setUnpairedBase_num(unpairedNum(preLine1, preLine4) - (priRNA.getTerminalLoopSize() - midBase.length()));
        preRNA.setUnpairedBase_rate(unpairedRate(preRNA.getUnpairedBase_num(), preRNA.getPair_num()));
        //preRNA internal loop, do not include terminal loop
        preRNA.setInternalLoop_num(internalLoopNum(preLine2));
        preRNA.setInternalLoopSize(internalLoopSize(preLine1, preLine4));

        //set lowerStem
        if (basalEnd < upperStart) {
            lowerLine1 = priLine1.substring(basalEnd, upperStart);
            lowerLine2 = priLine2.substring(basalEnd, upperStart);
            lowerLine3 = priLine3.substring(basalEnd, upperStart);
            lowerLine4 = priLine4.substring(basalEnd, upperStart);
            //lower stem length
            preRNA.setLowerStemSize(upperStart-basalEnd);
            //lower stem unpaired base
            preRNA.setLowerStemUnpairedBase_num(unpairedNum(lowerLine1, lowerLine4));
            preRNA.setLowerStemUnpairedBase_rate(unpairedRate(preRNA.getLowerStemUnpairedBase_num(), pairNum(lowerLine2)));
            //lower stem internal loop
            preRNA.setLowerStemInternalLoop_num(internalLoopNum(lowerLine2));
            preRNA.setLowerStemInternalLoopSize(internalLoopSize(lowerLine1, lowerLine4));
        }

        preRNA.setUpperStemSize(upperEnd - upperStart + 1);

        //set topStem
        int stemEnd = strIndex[priRNA.getTerminalLoopStart()-1];//count from 1
        if (upperEnd + 1 < stemEnd) {
            topLine1 = priLine1.substring(upperEnd + 1, stemEnd);
            topLine2 = priLine2.substring(upperEnd + 1, stemEnd);
            topLine3 = priLine3.substring(upperEnd + 1, stemEnd);
            topLine4 = priLine4.substring(upperEnd + 1, stemEnd);
            //top stem length
            preRNA.setTopStemSize(stemEnd - upperEnd - 1);
            //top stem unpaired base
            preRNA.setTopStemUnpairedBase_num(unpairedNum(topLine1, topLine4));
            preRNA.setTopStemUnpairedBase_rate(unpairedRate(preRNA.getTopStemUnpairedBase_num(), pairNum(topLine2)));
            //top stem internal loop
            preRNA.setTopStemInternalLoop_num(internalLoopNum(topLine2));
            preRNA.setTopStemInternalLoopSize(internalLoopSize(topLine1, topLine4));
        }
        preRNA.setFeatureSet();


        //-----------------set miRNA
        miLine1 = priLine1.substring(upperStart, upperEnd + 1);
        miLine2 = priLine2.substring(upperStart, upperEnd + 1);
        miLine3 = priLine3.substring(upperStart, upperEnd + 1);
        miLine4 = priLine4.substring(upperStart, upperEnd + 1);
        
        miRNA.setSeq(priRNA.getSeq().substring(start, start + size));
        miRNA.setLength(miRNA.getSeq().length());
        miRNA.setGC_content(contentGC(miRNA.getSeq()));
//        miRNA.setA_content(contentNT(miRNA.getSeq(), "[Aa]"));
//        miRNA.setU_content(contentNT(miRNA.getSeq(), "[Uu]"));
//        miRNA.setG_content(contentNT(miRNA.getSeq(), "[Gg]"));
//        miRNA.setC_content(contentNT(miRNA.getSeq(), "[Cc]"));
        miRNA.setA_content(contentNT(miRNA.getSeq(),'A')+contentNT(miRNA.getSeq(),'a'));
        miRNA.setU_content(contentNT(miRNA.getSeq(),'U')+contentNT(miRNA.getSeq(),'u'));
        miRNA.setG_content(contentNT(miRNA.getSeq(),'G')+contentNT(miRNA.getSeq(),'g'));
        miRNA.setC_content(contentNT(miRNA.getSeq(),'C')+contentNT(miRNA.getSeq(),'c'));
        //GU pair number
        miRNA.setGU_num(pairGUnum(miLine2, miLine3));
        //pair number
        miRNA.setPair_num(pairNum(miLine2));
        //miRNA unpaired base
        miRNA.setUnpairedBase_num(unpairedNum(miLine1, miLine4));
        miRNA.setUnpairedBase_rate(unpairedRate(miRNA.getUnpairedBase_num(), miRNA.getPair_num()));
        //miRNA internal loop
        miRNA.setInternalLoop_num(internalLoopNum(miLine2));
        miRNA.setInternalLoopSize(internalLoopSize(miLine1, miLine4));
        //first base
        miRNA.setFirstBase(priRNA.getSeq().charAt(start));
        miRNA.setMiStart(start + 1);//count from 1
        miRNA.setMiEnd(start + size);//count from 1
        miRNA.setStart(start + priRNA.getStart());//count from 1
        miRNA.setEnd(miRNA.getStart() + size - 1);//count from 1
        miRNA.setLength(size);
        miRNA.setStrand(strand);
        //miRNA dangle
        if (strand == 5) {
            dangle = dangleSeq(priLine3.length()-upperStart-1,reverse(priLine4),reverse(priLine3));
        } else {
            dangle = dangleSeq(upperEnd,priLine1,priLine2);
        }

        miRNA.setDangleBaseOne(dangle.charAt(0));
        miRNA.setDangleBaseTwo(dangle.charAt(1));
        //stability
        miRNA.setStability(stability(miLine2, miLine3, strand));
        //miRNA id
        miRNA.setName(priRNA.getName());
        miRNA.setId(miRNA.getName() + "_MIR_" + miRNA.getStart()
                + "-" + miRNA.getLength());

        miRNA.setFeatureSet();

        //store product in priRNA
//        preRNA.addProduct(miRNA);
//        priRNA.addProduct(preRNA);
    }

    //###################tools#####################################

    /**
     * parse the hairpin structure of a pri-miRNA
     * @param priRNA
     */
    public boolean parsePriStr(PriMiRNA priRNA){
        
        m = pattern.matcher(priRNA.getStr());

        if (m.matches()) {
            priRNA.setBasalSegSize( Math.max(m.group(1).length(), m.group(5).length()));
            priRNA.setBasalSegEnd(priRNA.getBasalSegSize());//count from 1
            priRNA.setBasalBaseNum(m.group(1).length() + m.group(5).length());

            priRNA.setTerminalLoopSeq(priRNA.getSeq().substring(m.start(3), m.end(3)));
            priRNA.setTerminalLoopSize(priRNA.getTerminalLoopSeq().length());
            priRNA.setTerminalLoopStart(m.start(3)+1);//count from 1
            priRNA.setTerminalLoopEnd(m.end(3));//count from 1


            int end5 = m.end(2) + (m.end(3) - m.start(3)) / 2; //5' end, count from 1
            int start3 = m.start(4) - (m.end(3) - m.start(3)) / 2; //3' start, count from 0
            priRNA.setMidBase("");
            if (end5 < start3)
                priRNA.setMidBase(priRNA.getSeq().charAt(end5)+"");//top base on the terminal loop
                        //include terminal loop
            priRNA.setSeq5(priRNA.getSeq().substring(0, end5));
            priRNA.setStr5(priRNA.getStr().substring(0, end5));
            priRNA.setSeq3(priRNA.getSeq().substring(start3));
            priRNA.setStr3(priRNA.getStr().substring(start3));

            priRNA.setStrIndex(str2Plot(priRNA));
            
            return true;
        }
        return false;
    }

    /**
     * transform bracket-dot notation string structure to text plot and
     * store index of the bases on the plot
     * @param priRNA
     * @return
     */
    public int[] str2Plot(PriMiRNA priRNA) {
        StringBuilder pril1=new StringBuilder();
        StringBuilder pril2=new StringBuilder();
        StringBuilder pril3=new StringBuilder();
        StringBuilder pril4=new StringBuilder();
        String seq5=priRNA.getSeq5();
        String str5=priRNA.getStr5();
        String seq3=reverse(priRNA.getSeq3());
        String str3=reverse(priRNA.getStr3());
        int size=priRNA.getLength();

        int i = 0, j = 0, n = 0;
        int[] index = new int[size];
        while (i < seq5.length() && j < seq3.length()) {
            if (str5.charAt(i) == '.') {
                pril1.append(seq5.charAt(i));
                pril2.append(' ');
                pril3.append(' ');
                index[i] = n;
                i += 1;
                if (str3.charAt(j) == '.') {
                    pril4.append(seq3.charAt(j));
                    index[size - 1 - j] = n;
                    j += 1;
                } else {
                    pril4.append(' ');
                }
            } else {
                pril1.append(' ');
                if (str3.charAt(j) == '.') {
                    pril2.append(' ');
                    pril3.append(' ');
                    pril4.append(seq3.charAt(j));
                    index[size - 1 - j] = n;
                    j += 1;
                } else {
                    pril2.append(seq5.charAt(i));
                    pril3.append(seq3.charAt(j));
                    pril4.append(' ');
                    index[i] = n;
                    index[size - 1 - j] = n;
                    i += 1;
                    j += 1;
                }
            }
            n++;
        }

        priRNA.setPriLine1(pril1.toString());
        priRNA.setPriLine2(pril2.toString());
        priRNA.setPriLine3(pril3.toString());
        priRNA.setPriLine4(pril4.toString());

        return index;
    }

    /**
     * reverse a string
     * @param String str
     * @return String reverse str
     */
    public String reverse(String str) {
        StringBuilder reStr = new StringBuilder();
        for (int i = str.length() - 1; i >= 0; i--) {
            reStr.append(str.charAt(i));
        }
        return reStr.toString();
    }

    /**
     * transform structure plot to seq
     * @param String s1: the first line of one side of the plot
     * @param String s2: the second line of one side of the plot
     * @return
     */
    public String plot2Seq(String s1, String s2) {
        StringBuilder seq = new StringBuilder();

        for (int i = 0; i < s1.length(); i++) {
            seq.append(s1.charAt(i)).append(s2.charAt(i));
            
        }
        return seq.toString().replace(" ", "");
    }

    /**
     * calculate the GC content
     * @param String seq
     * @return float: the GC content of the seq
     */
//    public float contentGC(String seq) {
//        String s = seq.replaceAll("[CGcg]", "");
//        if (seq.length() == 0) {
//            return 0;
//        }
//        return 1 - (float) s.length() / seq.length();
//
//    }
    public float contentGC(String seq){
        return contentNT(seq,'G')+contentNT(seq,'C')+contentNT(seq,'g')+contentNT(seq,'c');
    }

    /**
     * calculate the nt(a,t,g,c,u) content
     * @param String seq
     * @param String nt: a,t,g,c,u
     * @return float: the content of a nt
     */
//    public float contentNT(String seq, String nt) {
//        String s = seq.replaceAll(nt, "");
//        if (seq.length() == 0) {
//            return 0;
//        }
//        return 1 - (float) s.length() / seq.length();
//    }
    public float contentNT(String seq, char nt){       
        return (float)numOfNT(seq,nt)/seq.length();
    }
    private int numOfNT(String seq, char nt){
        char[] cs=seq.toCharArray();
        int n=0;
        for(char c:cs){
            if(c==nt){
                n++;
            }
        }
        return n;
    }

    /**
     * calculate the GU-pair number
     * @param String priLine2
     * @param String priLine3
     * @return int: the number of GU-pairs
     */
//    public int pairGUnum(String line2, String line3) {
//        String s = line2 + line3;
//        int nc = s.replaceAll("[Cc]", "").length();
//        int ng = s.replaceAll("[Gg]", "").length();
//        return nc - ng;
//    }
    public int pairGUnum(String line2, String line3){
        int n=line2.length()+line3.length();
        int nc=n-(numOfNT(line2,'C')+numOfNT(line2,'c')+numOfNT(line3,'C')+numOfNT(line3,'c'));
        int ng=n-(numOfNT(line2,'G')+numOfNT(line2,'g')+numOfNT(line3,'G')+numOfNT(line3,'g'));
        return nc-ng;
    }

    /**
     * calculate the pair number
     * @param String line2
     * @return int: the number of pairs
     */
//    public int pairNum(String line2) {
//        return line2.replaceAll("\\s+", "").length();
//
//    }
    public int pairNum(String line2){
        return line2.length()-numOfNT(line2,' ');
    }

    /**
     * calculate the unpaired bases number
     * @param String line1
     * @param String line4
     * @return int: the number of unpaired bases
     */
//    public int unpairedNum(String line1, String line4) {
//        return (line1 + line4).replaceAll("\\s+", "").length();
//    }
    public int unpairedNum(String line1, String line4){
        return line1.length()-numOfNT(line1,' ')+line4.length()-numOfNT(line4,' ');
    }

    /**
     * calculate the unpaired rate
     * @param int unpaired: the number of unpaired bases
     * @param int paired: the number of pairs
     * @return float: the rate of unpaired bases rate
     */
    public float unpairedRate(int unpaired, int paired) {
        int sum = unpaired + paired * 2;
        if (sum == 0) {
            return 0;
        }
        return (float) unpaired / sum;
    }

    /**
     * calculate the internal loop number
     * @param String line2
     * @return int: the number of internal loop number
     */
//    public int internalLoopNum(String line2) {
//        String l=line2.trim().replaceAll("\\s+", " ");
//        return l.length()-l.replace(" ", "").length();
//    }
    public int internalLoopNum(String line2){
        char[] l=line2.toCharArray();
        int n=0;
        int size=l.length;
        for(int i=0;i<size-1;i++){
            if(l[i]!=' ' && l[i+1]==' '){
                n+=1;
            }
        }
        if(l[size-1]==' '){
            n-=1;
        }
        return n;
    }

    /**
     * calculate the maximal internal loop size
     * @param String line1
     * @param String line4
     * @return int: the size of the maximal internal loop
     */
//    public int internalLoopSize(String line1, String line4) {
//        String l1=line1.replaceAll("^\\w+\\s+","").replaceAll("\\s+\\w+$", " ");
//        String l4=line4.replaceAll("^\\w+\\s+"," ").replaceAll("\\s+\\w+$", "");
//        String[] s = (l1 + l4).replaceAll("\\s+", " ").split(" ");
//        int loop = 0;
//        for (int i = 0; i < s.length; i++) {
//            loop = Math.max(loop, s[i].length());
//        }
//        return loop;
//    }
    public int internalLoopSize(String line1, String line4){
        return Math.max(maxGapFreeInternalSubString(line1), 
                maxGapFreeInternalSubString(line4));
    }
    private int maxGapFreeInternalSubString(String str){
        char[] c1=str.toCharArray();
        int n=c1.length;
        int b=0,e=0;
        for(int i=0;i<n;i++){
            if(c1[i]==' '){
                b=i;
                break;
            }
        }
        for(int i=n-1;i>0;i--){
            if(c1[i]==' '){
                e=i;
                break;
            }
        }
        int max=0;
        int s=0;
        for(int i=b;i<e;i++){
            if(c1[i]==' ' && c1[i+1]!=' '){
                s=1;
            }
            else if(c1[i]!=' ' && c1[i+1]!=' '){
                s+=1;
            }
            else if(c1[i]!=' ' && c1[i+1]==' '){
                if(s>max){
                    max=s;
                }
            }
        }
        return max;
    }

    /**
     * decide the strand of miRNA on primiRNA
     * @param priRNA
     * @param start--miRNA start position, count from 0
     * @param size--miRNA length
     * @return
     */
    public static int strand(PriMiRNA priRNA, int start, int size){
        if(start+size<=priRNA.getSeq5().length())
            return 5;
        else if(start+1>priRNA.getSeq5().length()+priRNA.getMidBase().length() && start+size<=priRNA.getLength())
            return 3;
        else return 0;
    }

    /**
     * generate the 3' 2nt overhang
     * @param end--miRNA end position
     * @param line1
     * @param line2
     * @return
     */
    public String dangleSeq(int end, String line1, String line2) {
        String dangle="";
        int i=1;
        StringBuilder d=new StringBuilder();
        while(dangle.length()<2 && end+i<line1.length()){
//            dangle=dangle+line1.charAt(end+i)+line2.charAt(end+i);
            d.append(line1.charAt(end+i)).append(line2.charAt(end+i));
            dangle=d.toString().replaceAll("\\s+", "");
            i++;
        }
        if(dangle.length()==1)
            dangle+=" ";
        else if(dangle.length()==0)
            dangle="  ";
        return dangle;
    }

    /**
     * calculate miRNA stability(rate of the number of hydrogen bonds at two sides of miRNA)
     * @param String line2
     * @param String line3
     * @param int strand: the strand of miRNA
     * @return float: the stability of miRNA
     */
//    public float stability(String line2, String line3, int strand) {
//        int u = (line2.substring(0, 4) + line3.substring(0, 4)).replaceAll("\\s+", "").replaceAll("[cC]", "xx").length();
//        int d = (line2.substring(line2.length() - 4) + line3.substring(line3.length() - 4)).replaceAll("\\s+", "").replaceAll("[cC]", "xx").length();
//        if (d == 0 || u == 0) {
//            return 0;
//        }
//        if (strand == 5) {
//            return (float) u / d;
//        } else if (strand == 3) {
//            return (float) d / u;
//        } else {
//            return 0;
//        }
//    }
    public float stability(String line2, String line3, int strand){
        int u=0;
        char[] l2=line2.toCharArray();
        char[] l3=line3.toCharArray();
        for(int i=0;i<4;i++){
            if(l2[i]==' ') continue;
            else if(l2[i]=='c') u+=2;
            else if(l2[i]=='C') u+=2;
            else u+=1;
        }
        for(int i=0;i<4;i++){
            if(l3[i]==' ') continue;
            else if(l3[i]=='c') u+=2;
            else if(l3[i]=='C') u+=2;
            else u+=1;
        }
        int d=0;
        int n2=l2.length;
        int n3=l3.length;
        for(int i=1;i<5;i++){
            if(l2[n2-i]==' ') continue;
            else if(l2[n2-i]=='c') d+=2;
            else if(l2[n2-i]=='C') d+=2;
            else d+=1;
        }
        for(int i=1;i<5;i++){
            if(l3[n3-i]==' ') continue;
            else if(l3[n3-i]=='c') d+=2;
            else if(l3[n3-i]=='C') d+=2;
            else d+=1;
        }
        if (d == 0 || u == 0) {
            return 0;
        }
        if (strand == 5) {
            return (float) u / d;
        } else if (strand == 3) {
            return (float) d / u;
        } else {
            return 0;
        }
    }

    /**
     * calculate pre-miRNA energy
     * @param priSeq
     * @param priStr
     * @param preSeq
     * @return
     */
    public float preMFE(String priStr, String preSeq, int start) {
        int end = start + preSeq.length();
        String preStr = priStr.substring(start, end);
        return MfeFold.cal(preSeq, preStr);
    }

    /**
     * pri-miRNA hairpin structure
     * @param priRNA
     * @return
     */
    public static String setStrPlot(PriMiRNA priRNA){
//        String line1="",line2="",line3="",line4="",line5="";
        String priLine1=priRNA.getPriLine1();
        String priLine2=priRNA.getPriLine2();
        String priLine3=priRNA.getPriLine3();
        String priLine4=priRNA.getPriLine4();
        StringBuilder sp=new StringBuilder();
//        line1=priLine1.substring(0, priLine1.length()-1);
        sp.append(priLine1.substring(0, priLine1.length()-1)).append("\n");
//        line2=priLine2.substring(0,priLine3.length()-1)+priLine1.substring(priLine1.length()-1);
        sp.append(priLine2.substring(0,priLine3.length()-1)).append(priLine1.substring(priLine1.length()-1)).append("\n");
//        line3=priLine2.substring(0, priLine2.length()-1).replaceAll("\\S", "|")+priRNA.getMidBase();
        sp.append(priLine2.substring(0, priLine2.length()-1).replaceAll("\\S", "|")).append(priRNA.getMidBase()).append("\n");
//        line4=priLine3.substring(0,priLine3.length()-1)+priLine4.substring(priLine4.length()-1);
        sp.append(priLine3.substring(0,priLine3.length()-1)).append(priLine4.substring(priLine4.length()-1)).append("\n");
//        line5=priLine4.substring(0, priLine4.length()-1);
        sp.append(priLine4.substring(0, priLine4.length()-1));
        
//        return line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n"+line5;
        return sp.toString();
    }


    /**
     * pri-miRNA hairpin structure, miRNA region is highlighted by uppercase letter
     * @param priRNA
     * @param upperStart
     * @param upperEnd
     * @param strand
     * @return
     */
    public String setStrPlot(PriMiRNA priRNA, int upperStart, int upperEnd, int strand){
//        String line1="",line2="",line3="",line4="",line5="";
        StringBuilder sp=new StringBuilder();
        if(strand==5){
            if(upperEnd+1<priLine1.length()-1){
//                line1=priLine1.substring(0, upperStart).toLowerCase()+
//                        priLine1.substring(upperStart, upperEnd+1).toUpperCase()+
//                        priLine1.substring(upperEnd+1, priLine1.length()-1).toLowerCase();
                sp.append(priLine1.substring(0, upperStart).toLowerCase())
                        .append(priLine1.substring(upperStart, upperEnd+1).toUpperCase())
                        .append(priLine1.substring(upperEnd+1, priLine1.length()-1).toLowerCase()).append("\n");
//                line2=priLine2.substring(0, upperStart).toLowerCase()+
//                        priLine2.substring(upperStart, upperEnd+1).toUpperCase()+
//                        priLine2.substring(upperEnd+1, priLine2.length()-1).toLowerCase()+
//                        priLine1.substring(priLine1.length()-1).toLowerCase();
                sp.append(priLine2.substring(0, upperStart).toLowerCase())
                        .append(priLine2.substring(upperStart, upperEnd+1).toUpperCase())
                        .append(priLine2.substring(upperEnd+1, priLine2.length()-1).toLowerCase())
                        .append(priLine1.substring(priLine1.length()-1).toLowerCase()).append("\n");
            }
            sp.append(priLine2.substring(0, priLine2.length()-1).replaceAll("\\S", "|"))
                    .append(priRNA.getMidBase().toLowerCase()).append("\n");
//            line4=priLine3.substring(0,priLine3.length()-1).toLowerCase()+
//                    priLine4.substring(priLine4.length()-1).toLowerCase();
            sp.append(priLine3.substring(0,priLine3.length()-1).toLowerCase())
                    .append(priLine4.substring(priLine4.length()-1).toLowerCase()).append("\n");
//            line5=priLine4.substring(0, priLine4.length()-1).toLowerCase();
            sp.append(priLine4.substring(0, priLine4.length()-1).toLowerCase());

        }
        else{
//            line1=priLine1.substring(0, priLine1.length()-1).toLowerCase();
            sp.append(priLine1.substring(0, priLine1.length()-1).toLowerCase()).append("\n");
//            line2=priLine2.substring(0,priLine2.length()-1).toLowerCase()+
//                    priLine1.substring(priLine1.length()-1).toLowerCase();
            sp.append(priLine2.substring(0,priLine2.length()-1).toLowerCase())
                    .append(priLine1.substring(priLine1.length()-1).toLowerCase()).append("\n");
            sp.append(priLine2.substring(0, priLine2.length()-1).replaceAll("\\S", "|"))
                    .append(priRNA.getMidBase().toLowerCase()).append("\n");
            if(upperEnd+1<priLine1.length()-1){
//                line4=priLine3.substring(0, upperStart).toLowerCase()+
//                        priLine3.substring(upperStart, upperEnd+1).toUpperCase()+
//                        priLine3.substring(upperEnd+1, priLine3.length()-1).toLowerCase()+
//                        priLine4.substring(priLine4.length()-1).toLowerCase();
                sp.append(priLine3.substring(0, upperStart).toLowerCase())
                        .append(priLine3.substring(upperStart, upperEnd+1).toUpperCase())
                        .append(priLine3.substring(upperEnd+1, priLine3.length()-1).toLowerCase())
                        .append(priLine4.substring(priLine4.length()-1).toLowerCase()).append("\n");
//                line5=priLine4.substring(0, upperStart).toLowerCase()+
//                        priLine4.substring(upperStart, upperEnd+1).toUpperCase()+
//                        priLine4.substring(upperEnd+1, priLine4.length()-1).toLowerCase();
                sp.append(priLine4.substring(0, upperStart).toLowerCase())
                        .append(priLine4.substring(upperStart, upperEnd+1).toUpperCase())
                        .append(priLine4.substring(upperEnd+1, priLine4.length()-1).toLowerCase());
            }

        }
//        line3=priLine2.substring(0, priLine2.length()-1).replaceAll("\\S", "|")+priRNA.getMidBase().toLowerCase();

//        return line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n"+line5;
        return sp.toString();
    }

   

    /**
     * gather all the features of a pri-miRNA and its pre-miRNA and miRNA 
     * to store in a HashMap
     * @return
     */
    public HashMap gatherFeatures(){

        HashMap features=new HashMap();

        features.putAll(priRNA.getFeatureSet());
        features.put("plot", priRNA.getPriPlot());

        features.putAll(preRNA.getFeatureSet());

        features.putAll(miRNA.getFeatureSet());


        return features;
    }


}
