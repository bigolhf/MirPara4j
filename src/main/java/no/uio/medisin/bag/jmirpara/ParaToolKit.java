/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import no.uio.medisin.bag.core.MfeFoldRNA;
import no.uio.medisin.bag.core.SimpleRNASequence;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author weibo
 */
public class ParaToolKit {

//    /**
//     * compare two number
//     * and return the greater one
//     * @param int a
//     * @param int b
//     * @return
//     */
//    public static int max(int a, int b) {
//        if (a >= b) {
//            return a;
//        } else {
//            return b;
//        }
//    }
//
    public static void foldRNA(SimpleRNASequence rna){
        rna.setEnergy(MfeFoldRNA.foldSequence(rna.getSeq(),rna.getStructureStr()));
        rna.setStructureString(MfeFoldRNA.getStructure());
    }
//
//    public static ArrayList<PriMiRNA> extractHairpin(SimRNA rna, int threshold){
//        ArrayList priList=new ArrayList<PriMiRNA>();
//        String str=rna.getStr();
//
//        String regSL="([\\(\\.]*\\()(\\.+)(\\)[\\)\\.]*)";
//        Pattern pattern=Pattern.compile(regSL);
//        Matcher m=pattern.matcher(str);
//
//        int count=1;
//
//        while(m.find()){
//
//            int start1=str.lastIndexOf(")", m.start(1))+1;//replace m.start(1)
//            String left=str.substring(start1, m.end(1));//the 5' side of the stemloop
//            String right=m.group(3);//the 3' side of the stemloop
//            int n=b2Num(left)+b2Num(right);//compare the two sides of the stem
//            int start=0,end=0,l=0,r=0;
//            //if str is like (..((...)).. return ..((...))..
//            if(n>0){
//                l=bIndex(left,"(",n)+1; //new start of left
//                left=left.substring(l); //new leftht
//                start=start1+l; //start of stemloop, count from 0
//                end=m.end(3);//count from 1
//            }
//            //if str is like ..((...))..) return ..((...))..
//            else if(n<0){
//                r=bIndex(right,")",n)-1;
//                right=right.substring(0,r+1);
//                start=start1;
//                end=m.start(3)+r+1;//count from 1
//            }
//            //if str is like ..((...)).. return ..((...))..
//            else{
//                start=start1;
//                end=m.end(3);//count from 1
//            }
//            String subId=rna.getId()+"-mir-"+count; //the id of the stemloop
//            String subSeq=rna.getSeq().substring(start, end); //seq of the stemloop
//            String subStr=left+m.group(2)+right; //structure of the stemloop
//
//            //filter the SL less than the threshold
//            if(subStr.length()<threshold)
//                continue;
//
//            //create a new primiRNA
//            PriMiRNA pri=new PriMiRNA(subId,subSeq);
//
//            //if the stemloop has two or more loop
//            if(hasTwoLoop(pri)) continue;
//
//            pri.setName(rna.getName());
//            pri.setStart(start+rna.getStart());
//            pri.setEnd(end-1+rna.getEnd());
//            
//            priList.add(pri);
//
//            count++;
//        }
//        return priList;
//    }
//
//    /**
//     * transform bracket-dot string to number and return the sum
//     * @param String str: structure with bracket-dot notation string
//     * @return int: the sum of the str( each '(' is 1, ')' is -1, '.' is 0
//     */
//    private static int b2Num(String str){
//        int num=0;
//        for(int i=0;i<str.length();i++){
//            if(str.charAt(i)=='(')
//                num+=1;
//            else if(str.charAt(i)==')')
//                num-=1;
//
//        }
//        return num;
//    }
//    /**
//     * find the index of the nth '(' or ')'
//     * @param String p: the structure string
//     * @param String s: "(" or ")"
//     * @param int n: the '(' or ')' number
//     * @return
//     */
//    private static int bIndex(String p,String s,int n){
//        int m=Math.abs(n);
//        int c=0;
//        if(s.equals("(")){
//            for(int i=0;i<m;i++){
//                c=p.indexOf(s,c)+1;
//            }
//            c=c-1;
//        }
//        if(s.equals(")")){
//            c=p.length()-1;
//            for(int i=0;i<m;i++){
//                c=p.lastIndexOf(s, c)-1;
//            }
//            c=c+1;
//        }
//
//        return c;
//    }
//
//    /**
//     * judge whether the structure of a sequence has two or more loops
//     * @param String seq: sequence to be tested
//     * @return boolean; if have two or more return true, or false
//     */
//    private static boolean hasTwoLoop(SimRNA rna){
//        foldRNA(rna);
//        int end5=rna.getStr().lastIndexOf("(");
//        int start3=rna.getStr().indexOf(")");
//        if(end5>=start3)
//            return true;
//        else return false;
//    }
//
//    public static void parsePriStr(PriMiRNA priRNA){
//        String reg = "(\\.*)(\\([\\.\\(]*\\()(\\.+)(\\)[\\.\\)]*\\))(\\.*)";
//        Pattern pattern = Pattern.compile(reg);
//        Matcher m = pattern.matcher(priRNA.getStr());
//
//        if (m.matches()) {
//            priRNA.setBasalSegSize( ParaToolKit.max(m.group(1).length(), m.group(5).length()));
//            priRNA.setBasalSegEnd(priRNA.getBasalSegSize());//count from 1
//            priRNA.setBasalBaseNum(m.group(1).length() + m.group(5).length());
//
//            priRNA.setTerminalLoopSeq(priRNA.getSeq().substring(m.start(3), m.end(3)));
//            priRNA.setTerminalLoopSize(priRNA.getTerminalLoopSeq().length());
//            priRNA.setTerminalLoopStart( m.start(3)+1);//count from 1
//            priRNA.setTerminalLoopEnd(m.end(3));//count from 1
//
//
//            int end5 = m.end(2) + (m.end(3) - m.start(3)) / 2; //5' end, count from 1
//            int start3 = m.start(4) - (m.end(3) - m.start(3)) / 2; //3' start, count from 0
//            priRNA.setMidBase("");
//            if (end5 < start3)
//                priRNA.setMidBase(priRNA.getSeq().charAt(end5)+"");//top base on the terminal loop
//                        //include terminal loop
//            priRNA.setSeq5(priRNA.getSeq().substring(0, end5));
//            priRNA.setStr5(priRNA.getStr().substring(0, end5));
//            priRNA.setSeq3(priRNA.getSeq().substring(start3));
//            priRNA.setStr3(priRNA.getStr().substring(start3));
//
//            priRNA.setStrIndex(str2Plot(priRNA));
//        }
//    }
//
//    
//    public static int[] str2Plot(PriMiRNA priRNA) {
//        String priLine1="";
//        String priLine2="";
//        String priLine3="";
//        String priLine4="";
//        String seq5=priRNA.getSeq5();
//        String str5=priRNA.getStr5();
//        String seq3=reverse(priRNA.getSeq3());
//        String str3=reverse(priRNA.getStr3());
//        int size=priRNA.getLength();
//
//        int i = 0, j = 0, n = 0;
//        int[] index = new int[size];
//        while (i < seq5.length() && j < seq3.length()) {
//            if (str5.charAt(i) == '.') {
//                priLine1 += seq5.charAt(i);
//                priLine2 += ' ';
//                priLine3 += ' ';
//                index[i] = n;
//                i += 1;
//                if (str3.charAt(j) == '.') {
//                    priLine4 += seq3.charAt(j);
//                    index[size - 1 - j] = n;
//                    j += 1;
//                } else {
//                    priLine4 += ' ';
//                }
//            } else {
//                priLine1 += ' ';
//                if (str3.charAt(j) == '.') {
//                    priLine2 += ' ';
//                    priLine3 += ' ';
//                    priLine4 += seq3.charAt(j);
//                    index[size - 1 - j] = n;
//                    j += 1;
//                } else {
//                    priLine2 += seq5.charAt(i);
//                    priLine3 += seq3.charAt(j);
//                    priLine4 += ' ';
//                    index[i] = n;
//                    index[size - 1 - j] = n;
//                    i += 1;
//                    j += 1;
//                }
//            }
//            n++;
//        }
//
//        priRNA.setPriLine1(priLine1);
//        priRNA.setPriLine2(priLine2);
//        priRNA.setPriLine3(priLine3);
//        priRNA.setPriLine4(priLine4);
//
//        return index;
//    }
//
//    public static ArrayList<HashMap> getFeatureMap(ArrayList<PriMiRNA> priList){
//        ArrayList featureList=new ArrayList<HashMap>();
//        HashMap features=new HashMap();
//        for(PriMiRNA pri : priList){
//            features.putAll(pri.getFeatureSet());
//            features.put("plot", pri.getPriPlot());
//            for(int i=0;i<pri.sizeOfProduct();i++){
//                PreMiRNA pre=pri.nextPreRNA();
//                features.putAll(pre.getFeatureSet());
//                for(int j=0;j<pre.SizeOfProduct();j++){
//                    MatMiRNA mi=pre.NextMiRNA();
//                    features.putAll(mi.getFeatureSet());
//                    featureList.add((HashMap)features.clone());
//
//                }
//            }
//        }
//        return featureList;
//    }
//
//    /**
//     * reverse a string
//     * @param String str
//     * @return String reverse str
//     */
//    public static String reverse(String str) {
//        String reStr = "";
//        for (int i = str.length() - 1; i >= 0; i--) {
//            reStr += str.charAt(i);
//        }
//        return reStr;
//    }
//
//    /**
//     * transform structure plot to seq
//     * @param String s1: the first line of one side of the plot
//     * @param String s2: the second line of one side of the plot
//     * @return
//     */
//    public static String plot2Seq(String s1, String s2) {
//        String seq = "";
//        int n=s1.length();
//        char[] ss=new char[n*2];
//        for (int i = 0; i < n; i++) {
////            seq = seq + s1.charAt(i) + s2.charAt(i);
//            ss[i*2]=s1.charAt(i);
//            ss[i*2+1]=s2.charAt(i);
//        }
//        seq=new String(ss);
//        seq = seq.replace(" ", "");
//        
//        return seq;
//    }
//
//    /**
//     * calculate the GC content
//     * @param String seq
//     * @return float: the GC content of the seq
//     */
//    public static float contentGC(String seq) {
//        String s = seq.replaceAll("[CGcg]", "");
//        if (seq.length() == 0) {
//            return 0;
//        }
//        return 1 - (float) s.length() / seq.length();
//
//    }
//
//    /**
//     * calculate the nt(a,t,g,c,u) content
//     * @param String seq
//     * @param String nt: a,t,g,c,u
//     * @return float: the content of a nt
//     */
//    public static float contentNT(String seq, String nt) {
//        String s = seq.replaceAll(nt, "");
//        if (seq.length() == 0) {
//            return 0;
//        }
//        return 1 - (float) s.length() / seq.length();
//    }
//
//    /**
//     * calculate the GU-pair number
//     * @param String priLine2
//     * @param String priLine3
//     * @return int: the number of GU-pairs
//     */
//    public static int pairGUnum(String line2, String line3) {
//        String s = line2 + line3;
//        int nc = s.replaceAll("[Cc]", "").length();
//        int ng = s.replaceAll("[Gg]", "").length();
//        return nc - ng;
//    }
//
//    /**
//     * calculate the pair number
//     * @param String line2
//     * @return int: the number of pairs
//     */
//    public static int pairNum(String line2) {
//        return line2.replaceAll("\\s+", "").length();
//
//    }
//
//    /**
//     * calculate the unpaired bases number
//     * @param String line1
//     * @param String line4
//     * @return int: the number of unpaired bases
//     */
//    public static int unpairedNum(String line1, String line4) {
//        return (line1 + line4).replaceAll("\\s+", "").length();
//    }
//
//    /**
//     * calculate the unpaired rate
//     * @param int unpaired: the number of unpaired bases
//     * @param int paired: the number of pairs
//     * @return float: the rate of unpaired bases rate
//     */
//    public static float unpairedRate(int unpaired, int paired) {
//        int sum = unpaired + paired * 2;
//        if (sum == 0) {
//            return 0;
//        }
//        return (float) unpaired / sum;
//    }
//
//    /**
//     * calculate the internal loop number
//     * @param String line2
//     * @return int: the number of internal loop number
//     */
//    public static int internalLoopNum(String line2) {
//        return line2.replaceAll("\\s+", " ").split(" ", -1).length - 1;
//    }
//
//    /**
//     * calculate the maximal internal loop size
//     * @param String line1
//     * @param String line4
//     * @return int: the size of the maximal internal loop
//     */
//    public static int internalLoopSize(String line1, String line4) {
//        String[] s = (line1 + line4).replaceAll("\\s+", " ").split(" ");
//        int loop = 0;
//        for (int i = 0; i < s.length; i++) {
//            loop = max(loop, s[i].length());
//        }
//        return loop;
//    }
//
//    /**
//     * decide the strand of miRNA on primiRNA
//     * @param priRNA
//     * @param start
//     * @return
//     */
//    public static int strand(PriMiRNA priRNA, int start){
//        if(start+1<priRNA.getSeq5().length())
//            return 5;
//        else if(start+1>priRNA.getSeq5().length()+priRNA.getMidBase().length())
//            return 3;
//        else return 0;
//    }
//
//    /**
//     * generate the 5' dangle of miRNA
//     * @param String line1
//     * @param String line2
//     * @param int strand: the strand of miRNA
//     * @return String: two bases sequence refers to the 5' dangle
//     */
//    public static String dangleSeq(String seq, String line1, String line2, int strand) {
//        String line = plot2Seq(line1, line2);
//        if (strand == 5) {
//            line = reverse(line);
//        }
//        int end = seq.indexOf(line);
//        if (end == 0) {
//            return "  ";
//        } else if (end == 1) {
//            return seq.charAt(0) + " ";
//        } else {
//            return seq.substring(end - 2, end);
//        }
//    }
//
//    /**
//     * calculate miRNA stability
//     * @param String line2
//     * @param String line3
//     * @param int strand: the strand of miRNA
//     * @return float: the stability of miRNA
//     */
//    public static float stability(String line2, String line3, int strand) {  
//        int u = (line2.substring(0, 4) + line3.substring(0, 4)).replaceAll("\\s+", "").length();
//        int d = (line2.substring(line2.length() - 4) + line3.substring(line3.length() - 4)).replaceAll("\\s+", "").length();
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
//
//    /**
//     * calculate pre-miRNA energy
//     * @param priSeq
//     * @param priStr
//     * @param preSeq
//     * @return
//     */
//    public static float preMFE(String priStr, String preSeq, int start) {
//        int end = start + preSeq.length();
//        String preStr = priStr.substring(start, end);
//        return MfeFold.cal(preSeq, preStr);
//    }
//
//    
//
//    
//
//    public static String setStrPlot(PriMiRNA priRNA, int upperStart, int upperEnd, int strand){
//        String line1="",line2="",line3="",line4="",line5="";
//        String l1=priRNA.getPriLine1();
//        String l2=priRNA.getPriLine2();
//        String l3=priRNA.getPriLine3();
//        String l4=priRNA.getPriLine4();
//        if(strand==5){
//
//
//            if(upperEnd+1<l1.length()-1){
//                line1=l1.substring(0, upperStart).toLowerCase()+
//                        l1.substring(upperStart, upperEnd+1).toUpperCase()+
//                        l1.substring(upperEnd+1, l1.length()-1).toLowerCase();
//                line2=l2.substring(0, upperStart).toLowerCase()+
//                        l2.substring(upperStart, upperEnd+1).toUpperCase()+
//                        l2.substring(upperEnd+1, l2.length()-1).toLowerCase()+
//                        l1.substring(l1.length()-1).toLowerCase();
//            }
//            line4=l3.substring(0,l3.length()-1).toLowerCase()+
//                    l4.substring(l4.length()-1).toLowerCase();
//            line5=l4.substring(0, l4.length()-1).toLowerCase();
//
//        }
//        else{
//            line1=l1.substring(0, l1.length()-1).toLowerCase();
//            line2=l2.substring(0,l2.length()-1).toLowerCase()+
//                    l1.substring(l1.length()-1).toLowerCase();
//            if(upperEnd+1<l1.length()-1){
//                line4=l3.substring(0, upperStart).toLowerCase()+
//                        l3.substring(upperStart, upperEnd+1).toUpperCase()+
//                        l3.substring(upperEnd+1, l3.length()-1).toLowerCase()+
//                        l4.substring(l4.length()-1).toLowerCase();
//                line5=l4.substring(0, upperStart).toLowerCase()+
//                        l4.substring(upperStart, upperEnd+1).toUpperCase()+
//                        l4.substring(upperEnd+1, l4.length()-1).toLowerCase();
//            }
//
//        }
//        line3=l2.substring(0, l2.length()-1).replaceAll("\\S", "|")+priRNA.getMidBase().toLowerCase();
//        
//        return line1+"\n"+line2+"\n"+line3+"\n"+line4+"\n"+line5;
//    }

}
