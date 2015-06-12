
package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * this class finds hairpins from a long sequence.
 * the order of calling methods after instantiated:
 * slidingWindow();
 * scanStemLoop();
 * noRepeat();
 * @author weibo
 */
public class StemLoopScanner {

    private int window=500;
    private int step=250;
    private int start=1;
    private int distance = 60;

    private SimpleSeq sequence;
    
    private ArrayList<SimpleSeq> fragmentList;
    private ArrayList<PriMiRNA> primiRNAList;

//    private static String regSL="([\\(\\.]*\\()(\\.+)(\\)[\\)\\.]*)";
    private static String regSL="(\\.*\\([\\.\\(]*\\()(\\.+)(\\)[\\.\\)]*\\)\\.*)";
    private static Pattern pattern = Pattern.compile(regSL);
    private Matcher m;

    
    /**
     * Empty Class Constructor
     */
    public StemLoopScanner(){

    }
    
    
    
    /**
     * Class constructor from input sequence
     * 
     * @param seq 
     */
    public StemLoopScanner(SimpleSeq seq){
        sequence = seq;
    }
    
    
    
    /**
     * Class constructor from subsequence of input sequence
     * 
     * @param seq
     * @param window
     * @param step
     * @param distance 
     */
    public StemLoopScanner(SimpleSeq seq, int window, int step, int distance){
        
        this.sequence = seq;
        this.window = window;
        this.step = step;
        this.distance = distance;
    }
    

    /**
     * 
     * slide window with given step size to generate a serial of segments for testing.
     * before calling, one may set the sliding window size, step increment and start position
     */
    public void slidingWindow(){
        fragmentList=new ArrayList<SimpleSeq>();
        int length=sequence.getLength();
        int n=(length-window)/step;
        if(n<0) n=0;
        if(length-n*step+window>19) n+=1;
        int end=0;start-=1;
        for(int i=0;i<=n;i++){
            
            if(start>=length) break;
            end = start + window;
            if(end>length) end=length;
            String id = sequence.getName() + "_" + (start+1) + "-"+end;
            String subseq = sequence.getSeq().substring(start,end);

            SimpleSeq frag = new SimpleSeq(id,subseq);
            frag.setStart(start+1);// count from 1
            frag.setEnd(end); //count from 1
            frag.setName(sequence.getId());
            fragmentList.add(frag);
            start+= step;
            
        }
        
    }
    

    /**
     * generate candidate pri-miRNAs from all the segments
     */
    public void scanStemLoop(){
        primiRNAList=new ArrayList<PriMiRNA>();
        
        int i=1;
        for(SimpleSeq seg : fragmentList){
            primiRNAList.addAll(scanStemLoopFrom(seg));
            System.out.print(Output.decimal(i*100.0/fragmentList.size())+"%"+Output.backspace(Output.decimal(i*100.0/fragmentList.size())+"%"));
            i++;
        }
//        segList=null; // release seglist space
    }
    
    
    public ArrayList<PriMiRNA> scanStemLoopFrom(SimpleSeq seq){
        SimpleRNASequence rna = new SimpleRNASequence(seq);
        foldRNA(rna);
        return extractHairpin(rna);
    }

    /**
     * fold a rna sequence
     * @param rna
     */
    public static void foldRNA(SimpleRNASequence rna){
        rna.setEnergy(MfeFoldRNA.fold(rna.getSeq(), rna.getStr()));
        rna.setStr(MfeFoldRNA.getStructure());

//        jmipara.MfeFold mf=new jmipara.MfeFold();
//        mf.setSequence(rna.getSeq());
//        mf.setStructure("");
//        mf.cal();
//        rna.setEnergy(mf.getEnergy());
//        rna.setStr(mf.getStructure());
    }

    /**
     * extract hairpins from a rna secondary structure.
     * setCutOff() can change the minimal length of the hairpin
     * @param rna
     */
    public ArrayList<PriMiRNA> extractHairpin(SimpleRNASequence rna){
        
        ArrayList<PriMiRNA> pris=new ArrayList<PriMiRNA>();
        
        String str=rna.getStr(); 
        m=pattern.matcher(str);

        while(m.find()){

            int start1=str.lastIndexOf(")", m.start(1))+1;//replace m.start(1)
            String left=str.substring(start1, m.end(1));//the 5' side of the stemloop
            String right=m.group(3);//the 3' side of the stemloop
            int n=b2Num(left)+b2Num(right);//compare the two sides of the stem
            int slStart=0,slEnd=0,l=0,r=0;
            //if str is like (..((...)).. return ..((...))..
            if(n>0){
                l=bIndex(left,"(",n)+1; //new start of left
                left=left.substring(l); //new leftht
                slStart=start1+l; //start of stemloop, count from 0
                slEnd=m.end(3);//count from 1
            }
            //if str is like ..((...))..) return ..((...))..
            else if(n<0){
                r=bIndex(right,")",n)-1;
                right=right.substring(0,r+1);
                slStart=start1;
                slEnd=m.start(3)+r+1;//count from 1
            }
            //if str is like ..((...)).. return ..((...))..
            else{
                slStart=start1;
                slEnd=m.end(3);//count from 1
            }
            String subId=rna.getName()+"_mir_"; //the id of the stemloop
            String subSeq=rna.getSeq().substring(slStart, slEnd); //seq of the stemloop
            String subStr=left+m.group(2)+right; //structure of the stemloop

            //filter the SL less than the threshold
            if(subStr.length()<distance)
                continue;

            //create a new primiRNA
            PriMiRNA pri=new PriMiRNA(subId,subSeq);

            //if the stemloop has two or more loop
            if(hasTwoLoop(pri)) continue;
            if(pattern.matcher(pri.getStr()).matches() ==false) continue;
            
            pri.setName(rna.getName());
            pri.setStart(slStart+rna.getStart());
            pri.setEnd(slEnd-1+rna.getStart());
            pri.setId(pri.getId()+pri.getStart()+"-"+pri.getEnd());
            
            pris.add(pri);

//            priList.add(pri);          
            
        }
        return pris;
    }

    /**
     * remove the repeat hairpin caused by the overlap region during sliding window
     */
    public void noRepeat(){
        HashMap map=new HashMap();
        for(PriMiRNA pri:primiRNAList)
            map.put(pri.getId(), pri);

        primiRNAList=new ArrayList<PriMiRNA>(map.values());
    }

    /**
     * transform bracket-dot string to number and return the sum
     * @param String str: structure with bracket-dot notation string
     * @return int: the sum of the str( each '(' is 1, ')' is -1, '.' is 0
     */
    private int b2Num(String str){
        int num=0;
        for(int i=0;i<str.length();i++){
            if(str.charAt(i)=='(')
                num+=1;
            else if(str.charAt(i)==')')
                num-=1;

        }
        return num;
    }
    /**
     * find the index of the nth '(' or ')'
     * @param String p: the structure string
     * @param String s: "(" or ")"
     * @param int n: the '(' or ')' number
     * @return
     */
    private int bIndex(String p,String s,int n){
        int m=Math.abs(n);
        int c=0;
        if(s.equals("(")){
            for(int i=0;i<m;i++){
                c=p.indexOf(s,c)+1;
            }
            c=c-1;
        }
        if(s.equals(")")){
            c=p.length()-1;
            for(int i=0;i<m;i++){
                c=p.lastIndexOf(s, c)-1;
            }
            c=c+1;
        }

        return c;
    }

    /**
     * judge whether the structure of a sequence has two or more loops
     * @param String seq: sequence to be tested
     * @return boolean; if have two or more return true, or false
     */
    public static boolean hasTwoLoop(SimpleRNASequence rna){
        foldRNA(rna);
        int end5=rna.getStr().lastIndexOf("(");
        int start3=rna.getStr().indexOf(")");
        if(end5>=start3){
            return true;
        }
        return false;
    }

    /**
     * @return the window
     */
    public int getWindow() {
        return window;
    }

    /**
     * @param window the window to set
     */
    public void setWindow(int window) {
        this.window = window;
    }

    /**
     * @return the step
     */
    public int getStep() {
        return step;
    }

    /**
     * @param step the step to set
     */
    public void setStep(int step) {
        this.step = step;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the sequence
     */
    public SimpleSeq getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(SimpleSeq sequence) {
        this.sequence = sequence;
    }

    /**
     * @return the cutoff
     */
    public int getDistance() {
        return distance;
    }

    /**
     * @param cutoff the cutoff to set
     */
    public void setDistance(int distance) {
        this.distance = distance;
    }

    /**
     * @return the priList
     */
    public ArrayList<PriMiRNA> getPriList() {
        return primiRNAList;
    }

    /**
     * @param priList the priList to set
     */
    public void setPriList(ArrayList<PriMiRNA> priList) {
        this.primiRNAList = priList;
    }

    /**
     * @return the segList
     */
    public ArrayList<SimpleSeq> getSegList() {
        return fragmentList;
    }

    /**
     * @param segList the segList to set
     */
    public void setSegList(ArrayList<SimpleSeq> segList) {
        this.fragmentList = segList;
    }


}
