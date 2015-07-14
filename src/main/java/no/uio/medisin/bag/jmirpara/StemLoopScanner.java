
package no.uio.medisin.bag.jmirpara;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static no.uio.medisin.bag.jmirpara.MiRNAPredictionPipeLine.logger;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * this class finds hairpin structure from an input string in Vienna Bracket Notation.
 * 
 * @see <a href=http://rna.tbi.univie.ac.at/help.html#A6>  http://rna.tbi.univie.ac.at/help.html</a>
 * 
 * 
 * the order of calling methods after instantiated:
 * slidingWindow();
 * scanStemLoop();
 * noRepeat();
 * @author weibo & simon rayner
 */
public class StemLoopScanner {
    static Logger logger = LogManager.getRootLogger();    

    private int window=500;
    private int step=250;
    private int start=1;
    private int minHairpinLength = 60;

    private SimpleSeq               simpleSeq;
    
    private ArrayList<SimpleSeq>    fragmentList;
    private ArrayList<PriMiRNA>     primiRNAList;

    /*
      private static String regSL="([\\(\\.]*\\()(\\.+)(\\)[\\)\\.]*)";
      
      This can be broken down into
     
      {'('} + {zero or more '.' or '('} + {stem loop} + {zero or more '.' or ')'} + {')'}
    
      see 
    
    
    */
//  
    
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
        simpleSeq = seq;
    }
    
    
    
    /**
     * Class constructor from subsequence of input sequence
     * Allows the user to instruct the class to break the sequence
     * into smaller fragments Useful for long sequences which are 
     * likely to return avoid complex full length structures
     * 
     * 
     * @param seq
     * @param window
     * @param step
     * @param distance 
     */
    public StemLoopScanner(SimpleSeq seq, int window, int step, int distance){
        
        this.simpleSeq   = seq;
        this.window     = window;
        this.step       = step;
        this.minHairpinLength   = distance;
        
    }
    

    /**
     * 
     * slide window with given step size to generate a serial of segments for testing.
     * before calling, one may set the sliding window size, step increment and start position
     * 
     * It seems this should be removed as this is already implemented in the Pipeline class
     * @author sr 
     */
    public void breakQuerySeqIntoOverlappingFragments(){
        
        fragmentList = new ArrayList<SimpleSeq>();
        int length=simpleSeq.getLength();
        int n=(length-window)/step;
        if(n<0) n=0;
        if(length-n*step+window>19) n+=1;
        int end=0;start=0;
        for(int i=0;i<=n;i++){
            
            if(start>=length) break;
            end = start + window;
            if(end>length) end=length;
            String id = simpleSeq.getName() + "_" + (start+1) + "-"+end;
            String subseq = simpleSeq.getSeq().substring(start,end);

            SimpleSeq frag = new SimpleSeq(id,subseq);
            frag.setStart(start+1);// count from 1
            frag.setEnd(end); //count from 1
            frag.setName(simpleSeq.getId());
            fragmentList.add(frag);
            start+= step;
            
        }
        
    }
    

    /**
     * 
     * generate a list of candidate pri-miRNAs from all the fragments
     * 
     */
    public ArrayList<PriMiRNA> foldAndScanForStemloopByFragments(){
        
        fragmentList = new ArrayList<SimpleSeq>();

        int length=simpleSeq.getLength();
        int n=(length-window)/step+1;
        if(n<1) n=1;
        if(length-((n-1)*step+window)>19) n+=1;  // I think this must be to stop scanning fragments that are closer together than one miRNA

        
        int end=0; start=0;
        for(int i=0;i<n;i++){
            
            if(start>=length) break;
            end = start + window;
            if(end>length) end=length;
            String id = simpleSeq.getName() + "_" + (start+1) + "-"+end;
            String subseq = simpleSeq.getSeq().substring(start,end);

            SimpleSeq frag = new SimpleSeq(id,subseq);
            frag.setStart(start+1);// count from 1
            frag.setEnd(end); //count from 1
            frag.setName(simpleSeq.getId());
            fragmentList.add(frag);
            start+= step;
            
        }
        
        primiRNAList = new ArrayList<PriMiRNA>();
        
        int i=1;
        for(SimpleSeq frag : fragmentList){
            logger.info("SL: scan region " + (frag.getStart()) + "-" + frag.getEnd() + "...");
            ArrayList<PriMiRNA> fragPriMiRNAList = foldAndScanSequenceForStemloop(frag);
            String se= "AAGCUGGCAUUCUAUAUAAGAGAGAAACUACACGCAGCGCCUCAUUUUGUGGGUCA"
              + "CCAUAUUCUUGGGAACAAGAGCUACAGCAUGGGGCAAAUCUUUCUGUUCCCAAUCCUCUGGGA"
              + "UUCUUUCCCGAUCACCAGUUGGACCCUGCGUUUGGAGCCAACUCAAACAAUCCAGAUUGGGAC"
              + "UUCAACCCCAACAAGGAUCACUGGCCAGAGGCAAAUCAGGUAGGAGCGGGAGCAUUCGGGCCA"
              + "GGGUUCACCCC";
            logger.info("SLEE:" + MfeFoldRNA.foldSequence(se));
            logger.info("   " + fragPriMiRNAList.size() + " pri-miRNA were found");
            primiRNAList.addAll(fragPriMiRNAList);
//            System.out.print(OutputMiRNAPredictions.decimal(i*100.0/fragmentList.size())+"%"+OutputMiRNAPredictions.backspace(OutputMiRNAPredictions.decimal(i*100.0/fragmentList.size())+"%"));
            i++;
        }

        return primiRNAList;
    }
    
    /**
     * Fold query sequence and scan resulting structure for hairpin structure
     * @param seq
     * @return 
     */
    public ArrayList<PriMiRNA> foldAndScanSequenceForStemloop(SimpleSeq seq){
        SimpleRNASequence rna = new SimpleRNASequence(seq);

        foldRNASequence(rna);

        return extractHairpinsFromBracketNotation(rna);
    }

    
    
    
    /**
     * 
     * fold an RNA sequence using RNAFold and return predicted structure
     * in Vienna Bracket Notation
     * 
     * @param rnaSequence : the sequence to be folded
     */
    public static void foldRNASequence(SimpleRNASequence rnaSequence){
                
        rnaSequence.setEnergy(MfeFoldRNA.foldSequence(rnaSequence.getSeq()));
        rnaSequence.setStructureString(MfeFoldRNA.getStructure());
        logger.info("Energy" + rnaSequence.getEnergy());

    }

    
    
    
    /**
     * 
     * extract hairpins from a rna secondary structure defined using the 
     * Vienna Bracket Notation
     * 
     * @param rnaStructure - SimpleRNASequence (with predicted secondary structure) 
     * 
     * @return ArrayList<@link{PriMiRNA}> list of PrimiRNA sequences
     * 
     */
    public ArrayList<PriMiRNA> extractHairpinsFromBracketNotation(SimpleRNASequence rnaStructure){
        
        ArrayList<PriMiRNA> primiRNAList=new ArrayList<PriMiRNA>();
        
        String str = rnaStructure.getStructureStr(); 
        m = pattern.matcher(str);

        while(m.find()){

            // strip off the hanging ends from the structure
            // ..(((..((((...)))).)))... 
            //          becomes
            //   (((..((((...)))).)))   
            int start1      = str.lastIndexOf(")", m.start(1)) + 1; //replace m.start(1)
            String left     = str.substring(start1, m.end(1));//the 5' side of the stemloop
            String right    = m.group(3);//the 3' side of the stemloop
            
            int n = b2Num(left)+b2Num(right);//compare the two sides of the stem
            int slStart=0, slEnd=0, l=0, r=0;
            
            //if str is like (..((...)).. return ..((...))..
            if(n>0){
                l=bIndex(left,"(",n)+1; //new start of left
                left=left.substring(l); //new left
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
            
            String subId    = rnaStructure.getName()+"_mir_"; //the id of the stemloop
            String subSeq   = rnaStructure.getSeq().substring(slStart, slEnd); //seq of the stemloop
            String subStr   = left + m.group(2) + right; //structure of the stemloop


            if(subStr.length() < minHairpinLength)
                continue;

            //create a new pri-miRNA
            PriMiRNA pri = new PriMiRNA(subId, subSeq);

            //if the stemloop has two or more loops
            if(hasMultipleLoops(pri)) continue;
            if(pattern.matcher(pri.getStructureStr()).matches() == false) continue;
            
            pri.setName(rnaStructure.getName());
            pri.setStart(slStart+rnaStructure.getStart());
            pri.setEnd(slEnd-1+rnaStructure.getStart());
            pri.setId(pri.getId()+pri.getStart()+"-"+pri.getEnd());
            
            primiRNAList.add(pri);

            
        }
        return primiRNAList;
    }

    
    
    
    /**
     * 
     * remove repeat hairpins caused by overlapping fragments
     * 
     */
    public void removeDuplicatePrimiRNAs(){
        HashMap map=new HashMap();
        for(PriMiRNA pri:primiRNAList)
            map.put(pri.getId(), pri);

        primiRNAList = new ArrayList<PriMiRNA>(map.values());
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
     * check if the structure of a sequence has two or more loops
     * 
     * @param  rna SimpleRNASequence seq: sequence to be tested
     * @return boolean; if have two or more return true, or false
     * 
     */
    public static boolean hasMultipleLoops(SimpleRNASequence rna){
        
        foldRNASequence(rna);
        int end5=rna.getStructureStr().lastIndexOf("(");
        int start3=rna.getStructureStr().indexOf(")");
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
        return simpleSeq;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(SimpleSeq sequence) {
        this.simpleSeq = sequence;
    }

    /**
     * @return the cutoff
     */
    public int getDistance() {
        return minHairpinLength;
    }

    /**
     * @param cutoff the cutoff to set
     */
    public void setDistance(int distance) {
        this.minHairpinLength = distance;
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
