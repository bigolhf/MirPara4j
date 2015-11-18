
package no.uio.medisin.bag.jmirpara.train;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import no.uio.medisin.bag.jmirpara.MfeFoldRNA;
//import no.uio.medisin.bag.jmirpara.train.MySVM;
//import no.uio.medisin.bag.jmirpara.train.ParaSummary;
//import no.uio.medisin.bag.jmirpara.train.PriMiRNA;
import libsvm.svm_node;
import libsvm.svm_problem;
import no.uio.medisin.bag.jmirpara.PriMiRNA;

/**
 * this class prepare the training data of svm
 * @author weibo
 */
public class TrainingData {
    
    private ArrayList<TrainingEntry> inData;
    private PriMiRNA pri;
    /*
    private ParaSummary ps;
    private MySVM svxm=new MySVM();*/
    private int level=1;

    private svm_problem trainSet_p;
    private svm_problem trainSet_n;

//    private HashMap p_fas=new HashMap<String,String>();
//    private HashMap n_fas=new HashMap<String,String>();
    
    public TrainingData(){
        
    }

    public TrainingData(int level){
        this.level=level;
    }

    public void constructMatrix(ArrayList<MiRBaseEntry> mirb, int subnum, String taxo){
        ArrayList<TrainingEntry> pos=positiveData(mirb,subnum);
        trainSet_p=parsePara(pos, taxo);
        ArrayList<TrainingEntry> neg=negativeData(pos,level);
        trainSet_n=parsePara(neg, taxo);

        //report
        int total=0;
        for(TrainingEntry t:pos) total+=t.NumOfMi();
        for(TrainingEntry t:neg) total+=t.NumOfMi();
        System.out.println("Total Training Entries: "+total);
        System.out.println("Level: "+level);
    }
    
    /**
     * parse parameters of all entries
     * @param in
     */
    private svm_problem parsePara(ArrayList<TrainingEntry> in, String taxo){
        int num=0;
        for(TrainingEntry t:in) num+=t.NumOfMi();

        svm_node[][] data=new svm_node[num][];
        int i=0;
        /*
        for(TrainingEntry mi:in){
            pri=mi.getPriRNA();
            pri.setPriPara();
            for(int j=0;j<mi.NumOfMi();j++){
                pri.process(mi.getStart(j), mi.getSize(j), mi.getStrand(j));
                ps=new ParaSummary(pri, pri.getPreRNA(), pri.getMiRNA());
                data[i++]=svm.predictMatrix(mi.getLabel(j), ps.paraMatrix(taxo));
            }
        }
        return svm.trainMatrix(data);
        */
        return new svm_problem();
    }

    public ArrayList<TrainingEntry> positiveData(ArrayList<MiRBaseEntry> milist, int subnum){
        System.out.print("Reads miRBase data ");

        ArrayList<TrainingEntry> pos=new ArrayList<TrainingEntry>();

        int num=milist.size();
        ArrayList<Integer> rmis=randomization(num,num);//randomizate the miRBase entries

        MfeFoldRNA doFold;

        int evNum=0; //number of experimental verified miRNAs
        int uevNum=0; //number of unexperimental verified miRNAs

        int positive=0; //number of positive mirna
        int n=0;
/*        
        for(Integer m : rmis){

            if(subnum>0 && positive>=subnum) break; //get subnum mirans

            MiRBaseEntry mi=milist.get(m);

            TrainingEntry te=new TrainingEntry();

            pri=new PriMiRNA(mi.getPriAccession(),mi.getPriSequence());
            //fold pri-miRNA
            doFold = new MfeFoldRNA(pri.getSeq());
            doFold.cal();
            pri.setStr(doFold.getStructure());
            pri.setEnergy(doFold.getEnergy());

            te.setPriRNA(pri);//store pri

            int end5=pri.getStr().lastIndexOf("(")+1;//5' end position of hairpin,count from 1
            int start3=pri.getStr().indexOf(")")+1;//3' start position of hairpin,count from 1

            for(int i=0;i<mi.getMatNum();i++){
                //the miRNA should be experimental verified miRNA
                if(mi.getEvidence(i).equals("experimental")){
                    evNum+=1;

                    //the miRNA should not have two or more loops
                    if(end5>=start3) continue;

                    int strand;
                    if(mi.getMatEnd(i)<=end5) strand=5; //lie in 5'
                    else if(mi.getMatStart(i)>=start3) strand=3; //lie in 3'
                    else continue; //the miRNA should not lie in loop area

                    te.addMiEntry(mi.getMatStart(i)-1, mi.getMatEnd(i)-mi.getMatStart(i)+1, strand, 1);
                    //for the fasta format mirna data
//                    te.setMiId(mi.getMatID(i));
//                    p_fas.put(te.getMiId(te.NumOfMi()-1), te.getMiSeq(te.NumOfMi()-1));

                    positive++;
                }
                else uevNum+=1;


            }
            //add positive data
            if(te.NumOfMi()>0)
                pos.add(te);
            
            n++;
            System.out.print(n+":"+mi.getPriID()+backspace(n+":"+mi.getPriID()));
        }
*/        
        System.out.println();
        //report
        System.out.println("Loads "+(evNum+uevNum)+" Entries");
        System.out.println("Experimental verified Entries: "+evNum);
        System.out.println("Non-experimental verified Entries: "+uevNum);
        System.out.println("Positive Entries: "+positive);
        return pos;
    }

    public ArrayList<TrainingEntry> negativeData(ArrayList<TrainingEntry> pos, int level){
        Random r=new Random();//random number producer
        int rS=0;
        ArrayList<TrainingEntry> nes=new ArrayList<TrainingEntry>();
        TrainingEntry ne;
        int strand=0;
        int negative=0;
        for(TrainingEntry po : pos){
            ne=new TrainingEntry();
            ne.setPriRNA(po.getPriRNA());
            int end5=po.getPriRNA().getStructureStr().lastIndexOf("(");//from 0
            int start3=po.getPriRNA().getStructureStr().indexOf(")");//from 0
            int endP5=end5+(start3-end5-1)/2;//5' endpoint of hairpin,count from 0
            int endP3=start3-(start3-end5-1)/2;//3' endpoint of hairpin, count from 0

            int num=po.NumOfMi(); //the number of postive mi
            int size=po.getPriRNA().getLength(); //the size of pri in the entry
            for(int n=0;n<level;n++){ //each level
                for(int i=0;i<num;i++){ //each miRNA in the entry
                    int miL=po.getSize(i); //the positive mi size
                    HashMap bases=new HashMap();
                    int flag=0;
                    while(bases.size()<size){ //the negative mi start
                        rS=r.nextInt(size); //get a random start
                        if(bases.containsKey(rS)) continue;
                        bases.put(rS, null);//store the positions which have been used
                        //the random start is at least 5-bp from the true start positions on the entry
                        int j;
                        for(j=0;j<num;j++){
                            if(Math.abs(po.getStart(i)-rS)<5)
                                break;
                        }
                        if(j<num) continue;

                        //if the random start at 5' strand
                        if((rS<endP5 && rS+miL-1<=endP5)){
                            strand=5;
                            flag=1;
                            break;
                        }
                        //if the random start at 3' strand
                        else if(rS>=endP3 && rS+miL<size){
                            strand=3;
                            flag=1;
                            break;
                        }
                        //else continue;
                    }
                    if(flag==1){
                        ne.addMiEntry(rS,miL,strand,0);
                        //for the fasta format mirna data
//                        ne.setMiId(po.getMiId(i)+"n"+n);
//                        n_fas.put(ne.getMiId(i), ne.getMiSeq(i));

                        negative++;
                    }
                    else
                        System.out.println("no negative data found on "+po.getPriRNA().getId());
                }
            }
            nes.add(ne);
        }
        System.out.println("Negative Entries: "+negative);
        return nes;
    }

    /**
     * generate m unrepeatable numbers between 0 and n
     * @param n
     * @param m
     * @return
     */
    public ArrayList<Integer> randomization(int n, int m){
        ArrayList<Integer> list=new ArrayList<Integer>();
        Random rand=new Random();
        boolean[] bo=new boolean[n];
        int num=0;
        for(int i=0;i<m;i++){
            do{
                num=rand.nextInt(n);
            }while(bo[num]);
            bo[num]=true;
            list.add(num);
        }
//        for(int i=0;i<m;i++)
//            list.add(i);
        
        return list;
    }

    

    /**
     * save the training matrix to file
     * @param filename
     * @throws IOException
     */
    public void saveTrainingData(String filename) throws IOException{
        BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
        BufferedWriter bwp=new BufferedWriter(new FileWriter(filename+".positive"));
        BufferedWriter bwn=new BufferedWriter(new FileWriter(filename+".negative"));
        String line;
        //positive
        for(int i=0;i<trainSet_p.l;i++){
            line=(trainSet_p.y)[i]+"\t";
            for(int j=0;j<trainSet_p.x[0].length;j++){
                line+=(trainSet_p.x)[i][j].index+":";
                line+=(trainSet_p.x)[i][j].value+"\t";
            }
            line=line.trim();

            bw.write(line.trim());
            bw.write("\n");
            bwp.write(line.trim());
            bwp.write("\n");
        }
        bwp.close();
        //negative
        for(int i=0;i<trainSet_n.l;i++){
            line=(trainSet_n.y)[i]+"\t";
            for(int j=0;j<trainSet_n.x[0].length;j++){
                line+=(trainSet_n.x)[i][j].index+":";
                line+=(trainSet_n.x)[i][j].value+"\t";
            }
            line=line.trim();

            bw.write(line.trim());
            bw.write("\n");
            bwn.write(line.trim());
            bwn.write("\n");
        }
        bwn.close();

        bw.close();
    }


//    public void printMi(BufferedWriter bw, HashMap mil) throws IOException{
//
//
//        Iterator keys=mil.keySet().iterator();
//        Iterator values=mil.values().iterator();
//        while(keys.hasNext()){
//            String key=keys.next().toString();
//            String value=values.next().toString();
//            bw.write(">"+key+"\n"+value+"\n");
//        }
//
//    }
//    public void saveFasta(String filename) throws IOException{
//        BufferedWriter bwp=new BufferedWriter(new FileWriter(filename+".positive.fasta"));
//        printMi(bwp,p_fas);
//        bwp.close();
//        BufferedWriter bwn=new BufferedWriter(new FileWriter(filename+".negative.fasta"));
//        printMi(bwn,n_fas);
//        bwp.close();
//    }

    public String backspace(String con){
        return repeatStr("\b",con.length());
    }

    public String repeatStr(String str, int num){
        String strs="";
        for(int i=0;i<num;i++)
            strs+=str;
        return strs;
    }



    /**
     * @return the trainSet
     */
    public svm_problem getPositiveTrainSet() {
        return trainSet_p;
    }
    public svm_problem getNegativeTrainSet() {
        return trainSet_n;
    }

    /**
     * @param trainSet the trainSet to set
     */
    public void setPositiveTrainSet(svm_problem trainSet) {
        this.trainSet_p = trainSet;
    }
    public void setNegativeTrainSet(svm_problem trainSet) {
        this.trainSet_n = trainSet;
    }

    /**
     * @return the level
     */
    public int getLevel() {
        return level;
    }

    /**
     * @param level the level to set
     */
    public void setLevel(int level) {
        this.level = level;
    }

    
}
