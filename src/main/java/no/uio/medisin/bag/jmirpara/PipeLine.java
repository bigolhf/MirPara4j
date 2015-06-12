
package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *
 * @author weibo
 */
public class PipeLine {
    
    static Logger logger = LogManager.getRootLogger();    

    
    private static final String[] knownModels = new String[] {"m", "p", "v", "o", "animal", "plant", "virus", "overall" };
    private static final String[] knownTests = new String[] {"svm", "hairpin", };

    private MiRParaConfiguation miRParaConfig;
    private String configFilename;

        
    private String inputFilename;
    private String outputFolder;
    
    
    private double version=3.0;
    private int window=500;
    private int step=250;
    private int start=1;
    private int distance=60;
    private double cutoff=0.8;
    private String model="overall";
    private int level=1;
    private File workingDir=new File(".");
    private String packageDir;
    private boolean append=false;
    private double progress;
//    private ArrayList<String> results=new ArrayList<String>();
    private ArrayList<String> results=new ArrayList<>();
    


    private ArrayList<SimpleSeq> seqList;
    private ArrayList<SimpleSeq> segList;
    private ArrayList<PriMiRNA> priList;
    private String[] last;
    private ArrayList<HashMap> featureList;
    private ArrayList<HashMap> trueList=new ArrayList<HashMap>();
    
    private String test;

    
    /**
     * class constructor
     * 
     */
    public PipeLine() {  
        miRParaConfig = new MiRParaConfiguation();
    }

    
    
    /**
     * predict the miRNAs in a set of input sequences
     * 
     * @throws IOException
     */
    public void oldRun() throws IOException{
        
        initialize();

        loadFastaSequenceData(new File(getInputFilename()));

        for(SimpleSeq seq : seqList){

            logger.info("Begin to extract possible pri-miRNAs from sequence "+seq.getId()+"...");
            findPrimiRNA(seq);
            seq.setSeq(null); //free space
            print("\n"+priList.size()+" putative pri-miRNAs are found\n");

            print("Begin to scan possible miRNAs from pri-miRNAs...");
            int i=0;
            trueList=new ArrayList<HashMap>();
            for(PriMiRNA pri : priList){
                findMiRNA(pri);
                priList.set(i, null); //free space
                i++;
                logger.info(Output.decimal(i*100.0/priList.size())+"%"+Output.backspace(Output.decimal(i*100.0/priList.size())+"%"));                
            }
            logger.info("\n"+trueList.size()+" miRNA candidates are found\n");


            reportResult(seq);

        }
        print("All have done!\nSee results at "+(workingDir.getCanonicalPath())+"\n");

    }
    
    /**
     * 
     * why is there a second run method?
     * which one should we use?
     * s.r. 9/3/2015
     * 
     * @throws IOException 
     */
    public void predictMiRNAsInInputSequences() throws IOException{
        initialize();
        ReadFastaFile rf = new ReadFastaFile(new File(getInputFilename()));
        while(rf.hasSeq()){
            SimpleSeq seq = rf.getOneSeq();  //each seq
            setOutfileName(workingDir, new File(getInputFilename()), seq.getId());
//            trueList=new ArrayList<HashMap>();
            trueList = new ArrayList<>();
            last = new String[0];
            serializeResults(seq);
            this.append = true;
            int totalNumOfMiRNA=0;
            logger.info("scanning for miRNAs in sequence " + seq.getId() + "...\n");
            int length = seq.getLength();
            int n=(length-window)/step+1; //num of window, num of step = n-1
            if(n<1) n=1;
            if(length-((n-1)*step+window)>19) n+=1;  // I think this must be to stop scanning fragments that are closer together than one miRNA
            int end=0, start=this.start-1;
            logger.info("will scan " + n + " fragments");
            for(int i=0;i<n;i++){
                
                if(start>=length) break;
                end = start+window;
                if(end>length) end = length;
                String id=seq.getName() + "_" + (start+1) + "-" + end;
                String subseq=seq.getSeq().substring(start,end);
                
                logger.info(".. scan region " + (start+1) + "-" + end + "...");
                SimpleSeq frag = new SimpleSeq(id,subseq);  //each frag
                frag.setStart(start+1); // we count from 1
                frag.setEnd(end); 
                frag.setName(seq.getId());
                
                StemLoopScanner sl = new StemLoopScanner();
                ArrayList<PriMiRNA> pris = sl.scanStemLoopFrom(frag);
                noRepeat(pris);
                logger.info(pris.size() + ".. pri-miRNA were found...");
                int before = trueList.size();
                for(PriMiRNA pri : pris){
                    findMiRNA(pri);                
                }
                int add = trueList.size() - before;
                if (add == 0)
                    logger.info(".. didn't find any miRNAs");
                else
                    if(add == 1)
                        logger.info(".. found 1 miRNA ...");
                    else
                        logger.info(".. found " + add + " miRNAs...");
                    
                //output the results in time to avoid memory leak
                if(trueList.size()>100){
                    totalNumOfMiRNA+=trueList.size();
                    serializeResults(seq);
                    trueList=new ArrayList<HashMap>();
                }
                progress = Double.parseDouble(Output.decimal((i+1)*100.0/n));
//                print(Output.decimal((i+1)*100.0/n)+"%"+Output.backspace(Output.decimal((i+1)*100.0/n)+"%"));
                print(Output.decimal((i+1)*100.0/n)+"%\n");
                start+= step;
            }
            logger.info("\n");
           
            serializeResults(seq);
            totalNumOfMiRNA+= trueList.size();
            

            if (totalNumOfMiRNA == 0)
                logger.info("didn't find any miRNAs");
            else
                if(totalNumOfMiRNA == 1)
                    logger.info("found a total of 1 miRNA candidate ...");
                else
                    logger.info("found a total of " + totalNumOfMiRNA + " miRNAs candidates...");
            
            this.append=false;

        }
        logger.info("Analysis complete. \n<Results are written to folder <" + (workingDir.getCanonicalPath()) + ">\n");
    }
    
    
    
    private void noRepeat(ArrayList<PriMiRNA> pris){
        //remove repeat pri
        Iterator it=pris.iterator();
        while(it.hasNext()){
            PriMiRNA pri=(PriMiRNA)(it.next());
            for(String id:last){
                if(pri.getId().equals(id)){
                    it.remove();
                    break;
                }
            }
        }
        
        //update last array;
        int n=pris.size();
        last=new String[n];
        int i=0;
        for(PriMiRNA pri:pris){
            last[i++]=pri.getId();
        }     
    }
    
    

    /**
     * processes the tab delimited output file (filename.tab) from an miRNA prediction run 
     * and summarizes the prediction results
     * 
     * Note: this might be a different file - still trying to work this out (sr: 11/06/2015)
     * 
     * line has format:
     *  pri_id pri_seq mi1_start, mi1_size, mi2_start, mi2_size
     * @throws IOException
     */
    public void testGivenMir() throws IOException{
        
        
        try{
            
            initialize();
            
            BufferedReader br=new BufferedReader(new FileReader(new File(getInputFilename())));
            String line;
            int total=0;
            int hit=0;
            
            while((line = br.readLine()) != null){
                if(line.equals("")) continue;
                //pri_id, pri_seq, mi1_start, mi1_size, mi2_start, mi2_size,...
                String[] entry=line.split("\t");

                PriMiRNA pri=new PriMiRNA(entry[0],entry[1]);
                if(StemLoopScanner.hasTwoLoop(pri)){
                   logger.info("W1: "+entry[0]+" is not a hairpin structure!");
                    continue;
                }
                MiGenesis parser=new MiGenesis(pri);
                parser.parsePrimiRNA();

                int end5=pri.getSeq5().length();
                int start3=end5+pri.getMidBase().length()+1;
                for(int i=2;i<entry.length;i+=2){
                    int miStart=Integer.parseInt(entry[i]);
                    int miSize=Integer.parseInt(entry[i+1]);
                    if(miStart<=start3 && miStart+miSize-1>=end5){
                            logger.info("W2: the given miRNA spans more than half the terminal loop of "+entry[0]+" , I cannot handle at present!");
                            continue;
                        }
                    parser.maturateMiRNA(miStart-1,miSize);
                    HashMap feat=parser.gatherFeatures();
                    //outputSV(feat);
                    SVMToolKit.judge(getModel(),feat,cutoff);
                    Boolean isMir=Boolean.parseBoolean(SVMToolKit.judgeResult());
                    if(isMir){
                        logger.info("Y: miRNA was detected at site "+miStart+" of "+entry[0]);
                        hit++;
                    }
                    else logger.info("N: miRNA was not detected at site "+miStart+" of "+entry[0]);
                    total++;
                }
            }
            br.close();
            System.out.println("cutoff="+cutoff);
            System.out.println("accuracy="+((double)hit/total));
        }
        catch(IOException ex){
            
        }
    }

    /**
     * read a SVM format matrix file and predict it
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void testProgram() throws FileNotFoundException, IOException{
        initialize();
        BufferedReader br=new BufferedReader(new FileReader(new File(getInputFilename())));
        BufferedWriter bw=new BufferedWriter(new FileWriter("test.out"));
        bw.write("given_label\tpreditc_label\tsvm_value\tprobability\n");
        String line;
        while((line = br.readLine()) != null){
            double[] results=SVMToolKit.predict_test(getModel(), line);
            StringBuilder out=new StringBuilder();
            for(int i=0;i<results.length;i++)
                out.append(results[i]).append("\t");
            out.replace(out.length()-1, out.length(), "\n");
            bw.write(out.toString());
        }
        bw.close();
        br.close();
    }


    /**
     * Load query sequences in FASTA format
     * 
     * @param fn
     * @throws IOException 
     */
    private void loadFastaSequenceData(File fn) throws IOException{
        
        seqList=new ReadFastaFile(fn).getSeqs();
    }

    
    
    /**
     * search for pri-miRNA in specified Simple Sequence
     * 
     * @param seq 
     */
    private void findPrimiRNA(SimpleSeq seq){
        
        StemLoopScanner sl = new StemLoopScanner(seq);
        
        sl.setWindow(window);
        sl.setStep(step);
        sl.setDistance(distance);
        sl.slidingWindow();
        sl.scanStemLoop();
        sl.noRepeat();
        
        priList = sl.getPriList();
        
        sl=null; // free space : sr: is this necessary?
    }

    
    /**
     * search for miRNA in specified pri-miRNA
     * 
     * @param priRNA
     * @throws IOException 
     */
    private void findMiRNA(PriMiRNA priRNA) throws IOException{ 
        MiGenesis mg=new MiGenesis(priRNA);
        try{
            mg.parsePrimiRNA();
        }catch(Exception e){
            System.out.print("Cannot parse the pri-miRNA ");
            System.out.println(priRNA.getId());
            System.out.println("sequence: "+priRNA.getSeq());
            System.out.println("structure:"+priRNA.getStr());
            return;
        }
        
        for(int w=20;w<25;w++){//miRNA size
            int i; //i refer to the start of miRNA, count from 0
           
            //5' strand
            for(i=0;i<priRNA.getTerminalLoopStart()-w;i++){//miRNA start point
                mg.maturateMiRNA(i,w);
                judgeIfMiRNA(mg.gatherFeatures());
            }
            //3' strand
            for(i=priRNA.getTerminalLoopEnd()+1;i<priRNA.getLength()-w;i++){//miRNA start point
                mg.maturateMiRNA(i,w);
                judgeIfMiRNA(mg.gatherFeatures());
            }

        }
        
        

    }



    private Boolean judgeIfMiRNA(HashMap feat) throws IOException{
        String result;        
        if(FeatureRange.rangeFilter(feat,getModel())){
            SVMToolKit.judge(getModel(),feat, cutoff);
            result=SVMToolKit.judgeResult();
            if(result.equals("TRUE"))
                trueList.add(feat);
        }
        else
           result="FALSE";
//        print(feat.get("miRNA_id")+"......"+result+Output.backspace(feat.get("miRNA_id")+"......"+result));
        return Boolean.parseBoolean(result);
    }
    
    
    
    /**
     * set output file name based on folder, input file and sequence name
     * 
     * @param dir
     * @param infile
     * @param seqname 
     */
    public void setOutfileName(File dir, File infile, String seqname){

        String basename = infile.getName().replaceAll("\\.\\w+", "");
//        if(dir.endsWith("/"))  dir=dir.substring(dir.length()-1);
        setOutputFolder(dir+"/"+basename+"_"+seqname);
        results.add(getOutputFolder());
        Output.fname=getOutputFolder();
    }

    
    
    /**
     * output
     * @param seq
     * @throws IOException
     */
    public void reportResult(SimpleSeq seq) throws IOException{
        Output.detail(trueList,seq,append);
        Output.overall(trueList,seq,append);
//        Output.distribution(trueList,seq);

    }
    
    
    
    /**
     * output prediction results for specified sequence
     * 
     * @param seq
     * @throws IOException 
     */
    public void serializeResults(SimpleSeq seq) throws IOException {
        
        Output.detail(trueList,seq,append);
        Output.overall(trueList,seq,append);
        
    }

    /**
     * output a support vector
     * @param feat
     */
    public void outputSV(HashMap feat){
        double[] feats=SVMToolKit.featValueOf(model,feat);
        StringBuilder pLine=new StringBuilder();
        for(int fs=1;fs<=feats.length;fs++)
                pLine.append("\t").append(fs).append(":").append(feats[fs-1]);
        System.out.println(pLine.toString());
    }

    /**
     * check we have everything we need before we start doing anything heavy
     * @throws IOException
     */
    private void initialize() throws IOException{
        //set the work directory
        //setWorkingDir(PathSet.getWorkingDir());
        //load RNAFold library
        //in JAR package
        miRParaConfig.ReadConfigurationFile(configFilename);
        
        logger.info("installation folder is " + miRParaConfig.getInstallationFolder());
        logger.info("model folder is " + Paths.get(miRParaConfig.getInstallationFolder(), miRParaConfig.getModelFolder()));

        packageDir=PathSet.getPackageDir();
        
//        PathSet.setLibDir(packageDir+"/lib");
        File mirbaseData=new File(miRParaConfig.getPathToModelData());
        
        //in NetBeans
//        packageDir="/home/wb/Desktop/program/jmi2_v13";
//        PathSet.setLibDir("/home/wb/Desktop/program/jmi2_v13/lib");
//        File mirbaseData=new File("/home/wb/Desktop/program/jmi2_v13/data/mirbase.dat");
        
        //load SVM model
        String modelName=packageDir+"/models/"+getModel()+"_"+getLevel()+".model";

        logger.info("loading model data file <" + miRParaConfig.getPathToModelData() + "");
        SVMToolKit.loadModel(miRParaConfig.getPathToModelData());

        logger.info("loading miRBase data file <" + miRParaConfig.getPathToModelData() + "");
        Output.loadMirBaseData(new File(miRParaConfig.getPathToMirbaseData()));     
        
    }

    
    
    
    /**
     * summarize run settings
     * 
     * @return String
     * 
     */
    public String reportRunSettings()
    {
        String summary = 
            StringUtils.repeat("*", 60) + "\n" 
         +  "filein      :\t" + getInputFilename() + "\n"
         +  "fileout     :\t" + getOutputFolder() + "\n\n"
         +  StringUtils.repeat("*", 60) + "\n\n" 
         +  miRParaConfig.reportParameters() + "\n\n"
         +  StringUtils.repeat("*", 60) + "\n\n";

        return summary;
    }
    
    
    
    /**
     * not sure what this is for
     * 
     * @param s 
     */     
    public void print(String s){
        System.out.print(s);
        System.out.flush();
    }

    
    /**
     * @return the window
     */
    public int getWindow() {
        return window;
    }

    
    /**
     * @param window the window to set
     * @throws IndexOutOfBoundsException if the window is less than zero
     */
    public void setWindow(int window) {
        if (window <= 0)
            throw new IndexOutOfBoundsException("window must be greater than 0");
        this.window = window;
    }

    
    /**
     * @return the step
     */
    public int getStep() {
        return step;
    }

    /**
     * @param step size
     * @throws IndexOutOfBoundsException if the step size is less than zero
     */
    public void setStep(int step) {
        if (step <= 0)
            throw new IndexOutOfBoundsException("step size must be greater than 0");
        this.step = step;
    }

    /**
     * @return the threshold
     */
    public int getDistance() {
        return distance;
    }

    /**
     * @param threshold the threshold to set
     * @throws IndexOutOfBoundsException if the distance is less than zero
     */
    public void setDistance(int distance) {
        if (distance <= 0)
            throw new IndexOutOfBoundsException("distance must be greater than 0");
        this.distance = distance;
    }

    /**
     * @return the model
     */
    public String getModel() {
        
        return model;
    }

    /**
     * @param model the model to set
     * @throws IllegalArgumentException if the model type is unknown
     */
    public void setModel(String model) {
        if(Arrays.asList(knownModels).contains(model.toLowerCase()) == false)
            throw new IllegalArgumentException("unrecognized model type");
        this.model = model;
        
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
        if (level < 1 || level > 20)
            throw new IndexOutOfBoundsException("level must be between 1 and 20");
        this.level = level;
    }

    /**
     * @return the workingDir
     */
    public File getWorkingDir() {
        return workingDir;
    }

    /**
     * @param workingDir the workingDir to set
     */
    public void setWorkingDir(String workingDir) {
        this.workingDir = new File(workingDir);
    }
    
    public void setWorkingDir(File workingDir){
        this.workingDir=workingDir;
    }

    /**
     * @return the cutoff
     */
    public double getCutoff() {
        return cutoff;
    }

    /**
     * @param cutoff the cutoff to set
     */
    public void setCutoff(double cutoff) {
        if (cutoff < 0.0 || cutoff > 1.0)
            throw new IndexOutOfBoundsException("cutoff must be between 0.0 and 1.0");
        this.cutoff = cutoff;
    }

    /**
     * @return the progress
     */
    public double getProgress() {
        return progress;
    }

    /**
     * @param progress the progress to set
     */
    public void setProgress(double progress) {
        this.progress = progress;
    }

    /**
     * @return the results
     */
    public ArrayList<String> getResults() {
        return results;
    }

    /**
     * @param results the results to set
     */
    public void setResults(ArrayList<String> results) {
        this.results = results;
    }

    /**
     * @return the test
     */
    public String getTest() {
        return test;
    }

    /**
     * @param test: set the test 
     */
    public void setTest(String test) {
        if(!Arrays.asList(knownTests).contains(test.toLowerCase())){
            this.test = test;
        }
        else
            logger.warn("unknown test " + test.toLowerCase());           
    }

    /**
     * @return the configFilename
     */
    public String getConfigFilename() {
        return configFilename;
    }

    /**
     * @param configFilename the configFilename to set
     */
    public void setConfigFilename(String configFilename) {
        this.configFilename = configFilename;
    }

    /**
     * @return the inputFilename
     */
    public String getInputFilename() {
        return inputFilename;
    }

    /**
     * @param inputFilename the inputFilename to set
     */
    public void setInputFilename(String inputFilename) {
        this.inputFilename = inputFilename;
    }

    /**
     * @return the outputFolder
     */
    public String getOutputFolder() {
        return outputFolder;
    }

    /**
     * @param outputFolder the outputFolder to set
     */
    public void setOutputFolder(String outputFolder) {
        this.outputFolder = outputFolder;
    }

}
