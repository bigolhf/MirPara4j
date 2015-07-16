
package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.yaml.snakeyaml.Yaml;


/**
 *
 * @author weibo
 */
public class MiRNAPredictionPipeLine {
    
    static Logger logger = LogManager.getRootLogger();    
    

    private static final int        MIN_MIRNA_LEN = 20;
    private static final int        MAX_MIRNA_LEN = 25;
    
    private static final String[]   knownModels = new String[] {"m", "p", "v", "o", "animal", "plant", "virus", "overall" };
    private static final String[]   knownTests = new String[] {"svm", "hairpin", };

    private MiRParaConfiguation     miRParaConfig;
    private String                  configFilename;

        
    private String                  inputFilename;
    private String                  outputFolder;
    private String                  outfilePrefix;

    static  Yaml                    yaml = new Yaml();
    static  String                  FileSeparator = System.getProperty("file.separator");
    
    private String                  configurationFile="";
    
    private String                  installationFolder="";
    private String                  modelFolder="";
    private String                  rnaFoldFolder="";
    private String                  dataFolder="";
    
    private String                  rnaFoldDll="";
    private String                  mirbaseDataFile="";
    
    private String                  pathToLibrary = "";
    private String                  pathToMirbaseData = "";
    
    private String                  pathToModelData = "";
    
    private int                     maxNumOfPredictionListEntries = 100;
    //
    // to add:
    //
    //  max size of prediction list
    //  min miRNA length
    //  max miRNA length
    //
    
    private int                     window      = 500;
    private int                     step        = 250;
    private int                     start       = 1;
    private int                     distance    = 60;
    private double                  cutoff      = 0.8;
    private String                  model       = "overall";
    private int                     level       = 1;
    private File                    workingDir  = new File(".");
    private String                  packageDir;
    private boolean                 append      = false;
    private double                  progress;
//    private ArrayList<String> results=new ArrayList<String>();
    private ArrayList<String>       results     = new ArrayList<>(); // <- is this necessary ?
    


    private ArrayList<SimpleSeq>    seqList;
    private ArrayList<SimpleSeq>    segList;
    private ArrayList<PriMiRNA>     priList;
    private String[]                last;
    private ArrayList<HashMap>      featureList;
    private ArrayList<HashMap>      predictionList=new ArrayList<HashMap>();
    
    private String test;

    public MiRNAPredictionPipeLine() {  
        miRParaConfig = new MiRParaConfiguation();
    }

    
    
    /**
     * 
     * This pipeline performs miRNA prediction on an input list of FASTA query sequences.
     * 
     * long sequences can be split into smaller overlapping fragments to
     * reduce run time and prevent prediction of long range secondary structures
     * which are irrelevant here since we assume miRNAs are located in short
     * hairpin structures
     * 
     * run parameters can be specified in the YAML file that is required as input
     * 
     * @throws IOException 
     * 
     */
    public void predictMiRNAsInQuerySequences() throws IOException{
        
        this.initializePipeline();
        
        ReadFastaFile rf = new ReadFastaFile(new File(getInputFilename()));
        
        while(rf.hasSeq()){
            
            SimpleSeq querySeq = rf.getOneSeq();  //each seq
            setOutfileName(workingDir, new File(getInputFilename()), querySeq.getId());            

            predictionList = new ArrayList<>();
            last = new String[0];
            serializeResults(querySeq);
            
            this.append = true;
            int totalNumOfMiRNA=0;
            
            logger.info("analyzing sequence " + querySeq.getId() + "...\n");
            int length = querySeq.getLength();
            int n=(length-window)/step+1; //num of window, num of step = n-1
            if(n<1) n=1;
            if(length -((n-1)*step + window)>19) n+=1;  // I think this must be to stop scanning fragments that are closer together than one miRNA
            int end=0, start=this.start-1;
            
            logger.info("-- sequence is " + length + " nt");
            logger.info("-- will break into " + n + " fragments");
            
            StemLoopScanner sl = new StemLoopScanner(querySeq, window, step, distance);
            ArrayList<PriMiRNA> priMiRNAsInFrag = sl.foldAndScanForStemloopByFragments();
            
            removeDuplicatePrimiRNAs(priMiRNAsInFrag);
            logger.info(priMiRNAsInFrag.size() + " pri-miRNAs were found in the query sequence...");
            
            int before = predictionList.size();
            logger.info(" ");
            logger.info("scanning pri-miRNAs for miRNAs...");
            for(PriMiRNA primiRNA : priMiRNAsInFrag){
                logger.info(".." + primiRNA.getId());
                findMiRNAsInPrimiRNA(primiRNA);                
            }
                
            if (totalNumOfMiRNA == 0)
                logger.info("didn't find any miRNAs");
            else
                if(totalNumOfMiRNA == 1)
                    logger.info("found a total of 1 miRNA candidate ...");
                else
                    logger.info("found a total of " + totalNumOfMiRNA + " miRNAs candidates...");
            
            //output the results now to avoid memory leak
            if(predictionList.size()>100){
                logger.info(" ");
                logger.info("flush prediction list");
                logger.info(" ");
                totalNumOfMiRNA += predictionList.size();
                serializeResults(querySeq);
                predictionList = new ArrayList<HashMap>();
            }
//                progress = Double.parseDouble(OutputMiRNAPredictions.decimal((i+1)*100.0/n));
//                print(Output.decimal((i+1)*100.0/n)+"%"+Output.backspace(Output.decimal((i+1)*100.0/n)+"%"));
//                print(OutputMiRNAPredictions.decimal((i+1)*100.0/n)+"%\n");
                start+= step;
            logger.info("\n");
           
            serializeResults(querySeq);
            totalNumOfMiRNA+= predictionList.size();
            

            
            this.append=false;

        }
        logger.info("Analysis complete. \n<Results are written to folder <" + (workingDir.getCanonicalPath()) + ">\n");
    }
    
    
    
    private void removeDuplicatePrimiRNAs(ArrayList<PriMiRNA> primiRNAList){

        Iterator itPriList = primiRNAList.iterator();
        while(itPriList.hasNext()){
            PriMiRNA pri=(PriMiRNA)(itPriList.next());
            for(String id:last){
                if(pri.getId().equals(id)){
                    itPriList.remove();
                    break;
                }
            }
        }
        
        //update last array;
        int n = primiRNAList.size();
        last = new String[n];
        int i=0;
        for(PriMiRNA pri:primiRNAList){
            last[i++] = pri.getId();
        }     
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
        /*
        StemLoopScanner sl = new StemLoopScanner(seq);
        
        sl.setWindow(window);
        sl.setStep(step);
        sl.setDistance(distance);
        sl.breakQuerySeqIntoOverlappingFragments();
        sl.foldAndScanByFragments();
        sl.noRepeat();
        
        priList = sl.getPriList();
        
        sl=null; // free space : sr: is this necessary?
        */
    }

    
    
    /**
     * search for miRNA in specified pri-miRNA
     * 
     * @param priRNA
     * @throws IOException 
     */
    private void findMiRNAsInPrimiRNA(PriMiRNA priRNA) throws IOException{ 
        CharacterizedPriMiRNA characterizedPrimiRNA = new CharacterizedPriMiRNA(priRNA);
        
        try{
            
            characterizedPrimiRNA.parsePrimiRNA();
            
        }catch(Exception e){
            
            System.out.print("Couldn't parse the pri-miRNA ");
            System.out.println(priRNA.getId());
            System.out.println("sequence: "+priRNA.getSeq());
            System.out.println("structure:"+priRNA.getStructureStr());
            return;
            
        }
        
        /*
          cycle through all possible miRNAs in the hairpin loop,
          characterize, and use as input to SVM for classification
        */
        for(int mirLen=MIN_MIRNA_LEN;mirLen<MAX_MIRNA_LEN;mirLen++){//miRNA size
            int miRStart; // count from 0
           
            //5' strand
            for(miRStart=0;miRStart<priRNA.getTerminalLoopStart()-mirLen;miRStart++){  //miRNA start point
                characterizedPrimiRNA.defineAndCharacterizePreMiPair(miRStart, mirLen);
                classifyPriPreMiRTriplet(characterizedPrimiRNA.characterizePriPreMiTriplet());
            }
            
            //3' strand
            for(miRStart=priRNA.getTerminalLoopEnd()+1;miRStart<priRNA.getLength()-mirLen;miRStart++){//miRNA start point
                characterizedPrimiRNA.defineAndCharacterizePreMiPair(miRStart,mirLen);
                classifyPriPreMiRTriplet(characterizedPrimiRNA.characterizePriPreMiTriplet());
            }

        }
        
        

    }


    
    /**
     * classify the input sequence using the specified model
     * 
     * @param HashMap of features characterizing the pri-miRNA, pre-miRNA, miRNA triplet
     * 
     * @return boolean of classification outcome
     * 
     * @throws IOException 
     */
    private Boolean classifyPriPreMiRTriplet(HashMap tripletFeatures) throws IOException{
        
        String classficationResult;  
        
        if(FeatureRange.featureInRange(tripletFeatures, getModel())){
            SVMToolKit.judge(getModel(),tripletFeatures, cutoff);
            classficationResult = SVMToolKit.judgeResult();
            if(classficationResult.equals("TRUE"))
                predictionList.add(tripletFeatures);
        }
        else
           classficationResult="FALSE";

        return Boolean.parseBoolean(classficationResult);
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

        //setOutputFolder(dir+"/"+basename+"_"+seqname);
        results.add(getOutputFolder());
        OutputMiRNAPredictions.outputFilePrefix = getOutputFolder();
    }

    
        
    public void SetOutfilePrefix(String seqName){
        OutputMiRNAPredictions.outputFilePrefix = this.outputFolder + FileSeparator + seqName;
    }
    
    /**
     * output prediction results for specified sequence
     * 
     * @param seq
     * @throws IOException 
     */
    public void serializeResults(SimpleSeq seq) throws IOException {
        
        OutputMiRNAPredictions.serializePredictionDetails(predictionList, seq, append);
        OutputMiRNAPredictions.serializePredictionSummary(predictionList, seq, append);
        OutputMiRNAPredictions.serializePredictionsAsHTML(predictionList, seq);
        
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
    private void initializePipeline() throws IOException{

        logger.info("reading configuration file");
        this.readConfigurationFile(configFilename);
        logger.info(this.reportConfiguration() + "\n");
        logger.info(this.reportParameters() + "\n");
        logger.info(this.reportRunSettings() + "\n");
        
        
        logger.info("loading model data file <" + this.getPathToModelData() + "");
        SVMToolKit.loadModel(this.getPathToModelData());

        logger.info("loading miRBase data file <" + this.getPathToModelData() + "");
        OutputMiRNAPredictions.loadMirBaseData(new File(this.getPathToMirbaseData()));     
        
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
    
    
    
        /**
     * 
     * 
     * @throws IOException 
     */
    public void predictMiRNAsInQuerySequencesWithSplit() throws IOException{
  
        
        initializePipeline();
        
        ReadFastaFile rf = new ReadFastaFile(new File(getInputFilename()));
        while(rf.hasSeq()){
            
            SimpleSeq querySeq = rf.getOneSeq();  //each seq
            
            setOutfileName(workingDir, new File(getInputFilename()), querySeq.getId());
            this.SetOutfilePrefix(querySeq.getId());

            predictionList = new ArrayList<>();
            last = new String[0];
            serializeResults(querySeq);
            
            this.append = true;
            int totalNumOfMiRNA=0;
            logger.info("scanning for miRNAs in sequence " + querySeq.getId() + "...\n");
            int length = querySeq.getLength();
            int n=(length-window)/step+1; //num of window, num of step = n-1
            if(n<1) n=1;
            if(length -((n-1)*step + window)>19) n+=1;  // I think this must be to stop scanning fragments that are closer together than one miRNA
            int end=0, start=this.start-1;
            
            logger.info("-- sequence is " + length + " nt");
            logger.info("-- will break into " + n + " fragments");
            for(int i=0; i<n; i++){
                
                if(start>=length) break;
                end = start + window;
                if(end>length) end = length;
                String id=querySeq.getName() + "_" + (start+1) + "-" + end;
                String subseq=querySeq.getSeq().substring(start,end);
                
                logger.info("PL: scan region " + (start+1) + "-" + end + "...");

                SimpleSeq frag = new SimpleSeq(id,subseq);  //each frag
                frag.setAbsStartInQuerySeq(start+1); // we count from 1
                frag.setAbsEndInQuerySeq(end); 
                frag.setName(querySeq.getId());
                StemLoopScanner sl = new StemLoopScanner();

                ArrayList<PriMiRNA> pris = sl.foldAndScanSequenceForStemloop(frag);
                removeDuplicatePrimiRNAs(pris);

                logger.info(pris.size() + ".. pri-miRNA were found...");
                int before = predictionList.size();
                for(PriMiRNA pri : pris){
                    findMiRNAsInPrimiRNA(pri);                
                }
                
                int add = predictionList.size() - before;
                if (add == 0)
                    logger.info(".. didn't find any miRNAs");
                else
                    if(add == 1)
                        logger.info(".. found 1 miRNA ...");
                    else
                        logger.info(".. found " + add + " miRNAs...");
                    
                //output the results in time to avoid memory leak
                if(predictionList.size()> getMaxNumOfPredictionListEntries()){
                    totalNumOfMiRNA+=predictionList.size();
                    serializeResults(querySeq);
                    predictionList=new ArrayList<HashMap>();
                }
                progress = Double.parseDouble(OutputMiRNAPredictions.decimal((i+1)*100.0/n));
//                print(Output.decimal((i+1)*100.0/n)+"%"+Output.backspace(Output.decimal((i+1)*100.0/n)+"%"));
                print(OutputMiRNAPredictions.decimal((i+1)*100.0/n)+"%\n");
                start+= step;
            }
            logger.info("\n");
           
            serializeResults(querySeq);
            totalNumOfMiRNA+= predictionList.size();
            

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
    
    

    
    
    /**
     * predict the miRNAs in a set of input sequences
     * 
     * @throws IOException
     */
    public void oldRun() throws IOException{
        
        initializePipeline();

        loadFastaSequenceData(new File(getInputFilename()));

        for(SimpleSeq seq : seqList){

            logger.info("Begin to extract possible pri-miRNAs from sequence "+seq.getId()+"...");
            findPrimiRNA(seq);
            seq.setSeq(null); //free space
            print("\n"+priList.size()+" putative pri-miRNAs are found\n");

            print("Begin to scan possible miRNAs from pri-miRNAs...");
            int i=0;
            predictionList=new ArrayList<HashMap>();
            for(PriMiRNA pri : priList){
                findMiRNAsInPrimiRNA(pri);
                priList.set(i, null); //free space
                i++;
                logger.info(OutputMiRNAPredictions.decimal(i*100.0/priList.size())+"%"+OutputMiRNAPredictions.backspace(OutputMiRNAPredictions.decimal(i*100.0/priList.size())+"%"));                
            }
            logger.info("\n"+predictionList.size()+" miRNA candidates are found\n");


            serializeResults(seq);

        }
        print("All have done!\nSee results at "+(workingDir.getCanonicalPath())+"\n");

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
            
            initializePipeline();
            
            BufferedReader br=new BufferedReader(new FileReader(new File(getInputFilename())));
            String line;
            int total=0;
            int hit=0;
            
            while((line = br.readLine()) != null){
                if(line.equals("")) continue;
                //pri_id, pri_seq, mi1_start, mi1_size, mi2_start, mi2_size,...
                String[] entry=line.split("\t");

                PriMiRNA pri=new PriMiRNA(entry[0],entry[1]);
                if(StemLoopScanner.hasMultipleLoops(pri)){
                   logger.info("W1: "+entry[0]+" is not a hairpin structure!");
                    continue;
                }
                CharacterizedPriMiRNA parser=new CharacterizedPriMiRNA(pri);
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
                    parser.defineAndCharacterizePreMiPair(miStart-1,miSize);
                    HashMap feat=parser.characterizePriPreMiTriplet();
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
        initializePipeline();
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
     * load run parameters and general program settings from configuration file
     * in YAML format
     * 
     * @param configFile : absolute path to configuration file
     * 
     */
    public void readConfigurationFile(String configFile){
        try{            
            setConfigurationFile(configFile);
            HashMap mirparaOptions = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));
            
            HashMap modelsOptions = (HashMap) mirparaOptions.get("model");
            this.setModelFolder((String) modelsOptions.get("folder"));

            HashMap pathsOptions = (HashMap) mirparaOptions.get("paths");
            this.setInstallationFolder((String) pathsOptions.get("installation_folder"));
            
            HashMap rnafoldOptions = (HashMap) mirparaOptions.get("rnafold");
            this.setRnaFoldFolder((String) rnafoldOptions.get("folder"));
            this.setRnaFoldDll((String) rnafoldOptions.get("dll"));
            
            HashMap dataOptions = (HashMap) mirparaOptions.get("data");
            this.setDataFolder((String) dataOptions.get("folder"));
            this.setMirbaseDataFile((String) dataOptions.get("data_file"));
            
            pathToMirbaseData = installationFolder + FileSeparator + this.getDataFolder() + FileSeparator + this.getMirbaseDataFile();
            pathToModelData = installationFolder + FileSeparator + this.getModelFolder() + FileSeparator 
              + this.getModel() + "_" + this.getLevel() + ".model";
            
            
            HashMap programParams = (HashMap) mirparaOptions.get("program_params");
            this.setMaxNumOfPredictionListEntries((Integer) programParams.get("prediction_list_size"));
            
            HashMap predictionOptions = (HashMap) mirparaOptions.get("prediction_params");
            this.setWindow((Integer) predictionOptions.get("window"));
            this.setStep((Integer) predictionOptions.get("step"));
            this.setStart((Integer) predictionOptions.get("start"));
            this.setDistance((Integer) predictionOptions.get("distance"));
            this.setCutoff((Double) predictionOptions.get("cutoff"));
            this.setLevel((Integer) predictionOptions.get("level"));
                       
        }
        catch(FileNotFoundException e){
            logger.fatal("Configuration file " + this.getConfigurationFile() + " not found");
        }
    }
    
    /**
     * summarize run parameters
     * 
     * @return String
     */
    public String reportParameters()
    {
        String summary = 
            "window size               :\t" + window + "\n"
          + "step size                 :\t" + step + "\n"
          + "start                     :\t" + start + "\n"
          + "distance                  :\t" + distance + "\n"
          + "cutoff                    :\t" + cutoff + "\n"
          + "model                     :\t" + model + "\n"
          + "level                     :\t" + level;
        
        return summary;
    }
    
    
    /**
     * report file and folder settings
     * 
     * @return String
     */
    public String reportConfiguration(){
        String summary = 
            "configuration file        :\t"  + configurationFile + "\n"
          + "query file                :\t"  + inputFilename + "\n"
          + "installation folder       :\t"  + installationFolder + "\n"
          + "model folder              :\t"  + modelFolder + "\n"
          + "RNAFold folder            :\t"  + rnaFoldFolder + "\n"
          + "data folder               :\t"  + dataFolder + "\n"
          + "rnaFold File              :\t"  + rnaFoldDll + "\n"
          + "mirbase data file         :\t"  + mirbaseDataFile + "\n"
          + "path to library           :\t"  + pathToLibrary + "\n"
          + "path to miRBase data      :\t"  + pathToMirbaseData + "\n"
          + "path to model data        :\t"  + pathToModelData + "\n"
          + "prediction list size      :\t"  + maxNumOfPredictionListEntries + "\n"
          + "\n";

    
        return summary;
        
    }
    
    /**
     * @return the installationFolder
     */
    public String getInstallationFolder() {
        return installationFolder;
    }

    /**
     * @param installationFolder the installationFolder to set
     */
    public void setInstallationFolder(String installationFolder) {
        this.installationFolder = installationFolder;
    }

    /**
     * @return the modelFolder
     */
    public String getModelFolder() {
        return modelFolder;
    }

    /**
     * @param modelFolder the modelFolder to set
     */
    public void setModelFolder(String modelFolder) {
        this.modelFolder = modelFolder;
    }

    /**
     * @return the rnaFoldFolder
     */
    public String getRnaFoldFolder() {
        return rnaFoldFolder;
    }

    /**
     * @param rnaFoldFolder the rnaFoldFolder to set
     */
    public void setRnaFoldFolder(String rnaFoldFolder) {
        this.rnaFoldFolder = rnaFoldFolder;
    }

    /**
     * @return the dataFolder
     */
    public String getDataFolder() {
        return dataFolder;
    }

    /**
     * @param dataFolder the dataFolder to set
     */
    public void setDataFolder(String dataFolder) {
        this.dataFolder = dataFolder;
    }

    /**
     * @return the rnaFoldDll
     */
    public String getRnaFoldDll() {
        return rnaFoldDll;
    }

    /**
     * @param rnaFoldDll the rnaFoldDll to set
     */
    public void setRnaFoldDll(String rnaFoldDll) {
        this.rnaFoldDll = rnaFoldDll;
    }

    /**
     * @return the mirbaseDataFile
     */
    public String getMirbaseDataFile() {
        return mirbaseDataFile;
    }

    /**
     * @param mirbaseDataFile the mirbaseDataFile to set
     */
    public void setMirbaseDataFile(String mirbaseDataFile) {
        this.mirbaseDataFile = mirbaseDataFile;
    }

    /**
     * @return the configurationFile
     */
    public String getConfigurationFile() {
        return configurationFile;
    }

    /**
     * @param configurationFile the configurationFile to set
     */
    public void setConfigurationFile(String configurationFile) {
        this.configurationFile = configurationFile;
    }

    /**
     * @return the pathToModelData
     */
    public String getPathToModelData() {
        return pathToModelData;
    }

    /**
     * @param pathToModelData the pathToModelData to set
     */
    public void setPathToModelData(String pathToModelData) {
        this.pathToModelData = pathToModelData;
    }

    /**
     * @return the pathToLibrary
     */
    public String getPathToLibrary() {
        return pathToLibrary;
    }

    /**
     * @return the pathToMirbaseData
     */
    public String getPathToMirbaseData() {
        return pathToMirbaseData;
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
     * @return the maxNumOfPredictionListEntries
     */
    public int getMaxNumOfPredictionListEntries() {
        return maxNumOfPredictionListEntries;
    }

    /**
     * @param maxNumOfPredictionListEntries the maxNumOfPredictionListEntries to set
     */
    public void setMaxNumOfPredictionListEntries(int maxNumOfPredictionListEntries) {
        this.maxNumOfPredictionListEntries = maxNumOfPredictionListEntries;
    }

}
