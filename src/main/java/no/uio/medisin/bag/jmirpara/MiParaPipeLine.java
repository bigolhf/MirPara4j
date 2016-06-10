
package no.uio.medisin.bag.jmirpara;

import no.uio.medisin.bag.core.FoldableRNASequence;
import no.uio.medisin.bag.core.PriMiRNA;
import no.uio.medisin.bag.core.FeatureRange;
import no.uio.medisin.bag.core.CharPriMiRNA;
import no.uio.medisin.bag.core.SimpleSequenceSet;
import no.uio.medisin.bag.core.SimpleSeq;
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
import no.uio.medisin.bag.jmirpara.utilities.MiRBaseEMBLDataFileParser;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.yaml.snakeyaml.Yaml;
import org.apache.commons.io.FilenameUtils;

/**
 * this class performs various types of processing against different datatypes
 * related to miRNAs
 * 
 * @author simon rayner
 */
public class MiParaPipeLine {
    
    static Logger                   logger                          = LogManager.getRootLogger();    
    

    private static final String     FILE_SEPARATOR                  = System.getProperty("file.separator");
    private static final String     DISCOUNT_PARAMETER_LIMIT        = "NA";
    
    public static final String      ACTION_TEST_SVM                 = "S";
    public static final String      ACTION_TEST_PREDICTION          = "P";
    public static final String      ACTION_PREDICT_MIRNAS           = "M";
    public static final String      ACTION_PREDICT_HAIRPINS         = "H";
    public static final String      ACTION_PERFORM_TRAINING         = "T";
    public static final String      ACTION_CHARACTERIZE_MIRBASE     = "C";

    private static final String     ID_PREDICTION_PARAMETERS        = "predictionParameters";
    private static final String     ID_PROGRAM_PARAMETERS           = "programParameters";
    private static final String     ID_REFDATA_PARAMETERS           = "refdataParameters";
    private static final String     ID_RNAFOLD_PARAMETERS           = "rnafoldParameters";
    private static final String     ID_PATHS                        = "filePaths";
    private static final String     ID_MODEL_PARAMETERS             = "modelParameters";
    private static final String     ID_MIRBASE_PARSE_PARAMETERS     = "mirbaseParameters";
    
    
    private static final String     ID_MIRBASE_DATAFILE             = "mirbaseDataFile";
    private static final String     ID_REFERENCE_DATAFOLDER         = "pathToReferenceData";
    private static final String     ID_RNAFOLD_FOLDER               = "rnafoldFolder";
    private static final String     ID_RNAFOLD_LIBRARY_FILE         = "rnafoldLibraryFile";
    private static final String     ID_INSTALLATION_FOLDER          = "installationFolder";
    private static final String     ID_MODEL_FOLDER                 = "modelFolder";


    private static final String     ID_PREDICTION_LIST_SIZE         = "predictionListSize";
    private static final String     ID_WINDOW_SIZE                  = "windowSize";
    private static final String     ID_STEP_SIZE                    = "stepSize";
    private static final String     ID_START_POSITION               = "startPosition";
    private static final String     ID_MIN_HAIRPIN_LENGTH           = "minHairpinLength";
    private static final String     ID_SCORE_CUTOFF                 = "predictionScoreCutoff";
    private static final String     ID_MODEL_LEVEL                  = "modelLevel";
    private static final String     ID_MIN_MIRNA_LENGTH             = "minMirnaLength";
    private static final String     ID_MAX_MIRNA_LENGTH             = "maxMirnaLength";
    private static final String     ID_PREDICTION_MODEL             = "model";
    private static final String     ID_HOSTS                        = "hosts";
    
    
    
    
    
    private static final int        MIN_MIRNA_LEN = 20;
    private static final int        MAX_MIRNA_LEN = 25;
    
    private static final String[]   KNOWN_MODELS = new String[] {"m", "p", "v", "o", "animal", "plant", "virus", "overall" };
    private static final String[]   KNOWN_TESTS = new String[] {"svm", "hairpin", };
    
    
    private String                  action;                     // type of action requested by user

    private String                  configFilename;

        
    private String                  inputFilename;
    private String                  outputFolder;

    static  Yaml                    yaml = new Yaml();
    
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

    
    miRParaPredParams               miRParaPredParams;
    private int                     window              = 120;
    private int                     step                = 5000;
    private int                     start               = 1;
    private int                     distance            = 60;
    private double                  cutoff              = 0.8;
    private String                  model               = "overall";
    private int                     level               = 1;
    private int                     mirnaMinLen         = 1;
    private int                     mirnaMaxLen         = 40;
    
    private File                    workingDir          = new File(".");
    private boolean                 append              = false;
    
    private String                  hosts               ="";
    //private double                  progress;

    private ArrayList<String>       results     = new ArrayList<>(); // <- is this necessary ?
    


    private ArrayList<SimpleSeq>    seqList;
    private ArrayList<PriMiRNA>     priList;
    private String[]                uniquePriMiRNAIDs;
    private ArrayList<HashMap>      predictionList = new ArrayList<>();
    
    private String test;

    public MiParaPipeLine() {  

    }

    
    
    /**
     * parse out and characterize pre-miRNA/miRNAs specified in an EMBL file from 
     * miRBase.
     * input can be filtered by host and whether feature is supported by experimental data
     * 
     * @throws IOException 
     */
    public void parseMiRBaseEMBLData() throws IOException, Exception{
        MiRBaseEMBLDataFileParser miRBaseEMBLParser = new MiRBaseEMBLDataFileParser();
        miRBaseEMBLParser.setEmblFileName(inputFilename);
        miRBaseEMBLParser.parseEMBLFile();
    }
    
    
    
    
    
    /**
     * predict secondary structures for a set of input sequences.
     * It is assumed they represent pre-miRNA/miRNA candidates and will
     * form a hairpin structure.
     * 
     * @throws IOException
     * @throws Exception 
     */
    public void predictHairpins() throws IOException, Exception{
        
        uniquePriMiRNAIDs = new String[0]; // not quite sure what this for. something to do with tracking duplicates....
        logger.info("predict hairpins for input sequence set");
        initializePipeline();
        verifyData();
        
        logger.info("input file is <" + this.getInputFilename() + ">");
        SimpleSequenceSet querySeqSet = new SimpleSequenceSet(new File(getInputFilename()));
        querySeqSet.readFromFastaFile();
        logger.info("read " + querySeqSet.getSeqs().size() + " sequences");
        
        String hairpinTSVFile = this.cleanPath(this.getOutputFolder() 
                + FILE_SEPARATOR + FilenameUtils.getBaseName(this.getInputFilename()) 
                + ".hairpins.tsv");
        String currSeqName = "";
        try{
            BufferedWriter bwPri = new BufferedWriter(new FileWriter(new File(hairpinTSVFile )));


            Boolean firstPriMiRNA = true;
            for(SimpleSeq querySeq: querySeqSet.getSeqs()){

                this.SetOutfilePrefix(querySeq.getId());
                logger.info("-- processing sequence <" + querySeq.getName() + "> : " + querySeq.getLength() + " nt");
                ArrayList<SimpleSeq> fragmentList = querySeq.splitSequence(step, window);
                logger.info("-- broken into " + fragmentList.size() + " fragments");
                for(SimpleSeq fragment:fragmentList){
                    currSeqName = fragment.getName();
                    logger.info("folding fragment <" + fragment + ">");
                    FoldableRNASequence foldableRNASeq = new FoldableRNASequence(fragment);
                    ArrayList<PriMiRNA> pris = foldableRNASeq.foldAndScanSequence();

                    removeDuplicatePrimiRNAs(pris);


                    for(PriMiRNA priMiRNA: pris){
                        CharPriMiRNA charPriMiRNA = new CharPriMiRNA(priMiRNA);
                        charPriMiRNA.characterize();
                        if(firstPriMiRNA){
                            bwPri.write(charPriMiRNA.getPriRNA().printFeatureSetKeysAsTSV() + "\n");
                            firstPriMiRNA=false;
                        }
                        bwPri.write(charPriMiRNA.getPriRNA().printFeatureSetValuesAsTSV().replace("\n", "") + "\n");
                    }
                }

            }   
            bwPri.close();
            
        }
        catch(IOException exIO){
            logger.info("IO error while predicting hairpins for sequence <" + currSeqName + ">");
            logger.error("IO error while predicting hairpins for sequence <" + currSeqName + ">");
            throw new IOException("IO error while predicting hairpins for sequence <" + currSeqName + ">\n"  + exIO);
        }
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
     * all run parameters can be specified in the YAML file that is required as input
     * 
     * @throws IOException 
     * 
     */
    public void predictMiRNAsInQuerySequences() throws IOException{
        
        initializePipeline();
        verifyData();
        verifyPredictionData();

        logger.info("loading model data file <" + this.getPathToModelData() + "");
        SVMToolKit.loadModel(this.getPathToModelData());

        logger.info("loading miRBase data file <" + this.getPathToModelData() + "");
        MiRNAPredictionsSerializer.loadMirBaseData(new File(this.getPathToMirbaseData()));     
        
        SimpleSequenceSet querySeqSet = new SimpleSequenceSet(new File(getInputFilename()));
        
        for(SimpleSeq querySeq: querySeqSet.getSeqs()){
            
            setOutfileName(workingDir, new File(getInputFilename()), querySeq.getId());            

            predictionList = new ArrayList<>();
            uniquePriMiRNAIDs = new String[0];
            serializeResults(querySeq); // this will just write out header information
            
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
            
            FoldableRNASequence foldableRNASeq = new FoldableRNASequence(querySeq);
            ArrayList<PriMiRNA> priMiRNAsInFrag = foldableRNASeq.foldAndScanSequence(); //this.getWindow(), this.getStep(), this.getDistance()
            
            removeDuplicatePrimiRNAs(priMiRNAsInFrag);
            logger.info(priMiRNAsInFrag.size() + " pri-miRNAs were found in the query sequence...");
            
            int before = predictionList.size();
            logger.info(" ");
            logger.info("scanning pri-miRNAs for miRNAs...");
            for(PriMiRNA primiRNA : priMiRNAsInFrag){
                logger.info(".." + primiRNA.getId());
                findMiRNAsInPrimiRNA(primiRNA);                
            }
                
            switch(totalNumOfMiRNA){
                case 0:
                    logger.info("didn't find any miRNAs");
                    break;
                    
                case 1:
                    logger.info("found a total of 1 miRNA candidate ...");
                    break;
                    
                default:
                    logger.info("found a total of " + totalNumOfMiRNA + " miRNAs candidates...");
            }
            
            //output the results now to avoid memory leak
            if(predictionList.size()>100){
                logger.info(" ");
                logger.info("flush prediction list");
                logger.info(" ");
                totalNumOfMiRNA += predictionList.size();
                serializeResults(querySeq);
                predictionList = new ArrayList<>();
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
    
    
    
    
    /**
     * Looks for identical pri-miRNA ids. These IDs are based on position, so this will
     * identify identical entries.
     * 
     * @param primiRNAList 
     */
    private void removeDuplicatePrimiRNAs(ArrayList<PriMiRNA> primiRNAList){
        logger.info("removing duplicate structures");
        
        Iterator itPriList = primiRNAList.iterator();
        while(itPriList.hasNext()){
            PriMiRNA pri=(PriMiRNA)(itPriList.next());
            for(String id:uniquePriMiRNAIDs){
                if(pri.getId().equals(id)){
                    itPriList.remove();
                    break;
                }
            }
        }
        
        //update unique ID array;
        int n = primiRNAList.size();
        uniquePriMiRNAIDs = new String[n];
        int i=0;
        for(PriMiRNA pri:primiRNAList){
            uniquePriMiRNAIDs[i++] = pri.getId();
        }     
    }
    
    

    /**
     * Load query sequences in FASTA format
     * 
     * @param fn
     * @throws IOException 
     */
    private void loadFastaSequenceData(File fn) throws IOException{
        
        seqList=new SimpleSequenceSet(fn).getSeqs();
        
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
        CharPriMiRNA characterizedPrimiRNA = new CharPriMiRNA(priRNA);
        
        try{
            
            characterizedPrimiRNA.characterize();
            
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
        MiRNAPredictionsSerializer.outputFilePrefix = getOutputFolder();
    }

    
        
    public void SetOutfilePrefix(String seqName){
        MiRNAPredictionsSerializer.outputFilePrefix = this.outputFolder + FILE_SEPARATOR + seqName;
    }
    
    /**
     * output prediction results for specified sequence
     * 
     * @param seq
     * @throws IOException 
     */
    public void serializeResults(SimpleSeq seq) throws IOException {
        
        MiRNAPredictionsSerializer.serializePredictionDetails(predictionList, seq, append);
        MiRNAPredictionsSerializer.serializePredictionSummary(predictionList, seq, append);
        MiRNAPredictionsSerializer.serializePredictionsAsHTML(predictionList, seq);
        
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
    private void initializePipeline() throws IOException, RuntimeException{

        logger.info("reading configuration file");
        try{
            this.readConfigurationFile(configFilename);
            logger.info(this.reportConfiguration() + "\n");
            logger.info(this.reportParameters() + "\n");
            logger.info(this.reportRunSettings() + "\n");
        }
        catch(IOException|RuntimeException ex){
            
        }
        
        
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
         +  this.reportParameters() + "\n\n"
         +  StringUtils.repeat("*", 60) + "\n\n";

        return summary;
    }
    
    
    
    
    public String reportAvailableActions(){
/*
    public static final String      ACTION_TEST_SVM                 = "S";
    public static final String      ACTION_TEST_PREDICTION          = "P";
    public static final String      ACTION_PREDICT_MIRNAS           = "M";
    public static final String      ACTION_PREDICT_HAIRPINS         = "H";
    public static final String      ACTION_PERFORM_TRAINING         = "T";
    public static final String      ACTION_CHARACTERIZE_MIRBASE     = "C";
        
        */        

        return    "known actions are:" + "\n"
                + "  ACTION_TEST_SVM             (S)" + "\n"
                + "  ACTION_TEST_PREDICTION      (P)" + "\n"
                + "  ACTION_PREDICT_MIRNAS       (M)" + "\n"
                + "  ACTION_PREDICT_HAIRPINS     (H)" + "\n"
                + "  ACTION_PERFORM_TRAINING     (T)" + "\n"
                + "  ACTION_CHARACTERIZE_MIRBASE (C)" + "\n";
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
        if(Arrays.asList(KNOWN_MODELS).contains(model.toLowerCase()) == false)
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
        if(!Arrays.asList(KNOWN_TESTS).contains(test.toLowerCase())){
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
    public void predictMiRNAsInQuerySequencesWithSplit() throws IOException, RuntimeException{
  
        
        initializePipeline();
        verifyData();
        verifyPredictionData();
        
        logger.info("loading model data file <" + this.getPathToModelData() + "");
        SVMToolKit.loadModel(this.getPathToModelData());

        logger.info("loading miRBase data file <" + this.getPathToModelData() + "");
        MiRNAPredictionsSerializer.loadMirBaseData(new File(this.getPathToMirbaseData()));     
        
        SimpleSequenceSet querySeqSet = new SimpleSequenceSet(new File(getInputFilename()));
        for(SimpleSeq querySeq: querySeqSet.getSeqs()){
        //while(querySeqSet.hasSeq()){
            
            //SimpleSeq querySeq = querySeqSet.getOneSeq();  //each seq
            setOutfileName(workingDir, new File(getInputFilename()), querySeq.getId());
            this.SetOutfilePrefix(querySeq.getId());

            predictionList = new ArrayList<>();
            uniquePriMiRNAIDs = new String[0];
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
            logger.info("-- will break into " + (length-window)/step+1 + " fragments");
            ArrayList<SimpleSeq> fragmentList = querySeq.splitSequence(step, window);
            for(SimpleSeq fragment:fragmentList){
                FoldableRNASequence sl = new FoldableRNASequence(fragment);

                ArrayList<PriMiRNA> pris = sl.foldAndScanSequence();
                
            }
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
                FoldableRNASequence sl = new FoldableRNASequence(frag);

                ArrayList<PriMiRNA> pris = sl.foldAndScanSequence();
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
                    predictionList=new ArrayList<>();
                }
//                progress = Double.parseDouble(OutputMiRNAPredictions.decimal((i+1)*100.0/n));
//                print(Output.decimal((i+1)*100.0/n)+"%"+Output.backspace(Output.decimal((i+1)*100.0/n)+"%"));
                print(MiRNAPredictionsSerializer.decimal((i+1)*100.0/n)+"%\n");
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
                logger.info(MiRNAPredictionsSerializer.decimal(i*100.0/priList.size())+"%"+MiRNAPredictionsSerializer.backspace(MiRNAPredictionsSerializer.decimal(i*100.0/priList.size())+"%"));                
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
                if(FoldableRNASequence.hasMultipleLoops(pri)){
                   logger.info("W1: "+entry[0]+" is not a hairpin structure!");
                    continue;
                }
                CharPriMiRNA parser=new CharPriMiRNA(pri);
                parser.characterize();

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
     * @throws Exception≈
     * @throws IOException
     * 
     */
    public void readConfigurationFile(String configFile) throws IOException, RuntimeException{
        
        logger.info("loading configuration file <" + ">");

        miRParaPredParams                                   = new miRParaPredParams();        
        HashMap                         mirparaOptions      = new HashMap();
        
        setConfigurationFile(configFile);
        try{            
            mirparaOptions = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));
        }
        catch(IOException e){
            logger.fatal("error loading configuration file <" + this.getConfigurationFile() + ">");
            throw new IOException("error loading configuration file <" + this.getConfigurationFile() + ">");
        }
        
        
        logger.info("loading model configuration");      
        HashMap modelsOptions = (HashMap) mirparaOptions.get(ID_MODEL_PARAMETERS);
        this.checkStringNotNull(this.getModelFolder(), ID_MODEL_FOLDER);
        this.setModelFolder((String) modelsOptions.get(ID_MODEL_FOLDER));


        logger.info("loading path configuration");  
        HashMap pathsOptions = (HashMap) mirparaOptions.get(ID_PATHS);
        this.checkStringNotNull(this.getInstallationFolder(), ID_INSTALLATION_FOLDER);
        this.setInstallationFolder((String) pathsOptions.get(ID_INSTALLATION_FOLDER));


        logger.info("loading RNAFold configuration");  
        HashMap rnafoldOptions = (HashMap) mirparaOptions.get(ID_RNAFOLD_PARAMETERS);
        this.checkStringNotNull(this.getRnaFoldFolder(), ID_RNAFOLD_FOLDER);
        this.setRnaFoldFolder((String) rnafoldOptions.get(ID_RNAFOLD_FOLDER));
        this.checkStringNotNull(this.getRnaFoldDll(), ID_RNAFOLD_LIBRARY_FILE);
        this.setRnaFoldDll((String) rnafoldOptions.get(ID_RNAFOLD_LIBRARY_FILE));


        logger.info("loading data configuration");  
        HashMap dataOptions = (HashMap) mirparaOptions.get(ID_REFDATA_PARAMETERS);
        this.checkStringNotNull(this.getDataFolder(), ID_REFERENCE_DATAFOLDER);
        this.setDataFolder((String) dataOptions.get(ID_REFERENCE_DATAFOLDER));
        
        this.checkStringNotNull(this.getMirbaseDataFile(), ID_MIRBASE_DATAFILE);
        this.setMirbaseDataFile((String) dataOptions.get(ID_MIRBASE_DATAFILE));
            
            
        logger.info("loading program parameters");
        
        HashMap programParams = (HashMap) mirparaOptions.get(ID_PROGRAM_PARAMETERS);
        String chk = "";
        
        chk = checkParameter("Integer", ID_PREDICTION_LIST_SIZE, Integer.toString((Integer) programParams.get(ID_PREDICTION_LIST_SIZE)), "1", "NA", logger);
        if(chk!=null)
            this.setMaxNumOfPredictionListEntries((Integer) programParams.get(ID_PREDICTION_LIST_SIZE));

        
        logger.info("loading miRBase options");
        HashMap miRBaseParseOptions = (HashMap) mirparaOptions.get(ID_MIRBASE_PARSE_PARAMETERS);
        hosts = (String) miRBaseParseOptions.get(ID_HOSTS);
        
        
        logger.info("loading other program options");  
        
        HashMap predictionOptions = (HashMap) mirparaOptions.get(ID_PREDICTION_PARAMETERS);

        miRParaPredParams.setModel((String)programParams.get(ID_PREDICTION_MODEL));

        chk = checkParameter("Integer", ID_MIN_MIRNA_LENGTH, Integer.toString((Integer) programParams.get(ID_MIN_MIRNA_LENGTH)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setMinMiRNAlen((Integer) predictionOptions.get(ID_MIN_MIRNA_LENGTH));
            
        chk = checkParameter("Integer", ID_MAX_MIRNA_LENGTH, Integer.toString((Integer) programParams.get(ID_MAX_MIRNA_LENGTH)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setMaxMiRNAlen((Integer) predictionOptions.get(ID_MAX_MIRNA_LENGTH));
            
        chk = checkParameter("Integer", ID_WINDOW_SIZE, Integer.toString((Integer) programParams.get(ID_WINDOW_SIZE)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setWindow((Integer) predictionOptions.get(ID_WINDOW_SIZE));
            
        chk = checkParameter("Integer", ID_STEP_SIZE, Integer.toString((Integer) programParams.get(ID_STEP_SIZE)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setStep((Integer) predictionOptions.get(ID_STEP_SIZE));
            
        chk = checkParameter("Integer", ID_START_POSITION, Integer.toString((Integer) programParams.get(ID_START_POSITION)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setStart((Integer) predictionOptions.get(ID_START_POSITION));
            
        chk = checkParameter("Integer", ID_MIN_HAIRPIN_LENGTH, Integer.toString((Integer) programParams.get(ID_MIN_HAIRPIN_LENGTH)), "1", "NA", logger);
        if(chk!=null)
            miRParaPredParams.setMinHairpinLength((Integer) predictionOptions.get(ID_MIN_HAIRPIN_LENGTH));
            
        chk = checkParameter("Double", ID_SCORE_CUTOFF, Integer.toString((Integer) programParams.get(ID_SCORE_CUTOFF)), "0.0", "1.0", logger);
        if(chk!=null)
            miRParaPredParams.setCutoff((Double) predictionOptions.get(ID_SCORE_CUTOFF));
            
        chk = checkParameter("Integer", ID_MODEL_LEVEL, Integer.toString((Integer) programParams.get(ID_MODEL_LEVEL)), "1", "20", logger);
        if(chk!=null)
            miRParaPredParams.setLevel((Integer) predictionOptions.get(ID_MODEL_LEVEL));
           
        
        logger.info("done");
            
            
    }
    

    
    
    /**
     * verify general data files before we begin
     * this requires checking file paths (input/output folders) and required files exist
     * 
     * @return
     * @throws IOException 
     */
    private Boolean verifyData() throws IOException{
        
        pathToMirbaseData   = cleanPath(this.getDataFolder()  + FILE_SEPARATOR + this.getMirbaseDataFile());
        if(new File(pathToMirbaseData).exists() == false){
            logger.error("path to miRBase data < " + pathToMirbaseData +"> not found");
            logger.info("path to miRBase data < " + pathToMirbaseData +"> not found");
            throw new IOException("path to miRBase data < " + pathToMirbaseData +"> not found");
        }

           
        
        return true;
    }
    
    
    
    
    /**
     * verifies prediction specific datafiles before we begin
     * 
     * @return
     * @throws IOException 
     */
    private Boolean verifyPredictionData() throws IOException{
        
        pathToModelData     = cleanPath(this.getDataFolder() + FILE_SEPARATOR + this.getModelFolder() + FILE_SEPARATOR
          + this.getModel() + "_" + this.getLevel() + ".model");
        if(new File(pathToModelData).exists() == false){
            logger.error("path to model data < " + pathToModelData +"> not found");
            logger.info("path to model data < " + pathToModelData +"> not found");
            throw new IOException("path to model data < " + pathToModelData +"> not found");
        }
        return true;
    }
    
    /**
     * check whether this parameter has been defined
     * 
     * @param var
     * @param parameter
     * @throws Exception 
     */
    private void checkStringNotNull(String var, String parameter) throws RuntimeException{
        if(var==null){            
            logger.info("<" + parameter + "> has not been defined");
            throw new RuntimeException("<" + parameter + "> has not been defined");            
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
     * Check the specified parameter is valid
     * For numerical parameters this requires not null and within the specified 
     * upper and lower limits
     * 
     * @param className
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return
     * @throws Exception 
     */
    protected String checkParameter(String className, String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws RuntimeException{
        int             iVal;
        double          dVal;
        Boolean         bVal;
                
        if(Integer.class.getName().contains(className)){
            iVal = checkIntegerParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Integer.toString(iVal);
        }
        
        if(Double.class.getName().contains(className)){
            dVal = checkDoubleParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Double.toString(dVal);
        }
        
        if(Boolean.class.getName().contains(className)){
            bVal = checkBooleanParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Boolean.toString(bVal);
        }
                
        return null;
    }
    
    
    
    
    
    /**
     * 
     * check the specified Integer parameter is valid, i.e., not null and within the specified range
     * if there is only an upper or lower limit, the other limit can be set to NGSStep.DISCOUNT_PARAMETER_LIMIT
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return 
     */
    private int checkIntegerParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws RuntimeException{
        
        
        int             iVal;
        
        
        try{
            iVal = Integer.parseInt(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> is not an integer");
            logger.error(paramID + " <" + param + "> is not an integer");
            throw new NumberFormatException(paramID + " <" + param + "> is not an integer");
        }        

        if(lowerLimit.equals(DISCOUNT_PARAMETER_LIMIT) && upperLimit.equals(DISCOUNT_PARAMETER_LIMIT))
            return iVal;

        if(upperLimit.equals(DISCOUNT_PARAMETER_LIMIT)){
            if (iVal < Integer.parseInt(lowerLimit)){
                logger.info(paramID + " <" + param + "> must be >= " + lowerLimit);
                logger.error(paramID + " <" + param + "> must be >= " + lowerLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be >= " + lowerLimit);
            }
            else 
                return iVal;
        }

        if(lowerLimit.equals(DISCOUNT_PARAMETER_LIMIT)){
            if (iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param + "> must be > " + upperLimit);
                logger.error(paramID + " <" + param + "> must be > " + upperLimit);
                throw new IllegalArgumentException(paramID+ " <" + param + "> must be > " + upperLimit);
            }
            else 
                return iVal;
        }

        if(iVal <Integer.parseInt(lowerLimit) 
                || iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                logger.error(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                throw new IllegalArgumentException(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);

        }
        else{
            return iVal;            
        }
        
    }
    
    
    /**
     * 
     * check the specified Double parameter is valid, i.e., not null and within the specified range
     * if there is only an upper or lower limit, the other limit can be set to NGSStep.DISCOUNT_PARAMETER_LIMIT
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return 
     */
    private double checkDoubleParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws RuntimeException{
        
        
        double             dVal;
        try{
            dVal = Double.parseDouble(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> is not an integer");
            logger.error(paramID + " <" + param + "> is not an integer");
            throw new NumberFormatException(paramID + " <" + param + "> is not an integer");
        }        

        if(lowerLimit.equals(DISCOUNT_PARAMETER_LIMIT) && upperLimit.equals(DISCOUNT_PARAMETER_LIMIT))
            return dVal;

        if(upperLimit.equals(DISCOUNT_PARAMETER_LIMIT)){
            if (dVal < Double.parseDouble(lowerLimit)){
                logger.info(paramID + " <" + param + "> must be >= " + lowerLimit);
                logger.error(paramID + " <" + param + "> must be >= " + lowerLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be >= " + lowerLimit);
            }
            else 
                return dVal;
        }

        if(lowerLimit.equals(DISCOUNT_PARAMETER_LIMIT)){
            if (dVal > Double.parseDouble(upperLimit)){
                logger.info(paramID + " <" + param + "> must be > " + upperLimit);
                logger.error(paramID + " <" + param + "> must be > " + upperLimit);
                throw new IllegalArgumentException(paramID+ " <" + param + "> must be > " + upperLimit);
            }
            else 
                return dVal;
        }

        if(dVal <Double.parseDouble(lowerLimit) 
                || dVal > Double.parseDouble(upperLimit)){
                logger.info(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                logger.error(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be > " + upperLimit);

        }
        else{
            return dVal;            
        }
        
    }
    
    
    
    
    
    
    /**
     * check the specified Boolean parameter is valid. i.e., can it be cast as a Boolean
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return
     * @throws Exception 
     */
    private boolean checkBooleanParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws RuntimeException{
        
        
        boolean             bVal;
        try{
            bVal = Boolean.parseBoolean(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> cannot be cast as Boolean");
            logger.error(paramID + " <" + param + "> cannot be cast as Boolean");
            throw new NumberFormatException(paramID + " <" + param + "> cannot be cast as Boolean");
        }        
            
        return bVal;            
            
    }
    
    
    /**
     * parse out actionString to determine the requested user action
     * 
     * @param actionString
     * @return
     * @throws RuntimeException 
     */
    public String findAction(String actionString) throws RuntimeException{
        
        switch(actionString.toUpperCase()){
            
            case ACTION_TEST_SVM:
                action = ACTION_TEST_SVM;
                break;

            case ACTION_TEST_PREDICTION:
                action = ACTION_TEST_PREDICTION;
                break;
                        
            case ACTION_PREDICT_MIRNAS:
                action = ACTION_PREDICT_MIRNAS;
                break;
            
            case ACTION_PREDICT_HAIRPINS:
                action = ACTION_PREDICT_HAIRPINS;
                break;
            
            case ACTION_PERFORM_TRAINING:
                action = ACTION_PERFORM_TRAINING;
                break;
                
            case ACTION_CHARACTERIZE_MIRBASE:
                action = ACTION_CHARACTERIZE_MIRBASE;
                break;
            
            default:
                logger.info("unrecognized action <" + actionString + ">");
                logger.error("unrecognized action <" + actionString + ">");
                throw new RuntimeException("unrecognized action <" + actionString + ">");
        
        }
        return action;
    }
    
    
    
    /**
     * strip out duplicate folder delimiter from a file path
     * 
     * @param path
     * @return 
     */
    final String cleanPath(String path){
        return path.replace(FILE_SEPARATOR + FILE_SEPARATOR, FILE_SEPARATOR);
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

    /**
     * @return the minMiRNAlen
     */
    public int getMinMiRNAlen() {
        return getMirnaMinLen();
    }

    /**
     * @param minMiRNAlen the minMiRNAlen to set
     */
    public void setMinMiRNAlen(int minMiRNAlen) {
        this.setMirnaMinLen(minMiRNAlen);
    }

    /**
     * @return the maxMiRNAlen
     */
    public int getMaxMiRNAlen() {
        return getMirnaMaxLen();
    }

    /**
     * @param maxMiRNAlen the maxMiRNAlen to set
     */
    public void setMaxMiRNAlen(int maxMiRNAlen) {
        this.setMirnaMaxLen(maxMiRNAlen);
    }

    /**
     * @return the action
     */
    public String getAction() {
        return action;
    }

    /**
     * @param action the action to set
     */
    public void setAction(String action) {
        this.action = action;
    }

    /**
     * @return the mirnaMinLen
     */
    public int getMirnaMinLen() {
        return mirnaMinLen;
    }

    /**
     * @param mirnaMinLen the mirnaMinLen to set
     */
    public void setMirnaMinLen(int mirnaMinLen) {
        this.mirnaMinLen = mirnaMinLen;
    }

    /**
     * @return the mirnaMaxLen
     */
    public int getMirnaMaxLen() {
        return mirnaMaxLen;
    }

    /**
     * @param mirnaMaxLen the mirnaMaxLen to set
     */
    public void setMirnaMaxLen(int mirnaMaxLen) {
        this.mirnaMaxLen = mirnaMaxLen;
    }

    /**
     * @return the hosts
     */
    public String getHosts() {
        return hosts;
    }

    /**
     * @param hosts the hosts to set
     */
    public void setHosts(String hosts) {
        this.hosts = hosts;
    }

}
