/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.File;
import java.io.IOException;
import org.apache.logging.log4j.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


import org.apache.logging.log4j.LogManager;


/**
 *
 * @author weibo
 */
public class JmiCMD {

    static String FileSeparator = System.getProperty("file.separator");
    static Options options = new Options();   
    

    private static String test="";
    //private static double version=4.0;
    private static String configFile;
    

    static Logger logger = LogManager.getRootLogger();    
    
    
    
    
    public JmiCMD() {
        
        
    }
    
    
    
    
    public static void main(String[] args) throws IOException {
        
        long start = System.currentTimeMillis();
        
        printBanner();
        
        MiRNAPredictionPipeLine pl = new MiRNAPredictionPipeLine();
        parseArguments(args, pl);

        if(test.equals("svm")){
            
            logger.info("\n\n\nSVM TEST");
            logger.info(new String(new char[80]).replace("\0", "x"));
            logger.info("\n\n\nSVM TEST");
            logger.info(new String(new char[80]).replace("\0", "x"));
            pl.testProgram();
            System.exit(0);
            
        }
        if(test.equals("hairpin")){
            
            pl.testGivenMir();
            System.exit(0);
            
        }
            String se= "AAGCUGGCAUUCUAUAUAAGAGAGAAACUACACGCAGCGCCUCAUUUUGUGGGUCA"
              + "CCAUAUUCUUGGGAACAAGAGCUACAGCAUGGGGCAAAUCUUUCUGUUCCCAAUCCUCUGGGA"
              + "UUCUUUCCCGAUCACCAGUUGGACCCUGCGUUUGGAGCCAACUCAAACAAUCCAGAUUGGGAC"
              + "UUCAACCCCAACAAGGAUCACUGGCCAGAGGCAAAUCAGGUAGGAGCGGGAGCAUUCGGGCCA"
              + "GGGUUCACCCC";
            logger.info("MLEE:" + MfeFoldRNA.foldSequence(se));
        
//      logger.info("STEM");
//      pl.predictMiRNAsInQuerySequences();
        logger.info("PIPE");
        pl.predictMiRNAsInQuerySequencesWithSplit();
        
        long end = System.currentTimeMillis();
        long time = end - start;
        
        logger.info("Completed: total time:" + time/1000+"s");
        
    }


    
    /**
     * print program banner
     * 
     */
    public static void printBanner(){
        logger.info("\n\n\n"
                + "    =======================================================================\n"
                + "    |    miRPara4j: performs microRNA prediction on input sequences        |\n"
                + "    =======================================================================\n\n\n");
        logger.info("The program can be run via a GUI or from the command line: \n"
                + "to start GUI: java -jar JmiPara2.jar or double click the package JmiPara2\n"
                + "to start CMD: java -cp JmiPara2.jar jmipara2.JmiCMD [options] input_file [output_directory]\n\n"
                + "to view results\n"
                + "start MirViewer: java -cp JmiPara2.jar jmipara2.MirViewer\n");
        
        logger.info("*** report bugs to simon.rayner@medisin.uio.no\n");
    }
    
    
    
        
    
    /**
     * parse program arguments
     * 
     * many of these arguments should probably be put inside the configuration file
     * (i.e., the prediction parameters, they probably don't need to be specified
     *  every run)
     * 
     * @param args
     * @param p
     * @throws IOException 
     */
    private static void parseArguments(String[] args, MiRNAPredictionPipeLine p) throws IOException{
        options.addOption("m", "model",         true, "model type: Animal (A), Plant (P), Virus (V) or Overall (O)");
        options.addOption("l", "level",         true, "level: the ratio of negative_data:positive_data (1 to 20: default 1)");
        options.addOption("s", "step",          true, "step size when splitting long sequences");
        options.addOption("w", "window",        true, "window size when splitting long sequences");
        options.addOption("d", "distance",      true, "minimum length for hairpin");
        options.addOption("b", "baseline",      true, "minimum classification score for prediction (0 to 1)");
        options.addOption("r", "root-folder",   true, "root folder for data");
        options.addOption("i", "input-file",    true, "FASTA input file for query sequences");
        options.addOption("o", "output-folder", true, "output folder for prediction results");
        options.addOption("t", "test",          true, "run test example ");
        options.addOption("c", "config",        true, "configuration file for miRPara");

        logger.trace("parsing args");
        CommandLineParser parser = new BasicParser();
        CommandLine cmd = null;

        
        try {
                cmd = parser.parse(options, args);
                if (cmd.hasOption("h"))
                    print_help();

                if (cmd.hasOption("m")) {
                    logger.info("model set to " + cmd.getOptionValue("m"));
                    p.setModel(cmd.getOptionValue("m"));	
                }                    

                if (cmd.hasOption("l")) {
                    logger.info("level set to " + cmd.getOptionValue("l"));
                    p.setLevel(Integer.parseInt(cmd.getOptionValue("l")));	
                }                    

                if (cmd.hasOption("s")){
                    logger.info("step set to " + cmd.getOptionValue("s"));
                    p.setStep(Integer.parseInt(cmd.getOptionValue("s")));
                }

                if (cmd.hasOption("w")){
                    logger.info("window set to " + cmd.getOptionValue("w"));
                    p.setWindow(Integer.parseInt(cmd.getOptionValue("w")));
                }

                if (cmd.hasOption("d")){
                    logger.info("distance set to " + cmd.getOptionValue("d"));
                    p.setDistance(Integer.parseInt(cmd.getOptionValue("d")));
                }

                if (cmd.hasOption("b")){
                    logger.info("baseline score set to " + cmd.getOptionValue("b"));
                    p.setCutoff(Double.parseDouble(cmd.getOptionValue("b")));
                }

                if (cmd.hasOption("t")){                    
                    logger.info("test set to " + cmd.getOptionValue("l"));
                    test = cmd.getOptionValue("t").toLowerCase();
                }
                
                if (cmd.hasOption("i")) {
                    logger.info("input file set to " + cmd.getOptionValue("i"));
                    p.setInputFilename(cmd.getOptionValue("i"));	
                }                    

                if (cmd.hasOption("o")){
                    logger.info("output folder is" + cmd.getOptionValue("o"));
                    p.setOutputFolder(cmd.getOptionValue("o"));
                }
                
                if (cmd.hasOption("r")) {
                    logger.info("root folder is" + cmd.getOptionValue("r"));
                    p.setInputFilename(cmd.getOptionValue("r") + FileSeparator + cmd.getOptionValue("i"));	
                    p.setOutputFolder(cmd.getOptionValue("r") + FileSeparator + cmd.getOptionValue("o"));
                    logger.info("input file now set to " + cmd.getOptionValue("i"));
                    logger.info("output folder now set to " + cmd.getOptionValue("o"));                    
                }                    
                
                if (cmd.hasOption("c")){
                    configFile = cmd.getOptionValue("c");
                    logger.info("Configuration file is " + configFile);
                    p.setConfigFilename(configFile);
                }
                else
                    throw new ParseException("no configuration file was specified") ;
                
                if(new File(p.getConfigFilename()).exists() == false){
                    throw new IOException("configuration file <" + p.getConfigFilename() + "> does not exist");
                }
                    
                
                if(new File(p.getInputFilename()).exists()== false)
                {
                    throw new IOException("input file <" + p.getInputFilename() + "> does not exist");
                }
                
                if(new File(p.getOutputFolder()).exists() == false)
                {
                    if (new File(p.getOutputFolder()).getParentFile().exists() == false)
                        throw new IOException("parent file <" + new File(p.getOutputFolder()).getParentFile() + "> does not exist");
                    
                    new File(p.getOutputFolder()).mkdir();
                    logger.info("create results folder <" + p.getOutputFolder() + ">");
                }

        } catch (ParseException e) {

            logger.fatal("Failed to parse command line properties", e);
            print_help();

        }

        
    }
            
    
    
    
    /**
     * print program usage
     * 
     */
    public static void print_help(){
        printBanner();
        HelpFormatter formatter = new HelpFormatter();

        formatter.printHelp("command line options", options);
        
        System.exit(0);
        
    }
            

}



