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
 * main class for performing miRNA prediction
 * user needs to provide 
 *   input file containing query sequence in fasta format
 *   destination folder for results
 *   configuration file in YAML format specifying program and run parameters 
 * 
 * @author Wei Bo and Simon Rayner
 * 
 */
public class MirParaCmd {

    static String               FileSeparator   = System.getProperty("file.separator");
    static Options              options         = new Options();   
    

    private static String       test            = "";
    private static String       configFile;

    static Logger               logger          = LogManager.getRootLogger();    
    
    public MirParaCmd() {
        
        
    }
    
    
    
    
    public static void main(String[] args) throws IOException, Exception {
        
        long start = System.currentTimeMillis();
        
        printBanner();
        
        MiParaPipeLine pl = new MiParaPipeLine();
        parseArguments(args, pl);
        switch(pl.getAction()){
            case MiParaPipeLine.ACTION_TEST_SVM:
                logger.info("\n\n\nSVM TEST");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("\n\n\nSVM TEST");
                logger.info(new String(new char[80]).replace("\0", "x"));
                pl.testProgram();
                break;

            case MiParaPipeLine.ACTION_TEST_PREDICTION:
                logger.info("\n\n\n");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("\n\n\ntest miRNA prediction");
                logger.info(new String(new char[80]).replace("\0", "x"));
                pl.testGivenMir();
                break;
                        
            case MiParaPipeLine.ACTION_PREDICT_MIRNAS:
                logger.info("\n\n\n");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("\n\n\nPredict miRNAs from input sequences");
                logger.info(new String(new char[80]).replace("\0", "x"));
                pl.predictMiRNAsInQuerySequencesWithSplit();
                break;
            
            case MiParaPipeLine.ACTION_PREDICT_HAIRPINS:
                logger.info("\n\n\n");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("\n\n\nPredict hairpins from input sequences");
                logger.info(new String(new char[80]).replace("\0", "x"));
                pl.predictHairpins();
                break;
            
            case MiParaPipeLine.ACTION_PERFORM_TRAINING:
                logger.info("\n\n\n");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("\n\n\ntrain models from input data");
                logger.info(new String(new char[80]).replace("\0", "x"));
                break;
            
            case MiParaPipeLine.ACTION_CHARACTERIZE_MIRBASE:
                logger.info("\n\n\n");
                logger.info(new String(new char[80]).replace("\0", "x"));
                logger.info("       train models from input data");
                logger.info(new String(new char[80]).replace("\0", "x"));
                pl.parseMiRBaseEMBLData();
                break;
            
            default:
                logger.info("unrecognized action <" + pl.getAction() + ">");
                logger.error("unrecognized action <" + pl.getAction() + ">");
                throw new RuntimeException("unrecognized action <" + pl.getAction() + ">");
        
        }

        
        
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
                + "start MirViewer: java -cp JmiPara2.jar jmipara2.MirViewer\n")
                ;
        logger.info("OS NAME      is <" + System.getProperty("os.name") + ">\n");
        logger.info("LIBRARY_PATH is <" + System.getProperty("java.library.path") + ">\n");
 
        
        logger.info("\n\n*** report bugs to simon.rayner@medisin.uio.no\n");
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
    private static void parseArguments(String[] args, MiParaPipeLine p) throws IOException{
        options.addOption("i", "input-file",    true,  "FASTA input file for query sequences");
        options.addOption("o", "output-folder", true,  "output folder for prediction results");
        options.addOption("t", "test",          true,  "run test example ");
        options.addOption("c", "config",        true,  "configuration file for miRPara");
        options.addOption("a", "action",        true,  "action to perform");
        options.addOption("h", "help  ",        false, "print help");

        logger.trace("parsing args");
        CommandLineParser parser = new BasicParser();
        CommandLine cmd = null;

        
        try {
                cmd = parser.parse(options, args);
                if (cmd.hasOption("h")){
                    print_help();
                    System.out.print(p.reportAvailableActions());
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
                
                if (cmd.hasOption("a")){
                    logger.info("requested action is " + cmd.getOptionValue("a"));
                    p.setAction(cmd.getOptionValue("a"));
                }
                    
                if (cmd.hasOption("c")){
                    configFile = cmd.getOptionValue("c");
                    logger.info("Configuration file is " + configFile);
                    p.setConfigFilename(configFile);
                }


                
                if(p.getConfigFilename()==null){
                    throw new ParseException("no configuration file was specified") ;                    
                }
                
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



