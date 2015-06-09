/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import org.apache.logging.log4j.Logger;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


import org.apache.logging.log4j.LogManager;
import org.yaml.snakeyaml.Yaml;

/**
 *
 * @author weibo
 */
public class JmiCMD {

    static Options options = new Options();   
    static Yaml yaml = new Yaml();

    private static String test=null;
    private static double version=4.0;
    private static String configFile;

    static Logger logger = LogManager.getRootLogger();    
    
    public JmiCMD() {
        
        options.addOption("h", "help", false, "show help.");

        options.addOption("m", "model",         true, "model type: Animal (A), Plant (P), Virus (V) or Overall (O)");
        options.addOption("l", "level",         true, "level: the ratio of negative_data:positive_data (1 to 20: default 1)");
        options.addOption("s", "step",          true, "step size when splitting long sequences");
        options.addOption("w", "window",        true, "window size when splitting long sequences");
        options.addOption("d", "distance",      true, "minimum length for hairpin");
        options.addOption("e", "cutoff",        true, "minimum classification score for prediction (0 to 1)");
        options.addOption("i", "input-file",    true, "FASTA input file for query sequences");
        options.addOption("o", "output-folder", true, "output folder for prediction results");
        options.addOption("t", "test",          true, "run test : svm (S) or hairpin (H)");
        
    }
    
    public static void main(String[] args) throws IOException {
        long start = System.currentTimeMillis();
        
        printBanner();
        
        PipeLine pl=new PipeLine();

        //readArguments(args,pl);
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
        
        //pl.run2();

        long end = System.currentTimeMillis();
        long time = end - start;
        logger.info("Completed: total time:" + time/1000+"s");
    }

    public static void exit_with_help(){
                printUsage();
		System.exit(0);
    }
    
    public static void printBanner(){
        logger.info("\n\n\n"
                + "    =======================================================================\n"
                + "    |    JmiPara "+version+" performs microRNA prediction on input sequences      |\n"
                + "    =======================================================================\n\n\n");
        logger.info("The program can be run via a GUI or from the command line: \n"
                + "to start GUI: java -jar JmiPara2.jar or double click the package JmiPara2\n"
                + "to start CMD: java -cp JmiPara2.jar jmipara2.JmiCMD [options] input_file [output_directory]\n\n"
                + "to view results\n"
                + "start MirViewer: java -cp JmiPara2.jar jmipara2.MirViewer\n");
        
        logger.info("*** report bugs to simon.rayner@medisin.uio.no\n");
    }
    
    public static void printUsage(){
        System.out.println("Usage: \n"
                + "[1]start GUI: java -jar JmiPara2.jar or double click the package JmiPara2\n"
                + "[2]start CMD: java -cp JmiPara2.jar jmipara2.JmiCMD [options] input_file [output_directory]\n"
                + "[3]start MirViewer: java -cp JmiPara2.jar jmipara2.MirViewer\n\n"
		+"CMD options:\n"
		+"\t-m model:     the model used in prediction, can be animal, plant, \n"
                +"\t              virus or overall (default overall)\n"
		+"\t-l level:     the value of negative_data/positive_data, 1~20(default 1)\n"
                +"\t-s step size: step size in spliting long sequences\n"
                +"\t-w window:    window size in spliting long sequences\n"
                +"\t-d distance:  the minimal length of a hairpin\n"
                +"\t-c cutoff:    the minimal possibility that testing seq is miRNA\n"
                +"\t-t test:      1: test a given svm matrix file\n"
                +"\t              2: given a hairpin sequence, test whether a given site \n"
                +"\t                 have miRNA with given size (default 0)\n"
                +"\t-h help\n\n"+
                 "input_file:        the input file, should be fasta format at this version\n"
                +"output_directory:  the output directory, default is current directory\n\n"
                +"[References]\n"
                + "MiRPara: a SVM-based software tool for prediction of most probable microRNA"
                + " coding regions in genome scale sequences. Wu Y., Wei B., Liu H., Li T., "
                + "Rayner S. BMC Bioinformatics. 2011 Apr 19; 12(1):107\n\n");
    }

    private static void exit_with_error(){
        System.out.println("use -h to see detail help");
        System.exit(1);
    }

    private static void parseArguments(String[] args, PipeLine p){
        options.addOption("m", "model",         true, "model type: Animal (A), Plant (P), Virus (V) or Overall (O)");
        options.addOption("l", "level",         true, "level: the ratio of negative_data:positive_data (1 to 20: default 1)");
        options.addOption("s", "step",          true, "step size when splitting long sequences");
        options.addOption("w", "window",        true, "window size when splitting long sequences");
        options.addOption("d", "distance",      true, "minimum length for hairpin");
        options.addOption("c", "cutoff",        true, "minimum classification score for prediction (0 to 1)");
        options.addOption("i", "input-file",    true, "FASTA input file for query sequences");
        options.addOption("o", "output-folder", true, "output folder for prediction results");
        options.addOption("t", "test",          true, "run test example ");
        options.addOption("f", "config",        true, "configuration file for miRPara");

        logger.trace("parsing args");
        CommandLineParser parser = new BasicParser();
        CommandLine cmd = null;

        
        try {
                cmd = parser.parse(options, args);
                if (cmd.hasOption("h"))
                    print_help();

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

                if (cmd.hasOption("c")){
                    logger.info("cutoff set to " + cmd.getOptionValue("c"));
                    p.setCutoff(Double.parseDouble(cmd.getOptionValue("c")));
                }

                if (cmd.hasOption("t")){                    
                    test = cmd.getOptionValue("t").toLowerCase();
                    logger.info("test set to " + cmd.getOptionValue("l"));
                }
                
                if (cmd.hasOption("f")){
                    configFile = cmd.getOptionValue("f");
                    logger.info("Configuration file is " + configFile);
                    p.setConfigFilename(cmd.getOptionValue("f"));
                }

        } catch (ParseException e) {

         logger.fatal("Failed to parse comand line properties", e);

         print_help();

        }

        
    }
            
    
    private static void print_help(){
        printBanner();
        HelpFormatter formatter = new HelpFormatter();

        formatter.printHelp("command line options", options);
        
        System.exit(0);
        
    }
            
    private static void readArguments(String[] argv, PipeLine p){
        int i;
        if(argv.length==0) exit_with_help();
        for(i=0;i<argv.length;i++){
            if (argv[i].charAt(0) != '-')	break;
            if(++i>=argv.length && !(argv[i-1].equals("-h"))){
                System.out.println("arguments error!");
                exit_with_error();
            }
            switch(argv[i-1].charAt(1)){
                case 'm': p.setModel(argv[i]);	break;
		case 'l': p.setLevel(Integer.parseInt(argv[i]));	break;
                case 's': p.setStep(Integer.parseInt(argv[i]));  break;
                case 'w': p.setWindow(Integer.parseInt(argv[i])); break;
                case 'd': p.setDistance(Integer.parseInt(argv[i])); break;
                case 'c': p.setCutoff(Double.parseDouble(argv[i])); break;
                case 'h': exit_with_help();
		default:
			System.err.println("unknown option");
			exit_with_error();
            }
	}
        if(!(p.getModel().equals("animal") || p.getModel().equals("plant")
                || p.getModel().equals("virus") || p.getModel().equals("overall"))){
            System.err.println("unknow taxonomy!");
            exit_with_error();
        }
        if(p.getLevel()<0){
            System.err.println("level should be greater than 0");
            exit_with_error();
        }
        if(p.getWindow()<0){
            System.err.println("window should not be less than 0");
            exit_with_error();
        }
        if(p.getStep()<0){
            System.err.println("step should not be less than 0");
            exit_with_error();
        }
        if(p.getDistance()<0){
            System.err.println("distance should not be less than 0");
            exit_with_error();
        }
        if(p.getCutoff()<0 || p.getCutoff()>1){
            System.err.println("cutoff should be a value between 0 to 1");
            exit_with_error();
        }
        // determine filenames, at least two filenames
        if(i>=argv.length){
            System.err.println("lack input file!");
            exit_with_error();
        }
        p.setFilename(new File(argv[i]));
        if(i<argv.length-1)
            p.setWorkingDir(argv[i+1]);
    }

}



