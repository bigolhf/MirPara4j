package no.uio.medisin.bag.jmirpara.train;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 *
 * @author weibo
 */
public class TrainingCMD {

    private String embl;
    private String taxo;
    private String out;
    private String taxonomy="overall";
    private int level=1;
    private int subset=0;


    /**
     * 
     * @param args
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {

//        jmipara.PathSet.setLibDir(jmipara.PathSet.getPackageDir());//set library path
        //in NetBeans
//        jmipara.PathSet.setLibDir("/home/weibo/NetBeansProjects/JmiPara/lib");

        TrainingCMD tc=new TrainingCMD();

        tc.readArgv(args);

        MiRBaseData mi=new MiRBaseData();
        mi.loadMiRBaseData(tc.embl, tc.taxo);

        TrainingData td=new TrainingData();
        td.setLevel(tc.level);
        
        td.constructMatrix(mi.getDataOf(tc.taxonomy),tc.subset,tc.taxonomy);
        //svm_problem data=td.getTrainSet();

        td.saveTrainingData(tc.out);
//        td.saveFasta(tc.out);

    }

    private static void exit_with_help(){
		System.out.print(
		 "Usage: java -jar JmiParaTrain [options] embl_file taxo_file [out_file]\n\n"
		+"options:\n"
		+"\t-t taxonomy:  the taxonamy of miRNAs, can be animal, plant, \n"
                +"\t              virus or overall (default overall)\n"
		+"\t-l level:     the value of negative_data/positive_data, should\n"
                +"\t              be greater than 0 (default 1)\n"
                +"\t-s subset:    random draw s samples as a subset\n"
                +"\t-h help\n\n"+
                 "embl_file: the file includes all miRNAs of miRBase stored in EMBL\n"
                +"           format, which can download from miRBase\n"
                +"taxo_file: the file includes all organisms of miRNAs in miRBase,\n"
                +"           which can download from miRBase\n"+
                 "out_file:  output the SVM training dataset\n"       
		);
		System.exit(1);
    }

    private void readArgv(String[] argv){
        int i;
        for(i=0;i<argv.length;i++){
            if (argv[i].charAt(0) != '-')	break;
            if(++i>=argv.length)  exit_with_help();
            switch(argv[i-1].charAt(1)){
                case 't': taxonomy = argv[i];	break;
		case 'l': level = Integer.parseInt(argv[i]);	break;
                case 's': subset = Integer.parseInt(argv[i]);  break;
                case 'h': exit_with_help();
		default:
			System.err.println("unknown option");
			exit_with_help();
            }
	}
        if(!(taxonomy.equals("animal") || taxonomy.equals("plant")
                || taxonomy.equals("virus") || taxonomy.equals("overall"))){
            System.err.println("unknow taxonomy!");
            exit_with_help();
        }
        if(level<0){
            System.err.println("level should be greater than 0");
            exit_with_help();
        }
        if(subset<0){
            System.err.println("subset number should not be less than 0");
        }
        // determine filenames, at least two filenames
        if(i>=argv.length-1){
            System.err.println("lack input file!");
            exit_with_help();
        }
        embl=argv[i];
        taxo=argv[i+1];
        if(i<argv.length-2)
            out=argv[i+2];
        else{

            out=taxonomy+"-"+level+".td";
        }


    }
}
