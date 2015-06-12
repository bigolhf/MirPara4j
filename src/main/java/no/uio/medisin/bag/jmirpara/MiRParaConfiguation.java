/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor. 
 */
package no.uio.medisin.bag.jmirpara;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import org.yaml.snakeyaml.Yaml;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *
 * @author sr
 */
public class MiRParaConfiguation {
    static Logger logger = LogManager.getRootLogger(); 
    static Yaml yaml = new Yaml();
    static String FileSeparator = System.getProperty("file.separator");
    
    private String configurationFile="";
    
    private String installationFolder="";
    private String modelFolder="";
    private String rnaFoldFolder="";
    private String dataFolder="";
    
    private String rnaFoldDll="";
    private String mirbaseDataFile="";
    
    private String pathToLibrary = "";
    private String pathToMirbaseData = "";
    
    private String pathToModelData = "";
    

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
    private ArrayList<String> results=new ArrayList<String>();
    

    public MiRParaConfiguation()
    {
    }
    
    public void ReadConfigurationFile(String configFile){
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
            "window size :\t" + window + "\n"
          + "step size   :\t" + step + "\n"
          + "start       :\t" + start + "\n"
          + "distance    :\t" + distance + "\n"
          + "cutoff      :\t" + cutoff + "\n"
          + "model       :\t" + model + "\n"
          + "level       :\t" + level;
        
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
     * @return the model
     */
    public String getModel() {
        return model;
    }

    /**
     * @param model the model to set
     */
    public void setModel(String model) {
        this.model = model;
    }

    /**
     * @return the pathToLibrary
     */
    public String getPathToLibrary() {
        return pathToLibrary;
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

    /**
     * @return the pathToMirbaseData
     */
    public String getPathToMirbaseData() {
        return pathToMirbaseData;
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
     * @return the distance
     */
    public int getDistance() {
        return distance;
    }

    /**
     * @param distance the distance to set
     */
    public void setDistance(int distance) {
        this.distance = distance;
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
        this.cutoff = cutoff;
    }
    
    
    
}
