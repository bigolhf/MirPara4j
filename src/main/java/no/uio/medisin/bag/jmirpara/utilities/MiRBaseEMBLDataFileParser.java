/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import no.uio.medisin.bag.core.MiRNAFeature;
import no.uio.medisin.bag.core.PreMiRNA;
import no.uio.medisin.bag.core.ShortPubMedEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author simonray
 */
public class MiRBaseEMBLDataFileParser {
    
    static Logger                       logger                      = LogManager.getLogger();
    
    private String                      emblFileName                = "";
    /*
     ID - identification             (begins each entry; 1 per entry)
     AC - accession number           (>=1 per entry)
     PR - project identifier         (0 or 1 per entry)
     DT - date                       (2 per entry)
     DE - description                (>=1 per entry)
     KW - keyword                    (>=1 per entry)
     OS - organism species           (>=1 per entry)
     OC - organism classification    (>=1 per entry)
     OG - organelle                  (0 or 1 per entry)
     RN - reference number           (>=1 per entry)
     RC - reference comment          (>=0 per entry)
     RP - reference positions        (>=1 per entry)
     RX - reference cross-reference  (>=0 per entry)
     RG - reference group            (>=0 per entry)
     RA - reference author(s)        (>=0 per entry)
     RT - reference title            (>=1 per entry)
     RL - reference location         (>=1 per entry)
     DR - database cross-reference   (>=0 per entry)
     CC - comments or notes          (>=0 per entry)
     AH - assembly header            (0 or 1 per entry)   
     AS - assembly information       (0 or >=1 per entry)
     FH - feature table header       (2 per entry)
     FT - feature table data         (>=2 per entry)    
     XX - spacer line                (many per entry)
     SQ - sequence header            (1 per entry)
     CO - contig/construct line      (0 or >=1 per entry) 
     bb - (blanks) sequence data     (>=1 per entry)
     // - termination line           (ends each entry; 1 per entry)
    
    */
    
    /*
    only the following subset of codes are used in release 21 of miRBase
    */
    private static final String     TERMINATION_LINE        =   "//";
    private static final String     ACCESSION_NUMBER        =   "AC";
    private static final String     COMMENT_LINE            =   "CC";
    private static final String     DESCRIPTION             =   "DE";
    private static final String     DB_XREF                 =   "DR";
    private static final String     FEATURE_TABLE_HEADER    =   "FH";
    private static final String     FEATURE_TABLE_DATA      =   "FT";
    private static final String     IDENTIFICATION          =   "ID";
    private static final String     REFERENCE_AUTHOR        =   "RA";
    private static final String     REFERENCE_COMMENT       =   "RC";
    private static final String     REFERENCE_LOCATION      =   "RL";
    private static final String     REFERENCE_NUMBER        =   "RN";
    private static final String     REFERENCE_TITLE         =   "RT";
    private static final String     REFERENCE_XF            =   "RX";
    private static final String     SEQUENCE_HEADER         =   "SQ";
    private static final String     SPACER_LINE             =   "XX";
    private static final String     EMPTY_LINE              =   "  ";
    
    private static final int        DATA_START_COLUMN       =   5;
    
    


    private ArrayList<PreMiRNA>     preMiRNAList            = new ArrayList<>();

    
    
    
    
    
    /**
     * 
     * @return
     * @throws IOException 
     */
    public int parseEMBLFile() throws IOException, Exception{
        
        String emblLine = "";
        MiRNAFeature miRNAFeature;
        PreMiRNA    preMiRNA = new PreMiRNA();
        ArrayList<String> referenceLines = new ArrayList<>();
        ArrayList<String> miRNALines = new ArrayList<>();
        ArrayList<String> dbxfLines = new ArrayList<>();
        ArrayList<String> commentLines = new ArrayList<>();
        ArrayList<String> seqLines = new ArrayList<>();
        try{
            BufferedReader brEMBL = new BufferedReader(new FileReader(new File(getEmblFileName())));
            while((emblLine = brEMBL.readLine())!=null){
                    switch (emblLine.trim().substring(0, 2)){
                        
                        case TERMINATION_LINE:
                            // add to pre-miRNA list
                            preMiRNAList.add(preMiRNA);
                            preMiRNA = new PreMiRNA();
                            break;
                        
                        //AC   MI0000001;
                        case ACCESSION_NUMBER:
                            String accNum = emblLine.split("\\s+")[1].trim();
                            accNum = accNum.substring(0,accNum.length()-1);
                            preMiRNA.setAccessionNumber(accNum);
                            break;
                            
                        case COMMENT_LINE:
                            commentLines.add(emblLine);
                            break;
                            
                        // DE   Caenorhabditis elegans let-7 stem-loop    
                        case DESCRIPTION:
                            String hostName = emblLine.split("\\s+")[1].trim() + " " + emblLine.split("\\s+")[2].trim();
                            preMiRNA.setHost(hostName);
                            break;
                            
                        case DB_XREF:
                            dbxfLines.add(emblLine);
                            break;
                            
                        // do nothing, the header line has no useful information    
                        case FEATURE_TABLE_HEADER: 
                            break;
                            
                        /*
                         Here we need to parse out miRNA entries. 
                            each miRNA starts with a FT line containing the key "miRNA", 
                            followed by additional lines specifying additional qualifiers.
                            the last miRNA entry is followed by a spacer line (XX) 
                        */    
                        case FEATURE_TABLE_DATA:
                            if(miRNALines.isEmpty()==false && emblLine.split("\\s+")[1].trim().equals("miRNA")==false 
                                    || emblLine.split("\\s+")[1].trim().equals("miRNA") && miRNALines.isEmpty()){
                                miRNALines.add(emblLine);
                            }else{
                                miRNAFeature = new MiRNAFeature();
                                miRNAFeature.parseMiRBaseMiRNAlines(miRNALines);
                                miRNALines.clear();
                                miRNALines.add(emblLine);
                                preMiRNA.addProduct(miRNAFeature);                               
                            }
                            
                            break;
                            
                        // ID   cel-let-7         standard; RNA; CEL; 99 BP    
                        case IDENTIFICATION:                            
                            preMiRNA.setName(emblLine.split("\\s+")[1].trim());
                            preMiRNA.setHost3Lettercode(emblLine.split("\\s+")[1].trim().split("-")[0].trim());
                            break;
                            
                        case REFERENCE_AUTHOR:
                            referenceLines.add(emblLine);
                            break;
                            
                        case REFERENCE_COMMENT:
                            referenceLines.add(emblLine);
                            break;
                            
                        case REFERENCE_LOCATION:
                            referenceLines.add(emblLine);
                            break;
                            
                        case REFERENCE_NUMBER:
                            referenceLines = new ArrayList<>();
                            referenceLines.add(emblLine);
                            break;
                            
                        case REFERENCE_TITLE:
                            referenceLines.add(emblLine);
                            break;
                            
                        case REFERENCE_XF:
                            referenceLines.add(emblLine);
                            break;
                            
                        case SEQUENCE_HEADER:
                        case EMPTY_LINE:
                            seqLines.add(emblLine);
                            break;
                            
                        /*
                            the action we take depends on how we got here
                            we can be ending a reference, a database x-ref, 
                            a comment line or an miRNA entry
                        */    
                        case SPACER_LINE: // 
                            if(referenceLines.size()>0){
                                ShortPubMedEntry shortPubMedEntry= new ShortPubMedEntry(referenceLines);
                                preMiRNA.addPubmedRef(shortPubMedEntry);
                                referenceLines.clear();
                                continue;
                            }
                            
                            if(miRNALines.size()>0){
                                miRNAFeature = new MiRNAFeature();
                                miRNAFeature.parseMiRBaseMiRNAlines(miRNALines);
                                preMiRNA.addProduct(miRNAFeature);
                                continue;
                            }
                            
                            if(dbxfLines.size()>0){
                                preMiRNA.setDbxrefs(parseDBXFlines(dbxfLines));
                                dbxfLines.clear();
                                continue;
                            }
                            
                            if(commentLines.size()>0){
                                preMiRNA.setNote(parseCommentlines(commentLines));
                                commentLines.clear();
                                continue;
                            }
                            
                            if(seqLines.size()>0){
                                parseSeqLines(seqLines);
                                seqLines.clear();
                            }
 
                            break;
                            
                        default:
                            break;
                    }
                }
            brEMBL.close();
        }
        catch(IOException exIO){
            logger.info("exception parsing miRBase EMBL file\nline<" 
                    + emblLine + ">\ncould not be parsed");
            logger.error("exception parsing miRBase EMBL file\nline<" 
                    + emblLine + ">\ncould not be parsed");
            throw new RuntimeException("exception parsing miRBase EMBL file\nline<" 
                    + emblLine + ">\ncould not be parsed");
        }
        
        return 0;
    }
    
    
    
    
    
    


    /**
     * concatenate the DB xRef lines
     * 
     * @param emblLines
     * @return 
     */
    private String parseDBXFlines(ArrayList<String> emblLines){
        String dbxfLine = "";
        for(String emblLine:emblLines){
            dbxfLine = dbxfLine.concat(emblLine.substring(DATA_START_COLUMN).trim() + "\t");
        }
        return dbxfLine.trim();
        
    }

    
    
    
    
    /**
     * concatenate the comment lines
     * 
     * @param emblLines
     * @return 
     */
    private String parseCommentlines(ArrayList<String> emblLines){
        String commentLine = "";
        for(String emblLine:emblLines){
            if(emblLine.substring(1, 2).equals(SEQUENCE_HEADER)) continue;
            commentLine = commentLine.concat(emblLine.substring(DATA_START_COLUMN).trim() + " ");
        }
        return commentLine.trim();
        
    }
    
    
    
    /**
     * parse out sequence from sequence lines
     * e.g.
     * SQ   Sequence 99 BP; 26 A; 19 C; 24 G; 0 T; 30 other;
     *      uacacugugg auccggugag guaguagguu guauaguuug gaauauuacc accggugaac        60
     *      uaugcaauuu ucuaccuuac cggagacaga acucuucga                               99
     * //
     * 
     * @param emblLines
     * @return 
     */
    private String parseSeqLines(ArrayList<String> emblLines){
        String seqLine = "";
        for(String emblLine:emblLines){
            seqLine = seqLine.concat(emblLine.substring(DATA_START_COLUMN, emblLine.length()-3).trim());
        }
        return seqLine.replaceAll("\\s+","");
        
    }
    
    
    /**
     * format of the miRBase ID line is
     *   ID   cel-let-7         standard; RNA; CEL; 99 BP.
     * 
     * the only useful information is the common name for the miRNA (cel-let-7 here)
     * Not certain what "standard" refers to 
     * 
     * @param line 
     */
    private String parseIDline(String line){
        
        return this.cutIDcode(line).split("\\\\s+")[0];
    }
    
    /**
     * remove ID code and leading spaces from input line
     * 
     * @param rawLine
     * @return 
     */
    private String cutIDcode(String rawLine){
        return rawLine.substring(2).trim();
    }

    /**
     * @return the emblFileName
     */
    public String getEmblFileName() {
        return emblFileName;
    }

    /**
     * @param emblFileName the emblFileName to set
     */
    public void setEmblFileName(String emblFileName) {
        this.emblFileName = emblFileName;
    }
}
