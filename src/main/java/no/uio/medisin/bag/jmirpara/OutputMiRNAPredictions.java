

package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Serialize a set of prediction results. 
 * 
 * 
 * @author Wei Bo & Simon Rayner
 */
public class OutputMiRNAPredictions {
    
    static Logger logger = LogManager.getRootLogger();    

    private static HashMap<String,String> mirbase;
    public static String outputFilePrefix;
    public static int fpos;
    public static int lpos;

    /**
     * 
     * @param predictionList
     * @param seq
     * @param appendData
     * 
     * @throws IOException 
     */
    public static void serializePredictionSummary(ArrayList<HashMap> predictionList, SimpleSeq seq, boolean appendData) throws IOException{
        
        logger.info("serialize summary predictions");
        BufferedWriter bwOverallTabbedData;
        
        if(appendData==false){
            
            bwOverallTabbedData=new BufferedWriter(new FileWriter(outputFilePrefix + ".tab"));
            bwOverallTabbedData.write(seq.getId() + "\n");
            bwOverallTabbedData.write(seq.getLength() + "\n");
            bwOverallTabbedData.write("pri_id\tpri_start\tpri_end\tpri_seq\tpri_str\tmir_id\tmir_start\tmir_size\tmir_seq\tmir_strand\thit\n");
            
        }
        else{
            
            bwOverallTabbedData = new BufferedWriter(new FileWriter(outputFilePrefix + ".tab", true));
            
        }
        
        sortList(predictionList);
        for(HashMap prediction:predictionList){
            StringBuilder entry = new StringBuilder();
            entry.append(prediction.get("priRNA_id")).append("\t").append(prediction.get("priRNA_start")).append("\t");
            entry.append(prediction.get("priRNA_end")).append("\t").append(prediction.get("priRNA_sequence")).append("\t");
            entry.append(prediction.get("priRNA_structure")).append("\t").append(prediction.get("miRNA_id")).append("\t");
            entry.append(prediction.get("miRNA_start")).append("\t").append(prediction.get("miRNA_size")).append("\t");
            entry.append(prediction.get("miRNA_sequence")).append("\t").append(prediction.get("strand")).append("\t");
            String hit=blastMirBase(prediction.get("miRNA_sequence").toString());
            entry.append(hit).append("\n");
            bwOverallTabbedData.write(entry.toString());
        }
        bwOverallTabbedData.close();
    }
    

    /**
     * the distribution of miRNA on the sequence
     * @param predictionList
     * @param seq
     * @throws IOException
     */
    public static void serializePredictionsAsHTML(ArrayList<HashMap> predictionList, SimpleSeq seq) throws IOException{

        logger.info("serialize prediction details as HTML");
        
        int seqW=1500, seqH=10, x0=5, y0=50, mirH=5;
        BufferedWriter bwLocationDataHtml = new BufferedWriter(new FileWriter(outputFilePrefix + "." + seq.getId() + ".html"));

        bwLocationDataHtml.write("<html>\n"
                + "<head>\n"
                + "<title>JmiPara 2.0 miRNA results</title>\n"
                + "<style type=\"text/css\">\n"
                + ".seq{width:"+seqW+"px;height:"+seqH+"px;background:#069;position:absolute;}\n"
                + ".mir{height:"+mirH+"px;position:absolute;}\n"
                + "</style>\n"
                + "<script type=text/javascript>\n"
                + "document.write(\"<style type='text/css'>#Tag {display:block;font:12px Tahoma,Verdana;background-color:#FFC;border:1px #000 solid;padding:3px;position:absolute;z-index:1000;visibility:hidden}</style>\");\n"
                + "document.write(\"<tt id='Tag' style='filter:blendtrans(duration=.2) revealTrans(duration=.1,transition=12) alpha(opacity=90,enabled=1);-moz-opacity:0.9'></tt>\");\n"
                + "document.onmouseover=ShowTag;\n"
                + "function ShowTag(e){\n"
                + "	e=e || window.event;\n"
                + "	var sPop = null;\n"
                + "	if(e){\n"
                + "		o=e.target;\n"
                + "		MouseX=e.pageX?e.pageX:clientX + document.body.scrollLeft - document.body.clientLeft;\n"
                + "		MouseY=e.pageY?e.pageY:clientY + document.body.scrollTop  - document.body.clientTop;\n"
                + "	}\n"
                + " 	else{\n"
                + "		o=event.srcElement;\n"
                + "		MouseX=event.x;\n"
                + "		MouseY=event.y;\n"
                + "	}	"
                + " 	if(o.title){\n"
                + "		o.pop=o.title;\n"
                + "		o.title=\"\";\n"
                + "	}\n"
                + " 	if(o.pop){\n"
                + "		o.pop=o.pop.replace(/\\n/g,\"<br>\");\n"
                + "	}\n"
                + "	if(o.pop!=sPop){\n"
                + "		sPop=o.pop;\n"
                + "		if(sPop){\n"
                + "			obj=(document.all)? Tag : document.getElementById(\"Tag\");\n"
                + "			obj.innerHTML=o.pop;\n"
                + "			iebody=document.body;\n"
                + "			objWidth=obj.offsetWidth;objHeight=obj.offsetHeight;\n"
                + "			popLeftAdjust=(MouseX+12+objWidth>iebody.clientWidth)?(-objWidth-24):0;\n"
                + "			popTopAdjust=(MouseY+12+objHeight>iebody.clientHeight)?(-objHeight-24):0;\n"
                + "			obj.style.left=MouseX+12+iebody.scrollLeft+popLeftAdjust;\n"
                + "			obj.style.top=MouseY+12+iebody.scrollTop+popTopAdjust;\n"
                + "			if(obj.filters && obj.filters.length!=0){\n"
                + "				obj.filters[1].apply();\n"
                + "				obj.style.visibility=\"visible\";\n"
                + "				obj.filters[1].play();\n"
                + "			}\n"
                + " 			else obj.style.visibility=\"visible\";\n"
                + "		}\n"
                + "		else{\n"
                + "			if(obj.filters && obj.filters.length!=0){\n"
                + "				obj.filters[0].apply();\n"
                + "				obj.style.visibility=\"hidden\";\n"
                + "				obj.filters[0].play();\n"
                + "			}\n"
                + " 			else obj.style.visibility=\"hidden\";\n"
                + "		}	}}\n"
                + "</script>\n"
                + "</head>\n"
                + "</body>\n"
                + "<h2>Distribution of miRNAs on sequence "+seq.getId()+"</h2>\n"
                + "<div class=\"seq\" style=\"top:"+y0+"px;left:"+x0+"px\"></div>\n");

        //sort predicted miRNAs by start and end position
        sortList(predictionList);
        double left,top,mirW;
        int retop=0;
        String title;
        for(int i=0;i<predictionList.size();i++){
            HashMap mir=predictionList.get(i);
            left=(Integer)mir.get("miStart")*seqW*1.0/seq.getLength()+x0;
            mirW=(Integer)mir.get("miRNA_size")*seqW*1.0/seq.getLength();

            if(i>0 && (Integer)mir.get("miStart")>=(Integer)predictionList.get(i-1).get("miEnd"))
                retop=i;
            top=(i-retop)*(mirH+5)+y0;
            title=mir.get("miRNA_id").toString();

            double rate=(Double)mir.get("probability");
            String color=getColor(new Double[]{0.0,0.0,255.0},new Double[]{255.0,0.0,0.0},rate);

            String hit=blastMirBase(mir.get("miRNA_sequence").toString());
            if(hit.equals(""))
                bwLocationDataHtml.write("<a href=#><div title=\""+title+"\" class=\"mir\" style=\"top:"+top+"px;left:"+left+"px;width:"+mirW+"px;background:"+color+"\"></div></a>\n");
            else
                bwLocationDataHtml.write("<a href=#><div title=\""+title+"<p>hits to miRBase entris:<br>"+hit+"\" class=\"mir\" style=\"top:"+top+"px;left:"+left+"px;width:"+mirW+"px;background:#1d0\"></div></a>\n");
        }


        bwLocationDataHtml.write("</body></html>");
        bwLocationDataHtml.close();
    }

    
    
   /**
    * detailed report showing feature values of each predicted miRNA
    * 
    * @param miRNAPredictionResults
    * @param seq : current Query Sequence
    * @param appendRecord: append this data to the current file
    * 
    * @throws IOException
    */
    public static void serializePredictionDetails(ArrayList<HashMap> miRNAPredictionResults, SimpleSeq seq, boolean appendRecord) throws IOException{
        
        logger.info("serialize detailed predictions");
        
        BufferedWriter bwTxtReportFile;
        BufferedWriter bwIndReportFile;
        File txtReportFile = new File(outputFilePrefix + ".txt");
        String str = "";
        
        if(appendRecord == false){
            
            bwTxtReportFile = new BufferedWriter(new FileWriter(txtReportFile, appendRecord));
            str = "miRPara4j Prediction Report\n\n";
            bwTxtReportFile.write(str);
            bwIndReportFile = new BufferedWriter(new FileWriter(outputFilePrefix + ".ind", appendRecord));
            fpos = str.length();
            lpos = 0;
            
        }
        else{
            
            bwTxtReportFile = new BufferedWriter(new FileWriter(outputFilePrefix + ".txt", appendRecord));
            bwIndReportFile = new BufferedWriter(new FileWriter(outputFilePrefix + ".ind", appendRecord));
            
        }
        
        
        String[] paraList = ModelSet.model("all");
        sortList(miRNAPredictionResults);
        for(HashMap result : miRNAPredictionResults){
            bwIndReportFile.write(lpos+"\t"+(fpos+1));
            str = repeat("=",30)+"> " + result.get("miRNA_id") + " <" + repeat("=",30) + "\n\n";
            bwTxtReportFile.write(str);
            fpos+= str.length();
            str = result.get("plot").toString()+"\n\n";
            bwTxtReportFile.write(str);
            fpos+= str.length();
            str = "probability=" + result.get("probability") + "\n\n";
            bwTxtReportFile.write(str);
            fpos+= str.length();
            for(String para : paraList){
                if(result.containsKey(para)){
                    str = padStringWithSpaces(para,30)+"\t\t"+result.get(para)+"\n";
                    bwTxtReportFile.write(str);
                    fpos+= str.length();
                }
            }
            str = "\n";
            bwTxtReportFile.write(str);
            fpos+= str.length();
            lpos+= 1;
            bwIndReportFile.write("\t" + fpos + "\n");
        }
        bwTxtReportFile.close();
        bwIndReportFile.close();
    }

    /**
     * load miRBase database data
     * 
     * @param file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static void loadMirBaseData(File file) throws FileNotFoundException, IOException{

        logger.info("load miRBase data");
        
        
        mirbase=new HashMap<String,String>();

        BufferedReader mib=new BufferedReader(new FileReader(file));
        String line=mib.readLine();
        while(line!=null){
            String[] pair=line.split("\t");
            if(mirbase.containsKey(pair[0]))
                mirbase.put(pair[0], mirbase.get(pair[0])+","+pair[1]+"("+pair[2]+")");
            else
                mirbase.put(pair[0], pair[1]+"("+pair[2]+")");
            line=mib.readLine();
        }
        mib.close();

        Iterator iter = mirbase.entrySet().iterator();
        while (iter.hasNext()) {
            Map.Entry entry = (Map.Entry) iter.next();
            String key = (String)entry.getKey();
            String val = ((String)entry.getValue()).replaceAll(",", "<br>");
            
            mirbase.put(key, val);
        }
 

    }

    /**
     * compare the prediction results to the miRBase data
     * 
     * @param mir
     * @return
     */
    private static String blastMirBase(String mir){
        return mirbase.containsKey(mir)?mirbase.get(mir):"";
    }

    /**
     * sort miRNAs by start and end
     * @param list
     */
    public static void sortList(ArrayList list){
        Collections.sort(list, new Comparator<HashMap>(){
            public int compare(HashMap o1,HashMap o2){
                double pc=0;
                pc=(Integer)o1.get("priRNA_start")-(Integer)o2.get("priRNA_start");
                if(pc==0) {
                double c=(Integer)o1.get("miRNA_start")-(Integer)o2.get("miRNA_start");
                if(c>0) return 1;
                else if(c<0) return -1;
                else{
                    c=(Integer)o1.get("miRNA_end")-(Integer)o2.get("miRNA_end");
                    if(c>0) return 1;
                    else if(c<0) return -1;
                    else return 0;
                }
                }
                else if(pc<0) return -1;
                else return 1;
            }
        });
    }

    
    /**
     * make the string the same length by padding with spaces
     * @param str
     * @param size
     * @return
     */
    public static String padStringWithSpaces(String str, int size){
        String strf=str;
        int num=size-strf.length();
        if(num>0)
            for(int i=0;i<num;i++)
                strf+=" ";
        return strf;
    }

    
    
    
    /**
     * repeat a string
     * @param str
     * @param num
     * @return
     */
    public static String repeat(String str, int num){
        String strs="";
        for(int i=0;i<num;i++)
            strs+=str;
        return strs;
    }

    /**
     * backspace a string
     * @param con
     * @return
     */
    public static String backspace(String con){
        return repeat("\b",con.length());
    }
    /**
     * reserve the first twofigures after the decimal point
     * @param num
     * @return
     */
    public static String decimal(double num){
        String n=String.valueOf(num);
        if(n.length()-n.indexOf(".")-1>2)
            return n.substring(0,n.indexOf(".")+3);
        else return n;
    }
    public static String decimal(double num, int d){
        String n=String.valueOf(num);
        if(n.length()-n.indexOf(".")-1>d)
            return n.substring(0,n.indexOf(".")+d+1);
        else return n;
    }

    public static String getColor(Double[] start, Double[] end, double rate){
        double[] rgb=new double[3];
        String color="#";
        for(int i=0;i<3;i++){
            rgb[i]=rate*(end[i]-start[i])+start[i];
            String c=Integer.toHexString((int)rgb[i]);
            if(c.length()==1)
                c="0"+c;
            color+=c;
        }
        return color;
    }
}
