

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

/**
 * Seems to contain a set of miRNA predictions. 
 * Could probably do with a better name
 * 
 * @author weibo
 */
public class Output {

    private static HashMap<String,String> mirbase;
    public static String fname;
    public static int fpos;
    public static int lpos;

    public static void overall(ArrayList<HashMap> res, SimSeq seq, boolean a) throws IOException{
        
        BufferedWriter ov;
        if(a==false){
            ov=new BufferedWriter(new FileWriter(fname+".tab"));
            ov.write(seq.getId()+"\n");
            ov.write(seq.getLength()+"\n");
            ov.write("pri_id\tpri_start\tpri_end\tpri_seq\tpri_str\tmir_id\tmir_start\tmir_size\tmir_seq\tmir_strand\thit\n");
        }
        else{
            ov=new BufferedWriter(new FileWriter(fname+".tab",true));
        }
        sortList(res);
        for(HashMap re:res){
            StringBuilder entry=new StringBuilder();
            entry.append(re.get("priRNA_id")).append("\t").append(re.get("priRNA_start")).append("\t");
            entry.append(re.get("priRNA_end")).append("\t").append(re.get("priRNA_sequence")).append("\t");
            entry.append(re.get("priRNA_structure")).append("\t").append(re.get("miRNA_id")).append("\t");
            entry.append(re.get("miRNA_start")).append("\t").append(re.get("miRNA_size")).append("\t");
            entry.append(re.get("miRNA_sequence")).append("\t").append(re.get("strand")).append("\t");
            String hit=blastMirBase(re.get("miRNA_sequence").toString());
            entry.append(hit).append("\n");
            ov.write(entry.toString());
        }
        ov.close();
    }
    

    /**
     * the distribution of miRNA on the sequence
     * @param resultList
     * @param seq
     * @throws IOException
     */
    public static void distribution(ArrayList<HashMap> resultList, SimSeq seq) throws IOException{
        
        int seqW=1500,seqH=10,x0=5,y0=50,mirH=5;
        BufferedWriter di=new BufferedWriter(new FileWriter(fname+".html"));

        di.write("<html>\n"
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

        //sort result miRNAs by start and end
        sortList(resultList);
        double left,top,mirW;
        int retop=0;
        String title;
        for(int i=0;i<resultList.size();i++){
            HashMap mir=resultList.get(i);
            left=(Integer)mir.get("miStart")*seqW*1.0/seq.getLength()+x0;
            mirW=(Integer)mir.get("miRNA_size")*seqW*1.0/seq.getLength();

            if(i>0 && (Integer)mir.get("miStart")>=(Integer)resultList.get(i-1).get("miEnd"))
                retop=i;
            top=(i-retop)*(mirH+5)+y0;
            title=mir.get("miRNA_id").toString();

            double rate=(Double)mir.get("probability");
            String color=getColor(new Double[]{0.0,0.0,255.0},new Double[]{255.0,0.0,0.0},rate);

            String hit=blastMirBase(mir.get("miRNA_sequence").toString());
            if(hit.equals(""))
                di.write("<a href=#><div title=\""+title+"\" class=\"mir\" style=\"top:"+top+"px;left:"+left+"px;width:"+mirW+"px;background:"+color+"\"></div></a>\n");
            else
                di.write("<a href=#><div title=\""+title+"<p>hits to miRBase entris:<br>"+hit+"\" class=\"mir\" style=\"top:"+top+"px;left:"+left+"px;width:"+mirW+"px;background:#1d0\"></div></a>\n");
        }


        di.write("</body></html>");
        di.close();
    }

   /**
    * detail feature values of miRNAs
    * @param resultList
    * @throws IOException
    */
    public static void detail(ArrayList<HashMap> resultList, SimSeq seq, boolean a) throws IOException{
        
        BufferedWriter de;
        BufferedWriter in;
        File txt=new File(fname+".txt");
        String str="";
        if(a==false){
            de=new BufferedWriter(new FileWriter(txt));
            str="JmiPara-2.0 miRNA Scan Report\n\n";
            de.write(str);
            in=new BufferedWriter(new FileWriter(fname+".ind"));
            fpos=str.length();
            lpos=0;
        }
        else{
            de=new BufferedWriter(new FileWriter(fname+".txt",true));
            in=new BufferedWriter(new FileWriter(fname+".ind",true));
        }
        String[] paraList=ModelSet.model("all");
        sortList(resultList);
        for(HashMap result : resultList){
            in.write(lpos+"\t"+(fpos+1));
            str=repeat("=",30)+"> "+result.get("miRNA_id")+" <"+repeat("=",30)+"\n\n";
            de.write(str);
            fpos+=str.length();
            str=result.get("plot").toString()+"\n\n";
            de.write(str);
            fpos+=str.length();
            str="probability="+result.get("probability")+"\n\n";
            de.write(str);
            fpos+=str.length();
            for(String para : paraList){
                if(result.containsKey(para)){
                    str=flush(para,30)+"\t\t"+result.get(para)+"\n";
                    de.write(str);
                    fpos+=str.length();
                }
            }
            str="\n";
            de.write(str);
            fpos+=str.length();
            lpos+=1;
            in.write("\t"+fpos+"\n");
        }
        de.close();
        in.close();
    }

    /**
     * load miRBase database data
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static void loadMirBaseData(File file) throws FileNotFoundException, IOException{
        //load mirbase data
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
     * compare the predicting results to the miRBase data
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
     * make the string the same length
     * @param str
     * @param size
     * @return
     */
    public static String flush(String str, int size){
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
