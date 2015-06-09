/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;

/**
 *
 * @author wb
 */
public class MirViewData {
    public String fname;
    public String name;
    public int maxy;
    public int maxx;
    public ArrayList<int[]> pris=new ArrayList<int[]>();
    public ArrayList<int[]> mirs=new ArrayList<int[]>();
    public ArrayList<int[]> hits=new ArrayList<int[]>();
    public ArrayList<String> hitNames=new ArrayList<String>();
    public ArrayList<int[]> index=new ArrayList<int[]>();
    public int[][] dis=new int[100][3];
    public int maxbin;
    public RandomAccessFile rf;
    
    public MirViewData(File file) throws FileNotFoundException, IOException{
        this.fname=getBaseFileName(file.getAbsolutePath());
        parseTabResults();
        parseIndResults();
        rf = new RandomAccessFile(fname+".txt", "r");
    }
    
    public void parseIndResults() throws FileNotFoundException, IOException{
        BufferedReader br=new BufferedReader(new FileReader(fname+".ind"));
        String line="";
        while((line=br.readLine())!=null){
            String[] item=line.split("\t");
            int i=Integer.parseInt(item[0]);
            int s=Integer.parseInt(item[1]);
            int e=Integer.parseInt(item[2]);
            index.add(new int[]{i,s,e});
        }
        br.close();
    }
    
    public void parseTabResults() throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(fname+".tab"));
        name=br.readLine();
        maxx=Integer.parseInt(br.readLine());
        br.readLine();
        String line = "";
        String pri="";
        int mirCount=0;
        while((line = br.readLine()) != null){
            
            String[] item=line.split("\t",11);
            int priStart=Integer.parseInt(item[1]);
            int priEnd=Integer.parseInt(item[2]);
            if(item[0].equals(pri)==false){
                maxy=Math.max(maxy, mirCount);
                pris.add(new int[]{priStart,priEnd});
                mirCount=1;
            }
            pri=item[0];
            int mirStart=priStart+Integer.parseInt(item[6])-1;
            int mirEnd=mirStart+Integer.parseInt(item[7])-1;
            int miry=mirCount;
            int[] mi=new int[]{mirStart,mirEnd,miry};
            mirs.add(mi);
            String hitsName=item[10];
            if(hitsName.equals("")==false){
                hits.add(mi);
                hitNames.add(hitsName);
            }
            
            mirCount+=1;
        }
        br.close();
        maxy=Math.max(maxy, mirCount);
        
        overall();
    }
    
    public static String summary(String fname) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(fname+".tab"));
        String name=br.readLine();
        int seqLength=Integer.parseInt(br.readLine());
        br.readLine();
        int mirCount=0;
        while(br.readLine() != null){
            mirCount++;
        }
        return name+"\t"+seqLength+"\t"+mirCount;
    }
    
    public void overall(){
        maxbin=0;
        int bin=maxx/100;
        int p=0;
        for(int i=0;i<100;i++){
            int s=i*bin;
            int e=(i+1)*bin;
            dis[i][0]=s;
            dis[i][1]=e;
            dis[i][2]=0;
            for(;p<mirs.size();p++){
                int ms=mirs.get(p)[0];
                if(ms>s && ms<=e){
                    dis[i][2]+=1;
                }
                else if(ms>e){
                    break;
                }
            }
            if(maxbin<dis[i][2]){
                maxbin=dis[i][2];
            }
        }
    }
    
    
    public String getDetailInfo(int m) throws IOException{
        int[] mi=index.get(m);
        int n=mi[2]-mi[1];       
        byte[] binfo=new byte[n];
        rf.seek(mi[1]);
        for(int i=0;i<n;i++){  
            binfo[i]=rf.readByte();  
        }  
        String info=new String(binfo);
        return info;
    }
    
    public static String getBaseFileName(String filename) {   
        if ((filename != null) && (filename.length() > 0)) {   
            int dot = filename.lastIndexOf('.');   
            if ((dot >-1) && (dot < (filename.length()))) {   
                return filename.substring(0, dot);   
            }   
        }   
        return filename;   
    }
    
    
}
