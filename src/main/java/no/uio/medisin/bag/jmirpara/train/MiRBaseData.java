

package no.uio.medisin.bag.jmirpara.train;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * this class loads the data from miRBase
 * @author weibo
 */
public class MiRBaseData {

    
    private HashMap taxo=new HashMap();
    private ArrayList<MiRBaseEntry> miList=new ArrayList<MiRBaseEntry>();


    /**
     * load the data from miRBase
     * @param emblfile
     * @param taxofile
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void loadMiRBaseData(String emblfile,String taxofile) throws FileNotFoundException, IOException{
        parseEMBL(emblfile);
        parseTaxo(taxofile);
    }

    
    /**
     * get the data of specified class: animal,plant,virus or overall
     * @param taxo
     * @return
     */
    public ArrayList<MiRBaseEntry> getDataOf(String taxo){
        ArrayList<MiRBaseEntry> mis=new ArrayList<MiRBaseEntry>();
        for(MiRBaseEntry mi : miList){
            if(taxo.equals("animal")){
                if(mi.getTaxonomy().contains("Metazoa"))
                    mis.add(mi);
            }
            else if(taxo.equals("plant")){
                if(mi.getTaxonomy().contains("Viridiplantae"))
                    mis.add(mi);
            }
            else if(taxo.equals("virus")){
                if(mi.getTaxonomy().contains("Viruses"))
                    mis.add(mi);
            }
            else mis.add(mi); // overall
        }
        return mis;
    }
    
    /**
     * parse EMBL file of mirna data from miRBase database
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void parseEMBL(String emblfile) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(emblfile));
        String line = "";
        line = br.readLine();
        while(line!=null){
            MiRBaseEntry miEntry=new MiRBaseEntry();
            String seq="";
            while(line.substring(0,2).equals("//") == false){
                
                if(line.startsWith("ID"))
                    miEntry.setPriID(line.substring(5, line.indexOf(" ", 5)));
                else if(line.startsWith("AC"))
                    miEntry.setPriAccession(line.substring(line.lastIndexOf(" ")+1,line.length()-1));
                else if(line.contains("FT   miRNA")){
                    miEntry.addMat();
                    String[] p=line.substring(line.lastIndexOf(" ")+1).split("\\.\\.");
                    miEntry.setMatStart(Integer.parseInt(p[0]));
                    miEntry.setMatEnd(Integer.parseInt(p[1]));
                }
                else if(line.contains("/accession="))
                    miEntry.setMatAcc(
                            line.substring(line.indexOf("\"")+1, line.lastIndexOf("\"")));
                else if(line.contains("/product="))
                    miEntry.setMatID(
                            line.substring(line.indexOf("\"")+1, line.lastIndexOf("\"")));
                else if(line.contains("/evidence=experimental"))
                    miEntry.setEvidence(line.substring(line.indexOf("=")+1));
                else if(line.startsWith("SQ"))
                    seq="";
                else if(line.startsWith(" "))
                    seq+=line.replaceAll("\\s+|\\d+", "");
                line=br.readLine();
            }
            miEntry.setPriSequence(seq);
            miList.add(miEntry);
            line=br.readLine();
        }
        br.close();
        //System.out.println(miList.size());
    }

    /**
     * add taxonomy information to entry
     * @param taxofile
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void parseTaxo(String taxofile) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(taxofile));
        String line = "";
        while((line=br.readLine())!=null){
            String[] en=line.split("\t");
            taxo.put(en[0].trim(), new String[]{en[2].trim(),en[3].trim()});
        }
        br.close();

        for(MiRBaseEntry mi : miList){
            String abbr=mi.getPriID().substring(0, mi.getPriID().indexOf("-"));
            String[] ta=(String[])taxo.get(abbr);
            mi.setOrgnism(ta[0]);
            mi.setTaxonomy(ta[1]);
        }

    }




}
