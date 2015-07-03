
package no.uio.medisin.bag.jmirpara;

import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author weibo
 */
public class FeatureRange {

//
//    Wu's range used in the Perl version of miRPara
//
//    private static final HashMap range=new HashMap();
//    static{
//        range.put("miRNA_size", new double[]{17,27});
//        range.put("preRNA_size", new double[]{41,250.96});
//        range.put("priRNA_size", new double[]{49,379});
//        range.put("basalSegment_size", new double[]{0,24});
//        range.put("lowerStem_size", new double[]{-1,116.42});
//        range.put("upperStem_size", new double[]{17,35});
//        range.put("topStem_size", new double[]{-1,103.85});
//        range.put("terminalLoop_size", new double[]{3,27.42});
//        range.put("miRNA_pair_number", new double[]{10,23});
//        range.put("preRNA_pair_number", new double[]{12,100.42});
//        range.put("priRNA_pair_number", new double[]{15,150});
//        range.put("preRNA_energy", new double[]{-172.2,-6.16});
//        range.put("priRNA_energy", new double[]{-307.7,-20});
//        range.put("miRNA_GC_content", new double[]{0.16,0.86});
//        range.put("priRNA_GC_content", new double[]{0.21,0.83});
//        range.put("miRNA_A_content", new double[]{0,0.57});
//        range.put("miRNA_C_content", new double[]{0,0.59});
//        range.put("miRNA_G_content", new double[]{0,0.7});
//        range.put("miRNA_U_content", new double[]{0,0.6});
//        range.put("preRNA_A_content", new double[]{0.05,0.43});
//        range.put("preRNA_C_content", new double[]{0.06,0.44});
//        range.put("preRNA_G_content", new double[]{0.08,0.49});
//        range.put("preRNA_U_content", new double[]{0.07,0.47});
//        range.put("priRNA_A_content", new double[]{0.05,0.42});
//        range.put("priRNA_C_content", new double[]{0.07,0.43});
//        range.put("priRNA_G_content", new double[]{0.1,0.5});
//        range.put("priRNA_U_content", new double[]{0.08,0.45});
//        range.put("miRNA_internalLoop_size", new double[]{0,12});
//        range.put("preRNA_internalLoop_size", new double[]{0,17.42});
//        range.put("lowerStem_internalLoop_size", new double[]{-1,16});
//        range.put("topStem_internalLoop_size", new double[]{-1,17.42});
//        range.put("miRNA_internalLoop_number", new double[]{0,6});
//        range.put("priRNA_internalLoop_number", new double[]{0,20.85});
//        range.put("lowerStem_internalLoop_number", new double[]{-1,13.42});
//        range.put("topStem_internalLoop_number", new double[]{-1,12});
//        range.put("miRNA_unpair_number", new double[]{0,26});
//        range.put("preRNA_unpair_number", new double[]{0,57.42});
//        range.put("priRNA_unpair_number", new double[]{0,80.92});
//        range.put("lowerstem_unpair_number", new double[]{-1,56.54});
//        range.put("topStem_unpair_number", new double[]{-1,52.85});
//        range.put("miRNA_unpair_rate", new double[]{0,0.51});
//        range.put("preRNA_unpair_rate", new double[]{0,0.47});
//        range.put("priRNA_unpair_rate", new double[]{0,0.46});
//        range.put("lowerStem_unpair_rate", new double[]{-1,0.64});
//        range.put("topStem_unpair_rate", new double[]{-1,0.75});
//        range.put("miRNA_G-U_number", new double[]{0,6});
//        range.put("miRNA_start", new double[]{1,258.38});
//        range.put("miRNA_end", new double[]{20,279.65});
//        range.put("upperStem_start", new double[]{1,117.42});
//        range.put("upperStem_end", new double[]{19,136.42});
//        range.put("miRNA_stability", new double[]{-1,5.21});
//        range.put("preRNA_internalLoop_number", new double[]{-1,14});
//        range.put("preRNA_GC_content", new double[]{0.2,0.85});
//        range.put("priRNA_internalLoop_size", new double[]{0,23.85});
//        range.put("preRNA_G-U_number", new double[]{0,9});
//        range.put("priRNA_G-U_number", new double[]{0,12});
//
//
//    }

    private static final HashMap range_overall      = new HashMap();
    private static final HashMap range_animal       = new HashMap();
    private static final HashMap range_plant        = new HashMap();
    private static final HashMap range_virus        = new HashMap();

    static{
	range_animal.put("priRNA_size", new double[]{57,140});
	range_animal.put("preRNA_size", new double[]{47,96.25});
	range_animal.put("miRNA_size", new double[]{18,25});
	range_animal.put("miRNA_start", new double[]{1,81});
	range_animal.put("miRNA_end", new double[]{22,101.25});
	range_animal.put("upperStem_start", new double[]{1,45});
	range_animal.put("upperStem_end", new double[]{22,67});
	range_animal.put("basalSegment_size", new double[]{0,17});
	range_animal.put("lowerStem_size", new double[]{0,42});
	range_animal.put("upperStem_size", new double[]{19,30});
	range_animal.put("topStem_size", new double[]{0,26});
	range_animal.put("terminalLoop_size", new double[]{3,17});
	range_animal.put("priRNA_energy", new double[]{-66.5,-17});
	range_animal.put("preRNA_energy", new double[]{-44.15,-10.875});
	range_animal.put("miRNA_stability", new double[]{0,4});
	range_animal.put("priRNA_GC_content", new double[]{0.28571427,0.750771613});
	range_animal.put("preRNA_GC_content", new double[]{0.2631579,0.75347428});
	range_animal.put("miRNA_GC_content", new double[]{0.26086956,0.77272725});
	range_animal.put("priRNA_A_content", new double[]{0.090617702,0.36986303});
	range_animal.put("priRNA_U_content", new double[]{0.13142228,0.418426665});
	range_animal.put("priRNA_G_content", new double[]{0.14627294,0.414699035});
	range_animal.put("priRNA_C_content", new double[]{0.111111104,0.38181818});
	range_animal.put("preRNA_A_content", new double[]{0.079260648,0.38095236});
	range_animal.put("preRNA_U_content", new double[]{0.125,0.4285714});
	range_animal.put("preRNA_G_content", new double[]{0.13109876,0.421258223});
	range_animal.put("preRNA_C_content", new double[]{0.104220256,0.38095236});
	range_animal.put("miRNA_A_content", new double[]{0,0.47619045});
	range_animal.put("miRNA_U_content", new double[]{0.047619045,0.5217391});
	range_animal.put("miRNA_G_content", new double[]{0.04545456,0.529352243});
	range_animal.put("miRNA_C_content", new double[]{0,0.5});
	range_animal.put("priRNA_pair_number", new double[]{18,49});
	range_animal.put("preRNA_pair_number", new double[]{16,33});
	range_animal.put("miRNA_pair_number", new double[]{13,22});
	range_animal.put("priRNA_unpair_number", new double[]{3,38});
	range_animal.put("priRNA_unpair_rate", new double[]{0.052631579,0.344000012});
	range_animal.put("preRNA_unpair_number", new double[]{2,27});
	range_animal.put("preRNA_unpair_rate", new double[]{0.043478262,0.358490556});
	range_animal.put("miRNA_unpair_number", new double[]{1,17});
	range_animal.put("miRNA_unpair_rate", new double[]{0.023255814,0.369453043});
	range_animal.put("lowerStem_unpair_number", new double[]{0,26});
	range_animal.put("lowerStem_unpair_rate", new double[]{0,0.5});
	range_animal.put("topStem_unpair_number", new double[]{0,15.25});
	range_animal.put("topStem_unpair_rate", new double[]{0,0.666666687});
	range_animal.put("priRNA_G-U_number", new double[]{0,10});
	range_animal.put("preRNA_G-U_number", new double[]{0,7});
	range_animal.put("miRNA_G-U_number", new double[]{0,5});
	range_animal.put("priRNA_internalLoop_number", new double[]{2,10});
	range_animal.put("priRNA_internalLoop_size", new double[]{1,12});
	range_animal.put("preRNA_internalLoop_number", new double[]{1,7});
	range_animal.put("preRNA_internalLoop_size", new double[]{1,10});
	range_animal.put("miRNA_internalLoop_number", new double[]{0,5});
	range_animal.put("miRNA_internalLoop_size", new double[]{0,6});
	range_animal.put("lowerStem_internalLoop_number", new double[]{0,6});
	range_animal.put("lowerStem_internalLoop_size", new double[]{0,11});
	range_animal.put("topStem_internalLoop_number", new double[]{0,4});
	range_animal.put("topStem_internalLoop_size", new double[]{0,8});
    }
    static{
	range_overall.put("priRNA_size", new double[]{57,260.02});
	range_overall.put("preRNA_size", new double[]{48,183});
	range_overall.put("miRNA_size", new double[]{18,25});
	range_overall.put("miRNA_start", new double[]{1,169});
	range_overall.put("miRNA_end", new double[]{22,189});
	range_overall.put("upperStem_start", new double[]{1,58});
	range_overall.put("upperStem_end", new double[]{22,79});
	range_overall.put("basalSegment_size", new double[]{0,16});
	range_overall.put("lowerStem_size", new double[]{0,53});
	range_overall.put("upperStem_size", new double[]{19,29});
	range_overall.put("topStem_size", new double[]{0,73});
	range_overall.put("terminalLoop_size", new double[]{3,18});
	range_overall.put("priRNA_energy", new double[]{-142.704,-17.296});
	range_overall.put("preRNA_energy", new double[]{-100.83,-10.598});
	range_overall.put("miRNA_stability", new double[]{0,3.67});
	range_overall.put("priRNA_GC_content", new double[]{0.2745098,0.745105131});
	range_overall.put("preRNA_GC_content", new double[]{0.25,0.746733337});
	range_overall.put("miRNA_GC_content", new double[]{0.25,0.77272725});
	range_overall.put("priRNA_A_content", new double[]{0.093207795,0.375024753});
	range_overall.put("priRNA_U_content", new double[]{0.13235295,0.42028987});
	range_overall.put("priRNA_G_content", new double[]{0.133319725,0.41250002});
	range_overall.put("priRNA_C_content", new double[]{0.111090717,0.378906404});
	range_overall.put("preRNA_A_content", new double[]{0.089277793,0.385972904});
	range_overall.put("preRNA_U_content", new double[]{0.128368962,0.4310345});
	range_overall.put("preRNA_G_content", new double[]{0.122444445,0.417730707});
	range_overall.put("preRNA_C_content", new double[]{0.10169494,0.38095236});
	range_overall.put("miRNA_A_content", new double[]{0.041633354,0.47619045});
	range_overall.put("miRNA_U_content", new double[]{0.047619045,0.5217391});
	range_overall.put("miRNA_G_content", new double[]{0.04545456,0.52380955});
	range_overall.put("miRNA_C_content", new double[]{0,0.5});
	range_overall.put("priRNA_pair_number", new double[]{19,103});
	range_overall.put("preRNA_pair_number", new double[]{16,71});
	range_overall.put("miRNA_pair_number", new double[]{13,22});
	range_overall.put("priRNA_unpair_number", new double[]{3,55});
	range_overall.put("priRNA_unpair_rate", new double[]{0.043478262,0.352007617});
	range_overall.put("preRNA_unpair_number", new double[]{2,41.02});
	range_overall.put("preRNA_unpair_rate", new double[]{0.03846154,0.363657765});
	range_overall.put("miRNA_unpair_number", new double[]{0,17});
	range_overall.put("miRNA_unpair_rate", new double[]{0,0.372549027});
	range_overall.put("lowerStem_unpair_number", new double[]{0,28});
	range_overall.put("lowerStem_unpair_rate", new double[]{0,0.5});
	range_overall.put("topStem_unpair_number", new double[]{0,37});
	range_overall.put("topStem_unpair_rate", new double[]{0,0.666666687});
	range_overall.put("priRNA_G-U_number", new double[]{0,11});
	range_overall.put("preRNA_G-U_number", new double[]{0,8});
	range_overall.put("miRNA_G-U_number", new double[]{0,5});
	range_overall.put("priRNA_internalLoop_number", new double[]{2,15});
	range_overall.put("priRNA_internalLoop_size", new double[]{1,14});
	range_overall.put("preRNA_internalLoop_number", new double[]{1,12});
	range_overall.put("preRNA_internalLoop_size", new double[]{1,12});
	range_overall.put("miRNA_internalLoop_number", new double[]{0,5});
	range_overall.put("miRNA_internalLoop_size", new double[]{0,6});
	range_overall.put("lowerStem_internalLoop_number", new double[]{0,7});
	range_overall.put("lowerStem_internalLoop_size", new double[]{0,11});
	range_overall.put("topStem_internalLoop_number", new double[]{0,10});
	range_overall.put("topStem_internalLoop_size", new double[]{0,11});
    }
    static{
	range_plant.put("priRNA_size", new double[]{69.04,398});
	range_plant.put("preRNA_size", new double[]{48,261.96});
	range_plant.put("miRNA_size", new double[]{19,24});
	range_plant.put("miRNA_start", new double[]{3,262.84});
	range_plant.put("miRNA_end", new double[]{23,285.78});
	range_plant.put("upperStem_start", new double[]{1.02,110});
	range_plant.put("upperStem_end", new double[]{21,130});
	range_plant.put("basalSegment_size", new double[]{0,16});
	range_plant.put("lowerStem_size", new double[]{0,106.86});
	range_plant.put("upperStem_size", new double[]{19,26});
	range_plant.put("topStem_size", new double[]{0.02,113});
	range_plant.put("terminalLoop_size", new double[]{3,21});
	range_plant.put("priRNA_energy", new double[]{-330.926,-18.261});
	range_plant.put("preRNA_energy", new double[]{-182.468,-9.155});
	range_plant.put("miRNA_stability", new double[]{2.50E-01,2.99E+00});
	range_plant.put("priRNA_GC_content", new double[]{0.243265574,0.723647035});
	range_plant.put("preRNA_GC_content", new double[]{0.222263747,0.734615406});
	range_plant.put("miRNA_GC_content", new double[]{0.218,0.762});
	range_plant.put("priRNA_A_content", new double[]{0.128434999,0.393496951});
	range_plant.put("priRNA_U_content", new double[]{0.1339286,0.433030113});
	range_plant.put("priRNA_G_content", new double[]{0.105721119,0.382808151});
	range_plant.put("priRNA_C_content", new double[]{0.10344827,0.3660714});
	range_plant.put("preRNA_A_content", new double[]{0.111196575,0.410519673});
	range_plant.put("preRNA_U_content", new double[]{0.130986374,0.443729981});
	range_plant.put("preRNA_G_content", new double[]{0.095254232,0.384455521});
	range_plant.put("preRNA_C_content", new double[]{0.085408513,0.378310832});
	range_plant.put("miRNA_A_content", new double[]{0.04545456,0.499565217});
	range_plant.put("miRNA_U_content", new double[]{0.047619045,0.52380955});
	range_plant.put("miRNA_G_content", new double[]{0.041666687,0.52380955});
	range_plant.put("miRNA_C_content", new double[]{0,0.47619045});
	range_plant.put("priRNA_pair_number", new double[]{23.02,171.9});
	range_plant.put("preRNA_pair_number", new double[]{17,101});
	range_plant.put("miRNA_pair_number", new double[]{14,23});
	range_plant.put("priRNA_unpair_number", new double[]{2.02,75});
	range_plant.put("priRNA_unpair_rate", new double[]{0.025683761,0.374889163});
	range_plant.put("preRNA_unpair_number", new double[]{0,62.98});
	range_plant.put("preRNA_unpair_rate", new double[]{0,0.401745434});
	range_plant.put("miRNA_unpair_number", new double[]{0,19});
	range_plant.put("miRNA_unpair_rate", new double[]{0,0.41641843});
	range_plant.put("lowerStem_unpair_number", new double[]{0,45.94});
	range_plant.put("lowerStem_unpair_rate", new double[]{0,0.490768337});
	range_plant.put("topStem_unpair_number", new double[]{0,52});
	range_plant.put("topStem_unpair_rate", new double[]{0,0.555412213});
	range_plant.put("priRNA_G-U_number", new double[]{0,14});
	range_plant.put("preRNA_G-U_number", new double[]{0,11.98});
	range_plant.put("miRNA_G-U_number", new double[]{0,5});
	range_plant.put("priRNA_internalLoop_number", new double[]{1.02,21});
	range_plant.put("priRNA_internalLoop_size", new double[]{1,20});
	range_plant.put("preRNA_internalLoop_number", new double[]{0,16});
	range_plant.put("preRNA_internalLoop_size", new double[]{0,19.98});
	range_plant.put("miRNA_internalLoop_number", new double[]{0,4});
	range_plant.put("miRNA_internalLoop_size", new double[]{0,6.96});
	range_plant.put("lowerStem_internalLoop_number", new double[]{0,12.98});
	range_plant.put("lowerStem_internalLoop_size", new double[]{0,12});
	range_plant.put("topStem_internalLoop_number", new double[]{0,13});
	range_plant.put("topStem_internalLoop_size", new double[]{0,17.98});
    }
    static{
	range_virus.put("priRNA_size", new double[]{60,116.06});
	range_virus.put("preRNA_size", new double[]{52,96.37});
	range_virus.put("miRNA_size", new double[]{18.21,24});
	range_virus.put("miRNA_start", new double[]{1.21,66.37});
	range_virus.put("miRNA_end", new double[]{22,87.58});
	range_virus.put("upperStem_start", new double[]{1,30.74});
	range_virus.put("upperStem_end", new double[]{22.21,51.74});
	range_virus.put("basalSegment_size", new double[]{0,9});
	range_virus.put("lowerStem_size", new double[]{0,29.74});
	range_virus.put("upperStem_size", new double[]{19.21,27.79});
	range_virus.put("topStem_size", new double[]{0,27.74});
	range_virus.put("terminalLoop_size", new double[]{3,14.79});
	range_virus.put("priRNA_energy", new double[]{-58.402,-21.205});
	range_virus.put("preRNA_energy", new double[]{-45.069,-16.821});
	range_virus.put("miRNA_stability", new double[]{0.00E+00,4.90E+00});
	range_virus.put("priRNA_GC_content", new double[]{0.371692288,0.744536513});
	range_virus.put("preRNA_GC_content", new double[]{0.37332627,0.73967692});
	range_virus.put("miRNA_GC_content", new double[]{0.33637678,0.81636362});
	range_virus.put("priRNA_A_content", new double[]{0.074672376,0.311774079});
	range_virus.put("priRNA_U_content", new double[]{0.122941538,0.356888625});
	range_virus.put("priRNA_G_content", new double[]{0.186714295,0.45454544});
	range_virus.put("priRNA_C_content", new double[]{0.15707691,0.392318524});
	range_virus.put("preRNA_A_content", new double[]{0.093578814,0.335296623});
	range_virus.put("preRNA_U_content", new double[]{0.119978807,0.390064963});
	range_virus.put("preRNA_G_content", new double[]{0.182136372,0.447605265});
	range_virus.put("preRNA_C_content", new double[]{0.159267792,0.39999998});
	range_virus.put("miRNA_A_content", new double[]{0.009545458,0.444999995});
	range_virus.put("miRNA_U_content", new double[]{0.042047115,0.5});
	range_virus.put("miRNA_G_content", new double[]{0.055119041,0.576818195});
	range_virus.put("miRNA_C_content", new double[]{0.045909102,0.561428552});
	range_virus.put("priRNA_pair_number", new double[]{19.21,42.37});
	range_virus.put("preRNA_pair_number", new double[]{17,34.58});
	range_virus.put("miRNA_pair_number", new double[]{13,22});
	range_virus.put("priRNA_unpair_number", new double[]{4,28.37});
	range_virus.put("priRNA_unpair_rate", new double[]{0.069271491,0.310507946});
	range_virus.put("preRNA_unpair_number", new double[]{3,22.79});
	range_virus.put("preRNA_unpair_rate", new double[]{0.059327731,0.338762498});
	range_virus.put("miRNA_unpair_number", new double[]{1,14});
	range_virus.put("miRNA_unpair_rate", new double[]{0.023255814,0.329918706});
	range_virus.put("lowerStem_unpair_number", new double[]{0,19});
	range_virus.put("lowerStem_unpair_rate", new double[]{0,0.48});
	range_virus.put("topStem_unpair_number", new double[]{0,15.16});
	range_virus.put("topStem_unpair_rate", new double[]{0,0.742500004});
	range_virus.put("priRNA_G-U_number", new double[]{0,9.37});
	range_virus.put("preRNA_G-U_number", new double[]{0,6});
	range_virus.put("miRNA_G-U_number", new double[]{0,5.79});
	range_virus.put("priRNA_internalLoop_number", new double[]{1.21,9});
	range_virus.put("priRNA_internalLoop_size", new double[]{1,9});
	range_virus.put("preRNA_internalLoop_number", new double[]{1.21,7.79});
	range_virus.put("preRNA_internalLoop_size", new double[]{1,9});
	range_virus.put("miRNA_internalLoop_number", new double[]{1,5});
	range_virus.put("miRNA_internalLoop_size", new double[]{1,5.58});
	range_virus.put("lowerStem_internalLoop_number", new double[]{0,4});
	range_virus.put("lowerStem_internalLoop_size", new double[]{0,8});
	range_virus.put("topStem_internalLoop_number", new double[]{0,4});
	range_virus.put("topStem_internalLoop_size", new double[]{0,8.37});
    }

    /**
     * check whether the parameter is within the defined range for the 
     * specified model
     * 
     * @param String para: the name of a parameter
     * @param double value: value of the parameter
     * @param model
     * 
     * @return boolean
     */
    private static boolean parameterInRange(String para, Object value, HashMap model){
        
        if(model.containsKey(para)){
            double v=Double.parseDouble(value.toString());
            double[] r=(double[])model.get(para);
            if(r[0]<=v && r[1]>=v)
                return true;
            else return false;
        }
        
        return true;
    }

    
    /**
     * Check whether 
     * @param feat
     * @param taxo
     * @return 
     */
    public static boolean featureInRange(HashMap feat, String taxo){
        Iterator paras=feat.keySet().iterator();
        Iterator values=feat.values().iterator();
        while(paras.hasNext()){
            String para=paras.next().toString();
//            double value=(Double)values.next();
            Object value=values.next();
            if(!parameterInRange(para, value, getRange(taxo))){
                return false;
            }
        }
        return true;
    }

    private static HashMap getRange(String taxo){
        HashMap range;

        if(taxo.equals("animal"))
            range=range_overall;
        else if(taxo.equals("plant"))
            range=range_plant;
        else if(taxo.equals("virus"))
            range=range_virus;
        else
            range=range_overall;

        return range;
    }

}
