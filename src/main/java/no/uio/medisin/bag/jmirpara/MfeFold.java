
package no.uio.medisin.bag.jmirpara;

/**
 *
 * @author weibo
 */
public class MfeFold {

    private static native float fold(String seq,float t);
    private static native void initIDs();
    static{
        System.loadLibrary("RNAFold");
	initIDs();
    }

    private static float temperature=37;
    private static String structure="";


    /**
     * fold the sequence
     */
    public static float cal(String sequence){
        structure="";
        return fold(sequence,temperature);
    }
    public static float cal(String sequence,String str){
        structure=str;
        return fold(sequence,temperature);
    }

    public static String getStructure(){
        return structure;
    }

}
