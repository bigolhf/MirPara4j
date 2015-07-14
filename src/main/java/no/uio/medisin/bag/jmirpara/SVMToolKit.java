
package no.uio.medisin.bag.jmirpara;

import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;
import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_problem;

/**
 * about SVM
 * @author weibo
 */
public class SVMToolKit {

    private static svm_model svmmodel;
    private static double result;



    /**
     * stores feature values of specified model into double[]
     * @param model
     * @param feat
     * @return
     */
    public static double[] featValueOf(String model, HashMap feat){
        String[] mn=ModelSet.model(model);
        double[] fl=new double[mn.length];
        for(int i=0;i<fl.length;i++){
            if(mn[i].equals("miRNA_firstBase") || mn[i].contains("overhang"))
                fl[i]=(double)nt2Num((Character)feat.get(mn[i]));
            else
                fl[i]=Double.parseDouble(feat.get(mn[i]).toString());
        }
        return fl;
    }

    /**
     * build svm_node matrix from double[]
     * @param feats
     * @return
     */
    public static svm_node[] buildNodes(double[] feats){
        svm_node[] x = new svm_node[feats.length];
	for(int j=0;j<x.length;j++){
            x[j] = new svm_node();
            x[j].index = j+1;
            x[j].value = feats[j];
	}
        return x;
    }

    /**
     * load one entry in the svm format file
     * @param feats
     * @return
     */
    private static svm_node[] loadNodes(String feats, double[] re){
        StringTokenizer st = new StringTokenizer(feats," \t\n\r\f:");

	re[0]=Double.parseDouble(st.nextToken());
	int m = st.countTokens()/2;
	svm_node[] x = new svm_node[m];
	for(int j=0;j<m;j++){
            x[j] = new svm_node();
            x[j].index = Integer.parseInt(st.nextToken());
            x[j].value = Double.parseDouble(st.nextToken());
	}
        return x;
    }

    /**
     * build svm_node matrix with label
     * @param label
     * @param feats
     * @return
     */
    public static svm_node[] buildNodes(int label,double[] feats){
        svm_node[] x = new svm_node[feats.length+1];
	for(int j=1;j<x.length;j++){
            x[j] = new svm_node();
            x[j].index = j;
            //x[j].index=j-1;//index start from 0
            x[j].value = feats[j-1];
	}
        x[0]=new svm_node();
        x[0].index=0;
        x[0].value=label;
        return x;
    }

    /**
     * build svm_problem matrix
     * @param x
     * @return
     */
    public static svm_problem trainMatrix(svm_node[][] x){
        svm_problem data=new svm_problem();
        data.l=x.length;
        double[] y=new double[data.l];
        svm_node[][] xx=new svm_node[data.l][];
        for(int i=0;i<y.length;i++){
            y[i]=x[i][0].value;
            xx[i]=new svm_node[x[i].length-1];
            for(int j=0;j<x[i].length-1;j++)
                xx[i][j]=x[i][j+1];
        }
        data.y=y;
        data.x=xx;

        return data;
    }


    public static void loadModel(String model) throws IOException{
        svmmodel=svm.svm_load_model(model);
        svmmodel= svm.svm_load_model("/Users/simonray/NetBeansProjects/MirPara4j/models/aqu_1.model");
    }

    /**
     * svm predict
     * @param model
     * @param feat
     */
    public static void judge(String model, HashMap feat, double cutoff){
        svm_node[] data=buildNodes(featValueOf(model,feat));
//        result=svm.svm_predict(svmmodel, data);
        int nr_class=svm.svm_get_nr_class(svmmodel);
//        int[] labels=new int[nr_class];
//	svm.svm_get_labels(svmmodel,labels);
	double[] prob_estimates= new double[nr_class];
        result = svm.svm_predict_probability(svmmodel,data,prob_estimates);
        if(prob_estimates[0]>cutoff)
            result=1;
        else
            result=0;
        feat.put("probability", prob_estimates[0]);
    }

    public static double[] predict_test(String model, String feat){
        double[] results=new double[4];
        svm_node[] data=loadNodes(feat, results);
        double[] dec=new double[1];
        results[1]=svm.svm_predict_values(svmmodel, data, dec);
        results[2]=dec[0];
        double[] prob=new double[2];
        svm.svm_predict_probability(svmmodel, data, prob);
        results[3]=prob[0];
        return results;
    }

    public static String judgeResult(){
        if(result==1)
            return "TRUE";
        else
            return "FALSE";
    }


    

    /**
     * denote the char type bases with number
     * @param char c: a,u,g,c
     * @return int
     */
    private static int nt2Num(char c){
        if(c=='A' || c=='a')
            return 1;
        else if(c=='C' || c=='c')
            return 2;
        else if(c=='G' || c=='g')
            return 3;
        else if(c=='U' || c=='u')
            return 4;
        else return 0;
    }

}
