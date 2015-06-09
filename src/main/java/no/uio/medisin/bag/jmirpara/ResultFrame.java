/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import javax.swing.BorderFactory;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EtchedBorder;

/**
 *
 * @author wb
 */
public class ResultFrame extends JFrame{
    private String title="Results";
    private File mirFile;
    private String[] seqs;
    private MirViewPanel mvp;
    private MirViewDisPanel mdp;
    private MirViewData mvd;
    
    private JLabel infoPanel;
    public JTextArea infoBox;
    private JPanel workPanel;
    private JPanel oneSeqPanel;
    private JComboBox seqsBox;

    
    public ResultFrame(){
//        MirViewPanel mvp=new MirViewPanel(mirFile);
//        this.add(mvp);
        this.setTitle(title);

        infoPanel=new JLabel();
        infoPanel.setPreferredSize(new Dimension(800,300));
        infoPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
        infoPanel.setLayout(new BorderLayout());
        infoBox=new JTextArea();
        infoBox.setEditable(false);
        infoBox.setFont(new Font("Courier New", Font.PLAIN, 12));
        infoPanel.add(new JScrollPane(infoBox));
        
            
        oneSeqPanel=new JPanel();
        oneSeqPanel.setLayout(new BorderLayout());
        oneSeqPanel.setBackground(Color.white);
        oneSeqPanel.setPreferredSize(new Dimension(800,300));
        workPanel=new JPanel();
        workPanel.setLayout(new BorderLayout());
        workPanel.add(oneSeqPanel);
        this.add(workPanel);
        this.add(infoPanel,BorderLayout.SOUTH);

        this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        this.pack();
//        this.setSize(800,350);
        this.setLocationRelativeTo(null);//in the center of screen
//        displayAtCenter(this);
        
        this.setVisible(true);
    }
    
    public void showResult(){
        try {
            mvd=new MirViewData(mirFile);
        } catch (FileNotFoundException ex) {
            System.err.println("file not found!");
        } catch (IOException ex) {
            System.err.println("file io error");
        }
        mvp=new MirViewPanel(this,mvd);
        mdp=new MirViewDisPanel(mvp,mvd);
        oneSeqPanel.removeAll();
        oneSeqPanel.add(mvp);
        oneSeqPanel.add(mdp,BorderLayout.SOUTH);
        oneSeqPanel.validate();
    }
    
    public void setResultList(){
        String[] seqsInfo=new String[seqs.length];
        for(int i=0;i<seqs.length;i++){
            String info="";
            try {
                info=MirViewData.summary(seqs[i]);
            } catch (FileNotFoundException ex) {
                System.out.println("file not found!");
            } catch (IOException ex) {
                System.out.println("cannot read the file");
            }
            
            String[] item=info.split("\t");
            seqsInfo[i]="Sequence: "+item[0]+"  length:"+item[1]+"  number of miRNAs: "+item[2];
        }
        
        seqsBox=new JComboBox(seqsInfo);
        seqsBox.addItemListener(new ItemListener(){

            public void itemStateChanged(ItemEvent ie) {
                mirFile=new File(seqs[seqsBox.getSelectedIndex()]);
                showResult();
            }
        });
        workPanel.add(seqsBox,BorderLayout.NORTH);
        workPanel.validate();
        mirFile=new File(seqs[0]);
    }

    /**
     * @return the mirFile
     */
    public File getMirFile() {
        return mirFile;
    }

    /**
     * @param mirFile the mirFile to set
     */
    public void setMirFile(File mirFile) {
        this.mirFile = mirFile;
    }

    /**
     * @return the seqs
     */
    public String[] getSeqs() {
        return seqs;
    }

    /**
     * @param seqs the seqs to set
     */
    public void setSeqs(String[] seqs) {
        this.seqs = seqs;
    }
    
}
