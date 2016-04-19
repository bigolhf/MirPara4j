/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.EtchedBorder;

/**
 *
 * @author wb
 */
public class JmiGUI extends JFrame implements ActionListener{
    private double version=4.0;
    
    private String title;
    private JPanel formPanel;
    private JPanel logPanel;
    private JTextField infileF;
    private JTextField outdirF;
    private JComboBox modelCB;
    private JComboBox levelCB;
    private JTextField lengthF;
    private JTextField cutoffF;
    private JTextField stepF;
    private JTextField windowF;
    private JButton defaultButton;
    private JButton scanButton;
    private JTextArea logArea;
    private JTextField processF;
    private JButton clearButton;
    private JButton helpButton;
    
    private File infile;
    private File outdir;
    private String model;
    private int level;
    private int length;
    private double cutoff;
    private int window;
    private int step;
    private MiRNAProcessingPipeLine pl;
    
    public JmiGUI(){
        
        title="miRPara4j";
        
        makeFormPanel();
        makeLogPanel();
        this.setLayout(new BorderLayout());
        this.add(formPanel,BorderLayout.WEST);
        this.add(logPanel);
        this.setTitle(title);
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.pack();
        this.setResizable(false);
//        this.setSize(800,350);
        this.setLocationRelativeTo(null);//in the center of screen
//        displayAtCenter(this);
        
        this.setVisible(true);
    }
    
    private void makeFormPanel(){
        
        formPanel=new JPanel();
        formPanel.setLayout(new BorderLayout(1,4));
        
        JPanel paraPanel=new JPanel();
        paraPanel.setLayout(new GridLayout(8,2,10,10));
        paraPanel.setPreferredSize(new Dimension(300,300));
        paraPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));

        JLabel infileLabel=new JLabel("Input file");
        infileLabel.setToolTipText("Choose input sequence file (FASTA format)");
        infileF=new JTextField();
        paraPanel.add(infileLabel);
        paraPanel.add(infileF);

        JLabel outdirLabel=new JLabel("Output directory ");
        outdirLabel.setToolTipText("Choose output directory to store result files");
        outdirF=new JTextField();
        paraPanel.add(outdirLabel);
        paraPanel.add(outdirF);

        JLabel modelLabel=new JLabel("Model");
        modelLabel.setToolTipText("Choose prediction model");
        String[] models=new String[]{"overall","animal","plant","virus"};
        modelCB=new JComboBox(models);
        paraPanel.add(modelLabel);
        paraPanel.add(modelCB);

        JLabel levelLabel=new JLabel("Level");
        levelLabel.setToolTipText("<html>Levels with lower numbers are less 'strict' -<br> i.e. you will get more positive predictions, but more false positives. <br>Using a higher level will filter out more of the false positives, <br>but you will miss some of the true positives. <br>Level 15 is a good compromise, but you can try different Levels <br>to see how it affects your predictions.</html>");
        String[] levels=new String[20];
        for(int i=0;i<20;i++) levels[i]=(i+1)+"";
        levelCB=new JComboBox(levels);
        paraPanel.add(levelLabel);
        paraPanel.add(levelCB);

        JLabel lengthLabel=new JLabel("Length");
        lengthLabel.setToolTipText("Specify the minimum length of pri-miRNA sequences");
        lengthF=new JTextField("60");
        paraPanel.add(lengthLabel);
        paraPanel.add(lengthF);
        
        JLabel cutoffLabel=new JLabel("Cutoff");
        cutoffLabel.setToolTipText("<html>Cutoffs with smaller numbers are less 'strict'. <br>The SVM probability range from -1 to +1, <br>indicated the posibilities to be TRUE negative or TRUE positive.<html>");
        cutoffF=new JTextField("0.8");
        paraPanel.add(cutoffLabel);
        paraPanel.add(cutoffF);

        JLabel windowLabel=new JLabel("Window");
        windowLabel.setToolTipText("The window size for the long sequence sliding scan");
        windowF=new JTextField("500");
        paraPanel.add(windowLabel);
        paraPanel.add(windowF);

        JLabel stepLabel=new JLabel("Step");
        stepLabel.setToolTipText("The step size for the long sequence sliding scan");
        stepF=new JTextField("250");
        paraPanel.add(stepLabel);
        paraPanel.add(stepF);  
        
        Box exeBox=Box.createHorizontalBox();
        defaultButton=new JButton("Default");
        defaultButton.addActionListener(this);
        scanButton=new JButton("Scan");
        scanButton.addActionListener(this);
        exeBox.add(Box.createHorizontalGlue());
        exeBox.add(defaultButton);
        exeBox.add(Box.createHorizontalStrut(5));
        exeBox.add(scanButton);
        exeBox.add(Box.createHorizontalStrut(11));
        
        formPanel.add(paraPanel);
        formPanel.add(exeBox,BorderLayout.SOUTH);
        formPanel.add(new JPanel(),BorderLayout.EAST);
        formPanel.add(new JPanel(),BorderLayout.NORTH);
        formPanel.add(new JPanel(),BorderLayout.WEST);
        
        MyMouseClickListener mmcl=new MyMouseClickListener();
        infileF.addMouseListener(mmcl);
        outdirF.addMouseListener(mmcl);
        lengthF.addMouseListener(mmcl);
        cutoffF.addMouseListener(mmcl);
        windowF.addMouseListener(mmcl);
        stepF.addMouseListener(mmcl);
    }
    
    private void makeLogPanel(){
        
        logPanel=new JPanel();
        logPanel.setLayout(new BorderLayout(1,4));
        logPanel.setPreferredSize(new Dimension(575,300));
        
        logArea=new JTextArea();
        logArea.setColumns(80);
        logArea.setLineWrap(true);
        logArea.setFont(new Font("Courier New", Font.PLAIN, 12));       
        logArea.setEditable(false);
        
        Box logStatusBox=Box.createHorizontalBox();
        processF=new JTextField();
        processF.setPreferredSize(new Dimension(200,24));
        logStatusBox.add(processF);
        clearButton=new JButton("clear");
        clearButton.addActionListener(this);
        logStatusBox.add(clearButton);
        helpButton=new JButton("help");
        helpButton.addActionListener(this);
        logStatusBox.add(helpButton);
        logPanel.add(new JScrollPane(logArea));
        logPanel.add(logStatusBox,BorderLayout.SOUTH);
        logPanel.add(new JPanel(),BorderLayout.NORTH);
        
        System.setOut(new PrintStream(new TextAreaOutputStream(logArea),true));
        System.setErr(new PrintStream(new TextAreaOutputStream(logArea),true));
        
        MirParaCmd.printBanner();
        MirParaCmd.print_help();

    }
    
    public static void main(String[] args){
//        try{
//        String plaf = UIManager.getSystemLookAndFeelClassName();
//        String gtkLookAndFeel = "com.sun.java.swing.plaf.gtk.GTKLookAndFeel";
//
//        if (plaf.contains("metal")){
//            UIManager.setLookAndFeel(gtkLookAndFeel);
//        }
//      
//        UIManager.setLookAndFeel(plaf);
//        }
//        catch(Exception e){
//            System.out.println("Could not set look and feel, exception : " + e.toString());
//        }
        
        JmiGUI mg=new JmiGUI();
        
    }

    public void actionPerformed(ActionEvent ae) {
       if(ae.getSource()==defaultButton){
            modelCB.setSelectedIndex(0);
            levelCB.setSelectedIndex(0);
            lengthF.setText("60");
            cutoffF.setText("0.8");
            windowF.setText("500");
            stepF.setText("250");           
        }
        else if(ae.getSource()==scanButton){
            pl=new MiRNAProcessingPipeLine();
            
            if(infile==null){
                infileF.setBackground(Color.red);
                System.err.println("lack input file!");
                return;
            }
            pl.setInputFilename(infile.toString());
            this.outdir=new File(outdirF.getText());
            pl.setWorkingDir(outdir);
            this.model=(String)(modelCB.getSelectedItem());
            pl.setModel(model);
            this.level=Integer.parseInt((String)(levelCB.getSelectedItem()));
            pl.setLevel(level);
            
            try{
                this.length=Integer.parseInt(lengthF.getText());
            }catch(Exception e){
                lengthF.setBackground(Color.red);
                System.err.println("Error:length should not be less than 0");
                return;
            }
            if(length<=0){
                lengthF.setBackground(Color.red);
                System.err.println("Error:length should not be less than 0");
                return;
            }
            pl.setDistance(length);
            
            try{
                this.cutoff=Double.parseDouble(cutoffF.getText());
            }catch(Exception e){
                cutoffF.setBackground(Color.red);
                System.err.println("Error:cutoff should be a value between 0 to 1");
                return;
            }
            if(cutoff<0 || cutoff>1){
                cutoffF.setBackground(Color.red);
                System.err.println("Error:cutoff should be a value between 0 to 1");
                return;
            }
            pl.setCutoff(cutoff);
            
            try{
                this.window=Integer.parseInt(windowF.getText());
            }catch(Exception e){
                windowF.setBackground(Color.red);
                System.err.println("Error:window size should not be less than 0");
                return;
            }
            if(window<=0){
                windowF.setBackground(Color.red);
                System.err.println("Error:window size should not be less than 0");
                return;
            }
            pl.setWindow(window);
            
            try{
                this.step=Integer.parseInt(stepF.getText());
            }catch(Exception e){
                stepF.setBackground(Color.red);
                System.err.println("Error:step size should not be less than 0");
                return;
            }
            if(step<=0){
                stepF.setBackground(Color.red);
                System.err.println("Error:step size should not be less than 0");
                return;
            }
            pl.setDistance(step);
            
            this.logArea.append("\nCommand: java -cp JmiPara2.jar jmipara2.JmiPara2CMD -m "
                    +model+" -l "+level+" -d "+length+" -c "+cutoff+" -w "+window
                    +" -s "+step+" "+infileF.getText()+" "+outdirF.getText()+"\n");
            this.logArea.append("\nrunning...\n");            

                    
            final Thread scanT=new Thread(){
                @Override
                public void run(){                  
                    try {
                        pl.predictMiRNAsInQuerySequences();
                    } catch (IOException ex) {
                        Logger.getLogger(JmiGUI.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            };
            scanT.start();
                
            final Thread processT=new Thread(){
                @Override
                public void run(){
                    double progress=0;
                    processF.setForeground(Color.blue);
                    while(progress<100){
                        JmiGUI.this.processF.setText("scanned "+progress+"%");
                    }
                    progress=100;
                    JmiGUI.this.processF.setText("scanned "+progress+"%");
                    while(scanT.isAlive()){
                        while(progress<100){
                            JmiGUI.this.processF.setText("scanned "+progress+"%");
                        }
                        progress=100;
                        JmiGUI.this.processF.setText("scanned "+progress+"%");
                    }
                    progress=100;
                    JmiGUI.this.processF.setText("scanned "+progress+"%");
                }
            };
            processT.start();
            
            Thread resultT=new Thread(){
                @Override
                public void run(){
                    try {
                        scanT.join();
                        processT.join();
                    } catch (InterruptedException ex) {
                        Logger.getLogger(JmiGUI.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    ArrayList<String> rs=pl.getResults();
                    ResultFrame rf=new ResultFrame();
                    rf.setSeqs(rs.toArray(new String[0]));
                    rf.setResultList();
                    rf.showResult();
                }
            };
            resultT.start();
            
        }
        else if(ae.getSource()==helpButton){
            MirParaCmd.printBanner();
            MirParaCmd.print_help();
        }
        else if(ae.getSource()==clearButton){
            logArea.setText("");
        }
       
    }
    
    class MyMouseClickListener extends MouseClickListener{
        
        @Override
        public void mouseSingleClicked(MouseEvent e){
            if(e.getButton()!=MouseEvent.BUTTON1) return;
            
            Object s=e.getSource();
            
            ((JTextField)s).setBackground(Color.white);
                   
            if(s==infileF){
                JFileChooser fileChooser=new JFileChooser();                      
                fileChooser.setCurrentDirectory(new File("."));           
                int result=fileChooser.showOpenDialog(JmiGUI.this);
                if(result==JFileChooser.APPROVE_OPTION){
                    File fileObj=fileChooser.getSelectedFile();                              
                    JmiGUI.this.setTitle("GRecScan--"+fileObj.getPath());
                    JmiGUI.this.infile=fileObj; 
                    JmiGUI.this.infileF.setText(fileObj.getAbsolutePath());
                    JmiGUI.this.outdirF.setText(fileObj.getParent());
                    
                }
            }
            else if(s==outdirF){
                JFileChooser fileChooser=new JFileChooser();                      
                fileChooser.setCurrentDirectory(new File("."));
                fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                int result=fileChooser.showOpenDialog(JmiGUI.this);
                if(result==JFileChooser.APPROVE_OPTION){
                    File fileObj=fileChooser.getSelectedFile();                              
                    JmiGUI.this.outdir=fileObj; 
                    JmiGUI.this.outdirF.setText(fileObj.getAbsolutePath());
                }
            }
            
        
        }
        
        @Override
        public void mouseDoubleClicked(MouseEvent e){
            Object s=e.getSource();
            if(s==infileF){
                infileF.requestFocus();
            }
            else if(s==outdirF){
                outdirF.requestFocus();
            }
        }
    }

    
}
