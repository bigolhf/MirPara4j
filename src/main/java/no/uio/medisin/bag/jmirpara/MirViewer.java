/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.UIManager;
import javax.swing.border.EtchedBorder;

/**
 *
 * @author wb
 */
public class MirViewer extends ResultFrame implements ActionListener{
    private String title="MirViewer";
    private MirViewPanel mvp;
    private MirViewDisPanel mdp;
    private MirViewData mvd;
    
    private JFileChooser fileChooser;
    
    private JLabel infoPanel;
    public JTextArea infoBox;
    private JPanel workPanel;
    
    private JMenu file=new JMenu("File");
    private JMenu help=new JMenu("Help");
    private JMenu analysis=new JMenu("Analysis");
    private JMenuItem openItem;
    private JMenuItem saveItem;
    private JMenuItem exitItem;
    private JMenuItem aboutItem;
    
    
    public MirViewer(){
        super();
        
        this.setTitle(title);  
        this.setJMenuBar(makeMenuBar());
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.pack();
        this.setVisible(true);
    }
    
    private JMenuBar makeMenuBar(){
        JMenuBar menubar=new JMenuBar();
        openItem=new JMenuItem("Open");
        openItem.addActionListener(this);
        
        saveItem=new JMenuItem("Save");
        saveItem.setIcon(new ImageIcon(""));
        saveItem.addActionListener(this);
        
        exitItem=new JMenuItem("Exit");
        exitItem.setIcon(new ImageIcon(""));
        exitItem.addActionListener(this);
        
        
        aboutItem=new JMenuItem("About");
        aboutItem.setIcon(new ImageIcon(""));
        aboutItem.addActionListener(this);
        
        file.add(openItem);
        file.add(saveItem);
        file.add(exitItem);
        help.add(aboutItem);
        
        menubar.add(file);
        menubar.add(analysis);
        menubar.add(help);
        
        return menubar;
    }

    public void actionPerformed(ActionEvent ae) {
        if(ae.getSource()==openItem){
            
            fileChooser=new JFileChooser();           
            
            fileChooser.setCurrentDirectory(new File("."));
            
            int result=fileChooser.showOpenDialog(this);
            if(result==JFileChooser.APPROVE_OPTION){
                File fileObj=fileChooser.getSelectedFile();                              
                this.setTitle("MirViewer--"+fileObj.getPath());
                this.setSeqs(new String[]{MirViewData.getBaseFileName(fileObj.getAbsolutePath())});
                this.setResultList();
                this.showResult();              
            }
        }
    }
    
    public static void main(String[] args){
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        MirViewer mv=new MirViewer();      
    }  
    
}
