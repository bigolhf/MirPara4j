/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import javax.swing.JPanel;

/**
 *
 * @author wb
 */
public class MirViewDisPanel extends JPanel implements MouseListener{
    private MirViewData mvd;
    private MirViewPanel mvp;
    private double uw;
    private double uh;
        
    public MirViewDisPanel(MirViewPanel mvp, MirViewData mvd){
        this.mvp=mvp;
        this.mvd=mvd;
            
        this.setLayout(new BorderLayout());
        this.setPreferredSize(new Dimension(600,100));
        this.addMouseListener(this);
    }
        
    public void paintComponent(Graphics g){
        super.paintComponent(g);
        int w=this.getWidth()-10;
        int h=this.getHeight()-10;
            
        BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        drawPlot(img);
        g.drawImage(img, 5, 5, this);
            
    }
        
    public void drawPlot(BufferedImage img){
        Graphics2D gg = (Graphics2D) img.getGraphics();
        int w=img.getWidth();
        int h=img.getHeight();
        gg.setColor(Color.white);
        gg.fillRect(0, 0, w, h);
        int ph=h-30;
        uw=w*1.0/mvd.maxx;
        uh=ph*1.0/mvd.maxbin;
        gg.setColor(Color.black);
        gg.drawLine(0, 0, 0, ph);
        gg.drawLine(0, ph, w, ph);
        gg.setColor(Color.GRAY);
        for(int i=0;i<mvd.dis.length;i++){
            int[] p=mvd.dis[i];
            gg.fill(new Rectangle2D.Double(p[0]*uw, ph-p[2]*uh, (p[1]-p[0])*uw, p[2]*uh));
            gg.draw(new Line2D.Double(p[1]*uw, ph, p[1]*uw, ph+3));
            if(i%10==0){
                gg.draw(new Line2D.Double(p[1]*uw, ph, p[1]*uw, ph+6));
                gg.drawString(p[1]+"", (float)(p[1]*uw), (float)(ph+20));
            }
        }
            
    }

    public void mouseClicked(MouseEvent me) {
        int x=me.getX();   
        int y=me.getY();
        if(x<=5 || x>=this.getWidth()-5) return;
        if(y<=5 || y>=this.getHeight()-35) return;
        x-=5;
        x/=uw;
        mvp.setScrollValues(x, 0);
    }

    public void mousePressed(MouseEvent me) {
            
    }

    public void mouseReleased(MouseEvent me) {
            
    }

    public void mouseEntered(MouseEvent me) {
            
    }

    public void mouseExited(MouseEvent me) {
            
    }
    
}
