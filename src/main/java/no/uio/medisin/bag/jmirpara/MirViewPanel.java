/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.ToolTipManager;

/**
 *
 * @author wb
 */
public class MirViewPanel extends JPanel implements AdjustmentListener, MouseListener, MouseMotionListener{
    private MirViewData mvd;
    private JFrame mv;
    private JPanel scalePanel;
    private JPanel seqPanel;
    private JScrollBar hscroll=new JScrollBar();
    private JScrollBar vscroll=new JScrollBar();
    
    private int width;
    private int height;
    private int uw=10;
    private int uh=10;
    private int startx;
    private int endx;
    private int starty;
    private int endy;
    
    int hextent = 0;
    int vextent = 0;
    
    public MirViewPanel(JFrame mv, MirViewData mvd){
        this.mvd=mvd;
        this.mv=mv;
        width=mvd.maxx;
        height=mvd.maxy;
        
        scalePanel=new ScaleCanvas(); 
        scalePanel.setPreferredSize(new Dimension(600,30));
        seqPanel=new MirCanvas();
        seqPanel.setPreferredSize(new Dimension(600,400));

        this.setLayout(new BorderLayout());
        this.add(scalePanel, BorderLayout.NORTH);
        this.add(seqPanel);
        this.add(vscroll, BorderLayout.EAST);
        this.add(hscroll,BorderLayout.SOUTH);
        hscroll.setOrientation(JScrollBar.HORIZONTAL);
        
        setScrollValues(0, 0);
        
        hscroll.addAdjustmentListener(this);
        vscroll.addAdjustmentListener(this);
        
        seqPanel.addMouseListener(this);
        seqPanel.addMouseMotionListener(this);
        ToolTipManager.sharedInstance().setInitialDelay(10);
    }

     public void setScrollValues(int x, int y) {
        
        int vWidth=seqPanel.getWidth();
        int vHeight=seqPanel.getHeight();

        endx=x+vWidth/uw;
        endy=y+vHeight/uh;
        
        hextent = vWidth / uw;
        vextent = vHeight / uh;

        if (hextent > width){
            hextent = width;
        }

        if (vextent > height){
            vextent = height;
        }

        if ((hextent + x) > width){
            x = width - hextent;
        }

        if ((vextent + y) > height){
            y = height - vextent;
        }

        if (y < 0){
            y = 0;
        }

        if (x < 0){
            x = 0;
        }

        hscroll.setValues(x, hextent, 0, width);
        vscroll.setValues(y, vextent, 0, height);
    }
    
    @Override
    public void adjustmentValueChanged(AdjustmentEvent evt) {
        int oldX = startx;
        int oldY = starty;
        
        int vWidth=seqPanel.getWidth();
        int vHeight=seqPanel.getHeight();

        if (evt.getSource() == hscroll){
            int x = hscroll.getValue();
            startx=x;
            endx=x+vWidth/uw;
        }

        if (evt.getSource() == vscroll){
            int offy = vscroll.getValue();
            starty=offy;
            endy=offy+vHeight/uh;     
        }


        int scrollX = startx - oldX;
        int scrollY = starty - oldY;

    
        // Make sure we're not trying to draw a panel
        // larger than the visible window
        if (scrollX > endx - startx){
            scrollX = endx - startx;
        }
        else if (scrollX < startx - endx){
            scrollX = startx - endx;
        }

        if (scrollX != 0 || scrollY != 0){
            repaint();
        }
    }
    
    public int findMir(int x, int y){
        int mp=-1;
        int mx=startx+x/uw+1;
        int my=starty+y/uh;
        int[][] dis=mvd.dis;
        int pos=0;
        if(width>1000){
        for(int i=0;i<dis.length;i++){
            if(dis[i][0]<=mx && mx<=dis[i][1]){
                if(i==0) pos=my-1;
                else{
                    pos-=dis[i-1][2];                    
                }
                if(pos<0) pos=0;
                break;
            }
            pos+=dis[i][2];
        }
        }
        ArrayList<int[]> mirs=mvd.mirs;
        for(int i=pos;i<mirs.size();i++){
            int[] mir=mirs.get(i);
            if(mir[0]<=mx && mx<=mir[1] && mir[2]==my){
                mp=i;
                break;
            }
        }
        return mp;
    }
    
    public void paintComponent(Graphics g){ 
        setScrollValues(startx, starty);
    }

    public void mouseClicked(MouseEvent me) {
//        if(me.getClickCount()!=2) return;
        int selected=findMir(me.getX(),me.getY());
        if(selected<0) return;
        String info="";
        try {
            info=mvd.getDetailInfo(selected);
        } catch (IOException ex) {
            System.out.println("cannot read file!");
        }
        ((ResultFrame)mv).infoBox.setText(info);
        ((ResultFrame)mv).infoBox.setCaretPosition(1);

    }

    public void mousePressed(MouseEvent me) {
        
    }

    public void mouseReleased(MouseEvent me) {
        
    }

    public void mouseEntered(MouseEvent me) {
        
    }

    public void mouseExited(MouseEvent me) {
        
    }

    public void mouseDragged(MouseEvent me) {
        
    }

    public void mouseMoved(MouseEvent me) {
        int selected=findMir(me.getX(),me.getY());
        ((MirCanvas)seqPanel).selected=selected;
        seqPanel.repaint();

        if(selected<0) {
            seqPanel.setToolTipText(null);           
            return;
        }
        
        int[] mi=mvd.mirs.get(selected);
        String info="<html><b>"+mvd.name + "_MIR_" + mi[0]
                + "-" + (mi[1]-mi[0]+1);
        int h=-1;
        ArrayList<int[]> hits=mvd.hits;
        for(int i=0;i<hits.size();i++){
            int[] hit=hits.get(i);
            if(mi==hit){
                h=i;
                break;
            }
        }
        if(h<0){
            info+="</b></html>";
        }
        else{
            String hi=mvd.hitNames.get(h);
            info+="</b><br>hits:<br><i>"+hi+"</i></html>";
        }        
        seqPanel.setToolTipText(info);
    }
    
    class ScaleCanvas extends JPanel{
        public void paintComponent(Graphics g){
            drawScale(g, startx, endx, getWidth(), getHeight());
        }

        // scalewidth will normally be screenwidth,
        public void drawScale(Graphics g, int startx, int endx, int width, int height){
            Graphics2D gg = (Graphics2D) g;
            Font font=new Font("Courier New", Font.PLAIN, 12);
            gg.setFont(font);


            // Fill in the background
            gg.setColor(new Color(238,238,238));
            gg.fillRect(0, 0, width, height);
            gg.setColor(Color.black);        

            FontMetrics fm = gg.getFontMetrics(font);
            int y = 20 - fm.getDescent();

            String label;
            for(int i=startx;i<=endx;i++){
                int v=i+1;
            
                if((v%10)==0){
                    label=String.valueOf(v);
                
                    gg.drawString(label, (i-startx+0.5f-label.length()*1f/2)*uw,y);
                    gg.drawLine((int)((i-startx+0.5)*uw), y, (int)((i-startx+0.5)*uw), y+fm.getDescent()*2);
                }
                else{
                    gg.drawLine((int)((i-startx+0.5)*uw), y+fm.getDescent(), (int)((i-startx+0.5)*uw), y+fm.getDescent()*2);
                }
            }
            gg.drawLine(0, y+fm.getDescent()*2, (endx-startx+1)*uw, y+fm.getDescent()*2);       

        }
    }
    
    class MirCanvas extends JPanel{
        public BufferedImage img;
        public int selected=-1;
        
        @Override
        public void paintComponent(Graphics g){
            super.paintComponent(g);
            draw(startx,endx,starty,endy);
            g.drawImage(img, 0, 0, this);
        }
        
        public void draw(int bx1, int bx2, int by1, int by2){
            int imgWidth=getWidth();
            int imgHeight=getHeight();
        
            img = new BufferedImage(imgWidth, imgHeight,
                    BufferedImage.TYPE_INT_RGB);
            Graphics2D gg = (Graphics2D) img.getGraphics();

            gg.setColor(Color.white);
            gg.fillRect(0, 0, imgWidth, imgHeight);
            if(by1<=0){
                gg.setColor(Color.LIGHT_GRAY);
                int end=width*uw>imgWidth?imgWidth:width*uw;
                gg.fill(new Rectangle2D.Double(0,0, end,uh));
            
                gg.setColor(Color.black);
                for(int[] p:mvd.pris){
                    if(p[1]<=bx1) continue;
                    if(p[0]>=bx2) break;
                    gg.fill(new Rectangle2D.Double((p[0]-bx1-1)*uw,0,(p[1]-p[0]+1)*uw,uh));
                }
            }
            gg.setColor(Color.blue);
            for(int[] p:mvd.mirs){
                if(p[1]<=bx1) continue;
                if(p[0]>=bx2) break;
                if(p[2]<by1) continue;
                gg.fill(new Rectangle2D.Double((p[0]-bx1-1)*uw,(p[2]-by1)*uh,(p[1]-p[0]+1)*uw,uh*0.8));
            }
            gg.setColor(Color.red);
            for(int[] p:mvd.hits){
                if(p[1]<=bx1) continue;
                if(p[0]>=bx2) break;
                if(p[2]<by1) continue;
                gg.fill(new Rectangle2D.Double((p[0]-bx1-1)*uw,(p[2]-by1)*uh,(p[1]-p[0]+1)*uw,uh*0.8));
            }
            if(selected>=0){
                gg.setColor(Color.yellow);
                int[] p=mvd.mirs.get(selected);
                gg.fill(new Rectangle2D.Double((p[0]-bx1-1)*uw,(p[2]-by1)*uh,(p[1]-p[0]+1)*uw,uh*0.8));
            }
            
        }
    }
    
}
