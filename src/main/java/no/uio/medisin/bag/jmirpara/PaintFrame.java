/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.jmirpara;

import java.awt.Color;

import java.awt.Frame;

import java.awt.Graphics;

public class PaintFrame extends Frame {

public PaintFrame() {

setBounds(200, 100, 400, 500);

setVisible(true);

}

public void paint(Graphics g) {

g.setColor(Color.BLUE);

g.fillOval(30, 30, 50, 50);

g.setColor(Color.GREEN);

g.fillRect(100, 100, 50, 50);

g.setColor(Color.RED);

g.drawLine(300, 200, 400, 500);

}

public static void main(String[] args) {

new PaintFrame();

}

}