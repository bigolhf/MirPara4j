/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Date;
import java.util.Timer;
import java.util.TimerTask;

/**
 *
 * @author wb
 */
public class MouseClickListener extends MouseAdapter{
    private static boolean flag=false;//用来判断是否已经执行双击事件
    private static int clickNum=0;//用来判断是否该执行双击事件

    @Override
    public void mouseClicked(MouseEvent e) {
        final MouseEvent me=e;//事件源

        MouseClickListener.flag=false;//每次点击鼠标初始化双击事件执行标志为false

        if (MouseClickListener.clickNum == 1) {//当clickNum==1时执行双击事件
            this.mouseDoubleClicked(me);//执行双击事件
            MouseClickListener.clickNum=0;//初始化双击事件执行标志为0
            MouseClickListener.flag=true;//双击事件已执行,事件标志为true
            return;
        }

        //定义定时器
        Timer timer=new Timer();

        //定时器开始执行,延时0.2秒后确定是否执行单击事件
        timer.schedule(new TimerTask() {
            private int n=0;//记录定时器执行次数
            public void run() {
                if(MouseClickListener.flag){//如果双击事件已经执行,那么直接取消单击执行
                    n=0;
                    MouseClickListener.clickNum=0;
                    this.cancel();
                    return;
                }
                if (n == 1) {//定时器等待0.2秒后,双击事件仍未发生,执行单击事件
                    mouseSingleClicked(me);//执行单击事件
                    MouseClickListener.flag = true;
                    MouseClickListener.clickNum=0;
                    n=0;
                    this.cancel();
                    return;
                }
                clickNum++;
                n++;
            }
        },new Date(),200);
    }

    /**
    * 鼠标单击事件
    * @param e 事件源参数
    */
    public void mouseSingleClicked(MouseEvent e){
        System.out.println("Single Clicked!");
    } 

    /**
    * 鼠标双击事件
    * @param e 事件源参数
    */
    public void mouseDoubleClicked(MouseEvent e){
        System.out.println("Doublc Clicked!");
    }
}
