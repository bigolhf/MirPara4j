/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.jmirpara;

import java.io.IOException;
import java.io.OutputStream;
import javax.swing.JTextArea;

/**
 *
 * @author wb
 */
public class TextAreaOutputStream extends OutputStream
{
    private static int BUFFER_SIZE = 8192;
    private JTextArea target;
    private byte[] buffer = new byte[BUFFER_SIZE];
    private int pos = 0;

    public TextAreaOutputStream(JTextArea target)
    {
        this.target = target;
    }

    @Override
    public void write(int b) throws IOException
    {
        buffer[pos++] = (byte)b;
        //Append to the TextArea when the buffer is full
        if (pos == BUFFER_SIZE) {
            flush();
        }
    }

    @Override
    public void flush() throws IOException
    {
        byte[] flush = null;
        if (pos != BUFFER_SIZE) {
            flush = new byte[pos];
            System.arraycopy(buffer, 0, flush, 0, pos);
        }
        else {
            flush = buffer;
        }
        
        int count=0;
        for(byte b:flush){
            if((byte)0x08==b){
                count+=1;
            }
        }
        String ns=new String(flush);
        //remove backspace charater '\b'
        ns=ns.substring(0, ns.length()-2*count);        
        target.append(ns);
        target.setCaretPosition(target.getText().length());
        
        
//        target.append(new String(flush));
        pos = 0;
    }
}
