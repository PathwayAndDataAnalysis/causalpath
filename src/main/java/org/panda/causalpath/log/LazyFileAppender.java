package org.panda.causalpath.log;

import ch.qos.logback.core.FileAppender;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class LazyFileAppender<E> extends FileAppender<E> {
    private boolean fileCreated = false;

    @Override
    public void start() {
        // Override the start method but do not create the file yet
        if (encoder == null) {
            addError("No encoder set for the appender named \"" + name + "\".");
            return;
        }
        started = true;
    }

    @Override
    public void append(E eventObject) {
        if (!fileCreated) {
            createFile();
        }
        super.append(eventObject);
    }

    private synchronized void createFile() {
        if (!fileCreated) {
            try {
                File file = new File(getFile());
                if (!file.getParentFile().exists()) {
                    file.getParentFile().mkdirs();
                }
                OutputStream os = new FileOutputStream(file, isAppend());
                setOutputStream(os);
                fileCreated = true;
                addInfo("File created: " + getFile());
            } catch (IOException e) {
                addError("Failed to create file [" + getFile() + "].", e);
            }
        }
    }
}
