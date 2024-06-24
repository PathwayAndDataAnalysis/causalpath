package org.panda.causalpath.log;

import java.io.File;

public class FileUtility {
    public static boolean fileExists(String filepath)
    {
        return new File(filepath).isFile();
    }

    public static boolean directoryExists(String filepath)
    {
        return new File(filepath).isDirectory();
    }
}
