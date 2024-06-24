package org.panda.causalpath.log;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
public class CPLogger
{
    public static boolean isInitialized = false;
    public static final String logbackConfigFilePath = "src/main/resources/logback.xml";
    private static String outputFilePath;
    public static Logger paramError;
    public static Logger paramWarning;
    public static Logger dataError;
    public static Logger dataWarning;
    public static void Initialize(String output)
    {
        isInitialized = true;
        outputFilePath = output;
        System.setProperty("InvalidParameterFilePath", outputFilePath);
        paramError = LoggerFactory.getLogger("InvalidParameterErrors");
        paramWarning = LoggerFactory.getLogger("InvalidParameterWarnings");
        dataError = LoggerFactory.getLogger("InvalidDataErrors");
        dataWarning = LoggerFactory.getLogger("InvalidDataWarnings");
    }

    public static String getParameterErrorLogFilePath()
    {
        return outputFilePath + "/invalid-parameter-errors.log";
    }
    public static boolean parameterErrorLogFileExists()
    {
         return FileUtility.fileExists(getParameterErrorLogFilePath());
    }

    public static String getParameterWarningLogFilePath()
    {
        return outputFilePath + "/invalid-parameter-warnings.log";
    }
    public static boolean parameterWarningLogFileExists()
    {
        return FileUtility.fileExists(getParameterWarningLogFilePath());
    }
}
