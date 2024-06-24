package org.panda.causalpath.exceptions;

public class BooleanUtil {
    public static boolean parseBoolWithException(String value)
    {
        if ("true".equals(value.toLowerCase())) return true;
        if ("false".equals(value.toLowerCase())) return false;

        throw new BooleanFormatException(value + "cannot be converted to a boolean");
    }
}
