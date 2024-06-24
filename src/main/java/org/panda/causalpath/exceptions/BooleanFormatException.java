package org.panda.causalpath.exceptions;

public class BooleanFormatException extends IllegalArgumentException
{
    public BooleanFormatException()
    {
        super();
    }

    public BooleanFormatException(String errorMessage)
    {
        super(errorMessage);
    }

    public BooleanFormatException(String errorMessage, Throwable err)
    {
        super(errorMessage, err);
    }

    public BooleanFormatException(Throwable cause) {
        super(cause);
    }
}
