package org.panda.causalpath.exceptions;

public class IllegalValueTransformationException extends IllegalArgumentException
{
    public IllegalValueTransformationException() { }

    public IllegalValueTransformationException(String s) {
        super(s);
    }

    public IllegalValueTransformationException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalValueTransformationException(Throwable cause) {
        super(cause);
    }
}
