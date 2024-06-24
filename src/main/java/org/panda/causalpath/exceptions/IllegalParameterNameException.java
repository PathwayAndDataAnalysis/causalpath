package org.panda.causalpath.exceptions;

public class IllegalParameterNameException extends  IllegalArgumentException {
    public IllegalParameterNameException() { }

    public IllegalParameterNameException(String s) {
        super(s);
    }

    public IllegalParameterNameException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalParameterNameException(Throwable cause) {
        super(cause);
    }
}
