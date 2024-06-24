package org.panda.causalpath.exceptions;

public class IllegalNetworkResourceException extends IllegalArgumentException {
    public IllegalNetworkResourceException() { }

    public IllegalNetworkResourceException(String s) {
        super(s);
    }

    public IllegalNetworkResourceException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalNetworkResourceException(Throwable cause) {
        super(cause);
    }
}
