package org.panda.causalpath.exceptions;

public class IllegalDataTypeException extends IllegalArgumentException {
    private String validValues;
    public IllegalDataTypeException() { }

    public IllegalDataTypeException(String s, String validValues) {
        super(s);
        this.validValues = validValues;
    }

    public IllegalDataTypeException(String message, Throwable cause, String validValues) {
        super(message, cause);
        this.validValues = validValues;
    }

    public IllegalDataTypeException(Throwable cause) {
        super(cause);
    }

    public String getValidValues()
    {
        return validValues;
    }
}
