package nicloj.header;

/** Raised for malformed, inconsistent or unsupported NIfTI data. */
public class NiftiError extends RuntimeException {
    private static final long serialVersionUID = 1L;

    public NiftiError(String message) { super(message); }

    public NiftiError(String message, Throwable cause) { super(message, cause); }
}
