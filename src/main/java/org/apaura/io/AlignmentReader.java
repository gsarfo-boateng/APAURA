package org.apaura.io;

import java.io.Closeable;
import java.io.IOException;

public interface AlignmentReader extends Closeable {
    
    AlignmentRecord next() throws IOException;
}