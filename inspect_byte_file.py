import struct
import numpy as np
import sys
import os

def inspect_byte_file(file_path, num_bytes=100):
    """Inspect the binary structure of a .byte file and display details."""
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    print(f"Inspecting file: {file_path}\n")
    
    try:
        with open(file_path, "rb") as f:
            # Read a small portion of the file to get an idea of the structure
            header = f.read(num_bytes)
            
            print(f"First {num_bytes} bytes (raw): {header[:num_bytes]}\n")
            
            # Try interpreting as ASCII text
            try:
                text_preview = header.decode("utf-8")
                print(f"ASCII Interpretation:\n{text_preview}\n")
            except UnicodeDecodeError:
                print("ASCII Interpretation: Not valid text\n")
            
            # Attempt to read as different numerical formats
            f.seek(0)
            raw_data = f.read()
            data_length = len(raw_data)

            # Check for standard floating-point formats
            if data_length % 4 == 0:
                try:
                    float_data = np.frombuffer(raw_data, dtype=np.float32)
                    print(f"Interpreted as 32-bit float (first 10 values): {float_data[:10]}")
                except:
                    print("Could not interpret as 32-bit float.")

            if data_length % 8 == 0:
                try:
                    double_data = np.frombuffer(raw_data, dtype=np.float64)
                    print(f"Interpreted as 64-bit float (first 10 values): {double_data[:10]}")
                except:
                    print("Could not interpret as 64-bit float.")

            # Check for integer formats
            try:
                int_data = np.frombuffer(raw_data, dtype=np.int32)
                print(f"Interpreted as 32-bit integers (first 10 values): {int_data[:10]}")
            except:
                print("Could not interpret as 32-bit integers.")

            try:
                short_int_data = np.frombuffer(raw_data, dtype=np.int16)
                print(f"Interpreted as 16-bit integers (first 10 values): {short_int_data[:10]}")
            except:
                print("Could not interpret as 16-bit integers.")

    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_byte_file.py <path_to_byte_file>")
    else:
        inspect_byte_file(sys.argv[1])


inspect_byte_file(r"C:\Users\Brayden.Holler\OneDrive - afacademy.af.edu\Documents\GitHub\499\Data\llipa.byte", num_bytes=1000)