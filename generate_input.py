
#Set up Gas Box
import ctypes
import pathlib

libname = pathlib.Path().absolute() / "functions.so"
c_lib = ctypes.CDLL(libname)

c_lib.