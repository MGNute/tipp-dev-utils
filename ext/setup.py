'''
To run this setup file, enter the subdirectory that this package sits in
(i.e. tipp_dev_utils/ext) and run the following:

    python setup.py build_ext --inplace

That will compile the module into a file called "nwops.<something>" where
the <something> depends on the OS and probably the compiler. On windows,
with Visual Studio, it is "nwops.cp36-win_amd64.pyd". From the main
tipp_dev_utils folder, you can then use:

    import ext.nwops

to import the module and use the functions. Try help(ext.nwops) for the
docstrings.
'''

# import platform, os
# if platform.system()=='Windows':
#     os.environ['VS90COMNTOOLS']=os.environ['VS140COMNTOOLS']

# from distutils.core import setup, Extension
from setuptools import setup, Extension, Command
import numpy as np
# import numpy.distutils.misc_util

c_ext = Extension("nwops", ["_nwops.c", "nwops.c"])

setup(
    ext_modules=[c_ext],
    # include_dirs=["."],
    include_dirs=np.get_include()
)
