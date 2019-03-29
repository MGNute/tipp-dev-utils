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
    include_dirs=np.get_include(),
)
