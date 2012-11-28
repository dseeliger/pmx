from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
import numpy as np


NGMX_XTC = Extension("NGMX_XTC",			
					include_dirs = [np.get_include()+'/numpy', "./C_libs/xdrLIB"],
					sources = ['C_libs/NGMX_XTC.c', 'C_libs/xdrLIB/xdrfile.c'])

NGMX_geometry = Extension("NGMX_geometry",			
					include_dirs = [np.get_include()+'/numpy', "./C_libs/LinearAssignmentSolver"],
					sources = ["C_libs/NGMX_geometry.c", "C_libs/LinearAssignmentSolver/lap.c"])



def configuration(parent_package='', top_path=None):
	config = Configuration(None, parent_package, top_path)
	return config

if __name__ == "__main__":
	setup(name="NGMX", version="0.5",
		author="David Kopfer",
		author_email="dkoepfe@mpibpc.mpg.de",
		license="BSD",
		ext_modules = [NGMX_XTC, NGMX_geometry],
		configuration=configuration)

