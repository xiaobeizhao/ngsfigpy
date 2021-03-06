
import sys
import os
import glob


try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup

if sys.version_info < (2, 7) or sys.version_info >= (3,0):
  raise NotImplementedError("Sorry, you need Python 2.7 use `xmisc`.")


def read(fname):
  return open(os.path.join(os.path.dirname(__file__), fname)).read()

  
def get_data_files(x,Fext=[".txt",".png",".pdf"]):
  data_files_ = []
  dirs_ = glob.glob('%s/'%x)
  while len(dirs_):
    for e in dirs_:
      files_ = glob.glob(e+'*')
      files_ = [ee for ee in files_ if os.path.isfile(ee) and os.path.splitext(ee)[1] in Fext ]
      data_files_.append((e, files_))
      dirs_=glob.glob('%s/*/'%e)
  return(data_files_)


config = {
  'name': 'ngsfig',
  'version': '0.1.1',
  'author': 'Xiaobei Zhao',
  'author_email': 'xiaobei_zhao(at)med.unc.edu',
  'description': ('Classes and functions for NGS data illustration'),
  'long_description': read('README.rst'),
  'license': 'GNU LGPL',
  'packages': ['ngsfig','ngsfig.utils','ngsfig.graphics'],
  'scripts': [],
  'data_files': get_data_files("data"),
  'url': 'https://github.com/xiaobeizhao/ngsfigpy',
  'download_url': '',
  'install_requires': ['nose',
                       'setuptools>=33.1.1',
                       'numpy>=1.11.3',
                       'pandas>=0.19.2',
                       'tabulate>=0.7.5',
                       'matplotlib>=2.0.0',
                       'xmisc'
                       ]
}

setup(include_package_data=True, **config)
