from setuptools import setup

setup(name='depmapomics',
      package_dir = {'': 'src'},
      install_requires=open("requirements.txt").readlines())