from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(name='depmapomics',
      description='Pipelines for the analysis of the DepMap omics data',
      url='https://github.com/broadinstitute/depmapomics',
      install_requires=open("requirements.txt").readlines(),
      packages=['depmapomics'],
      long_description=long_description,
      author='Jeremie Kalfon',
      author_email='jkobject@gmail.com',
      python_requires='>=3.8',
      package_data={'depmapomics': ['data/*']},
      )
