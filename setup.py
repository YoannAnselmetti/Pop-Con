from setuptools import setup

setup(
   name='Pop-Con',
   version='1.0.1',
   packages=['popcon'],
   classifiers=[
      "Development Status :: 4 - Beta",
      "Environment :: Console",
      "Intended Audience :: Science/Research",
      "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
      "Natural Language :: English",
      "Operating System :: POSIX :: Linux",
      "Programming Language :: Python",
      "Topic :: Scientific/Engineering :: Bio-Informatics",
      "Topic :: Scientific/Engineering :: Visualization"
   ],
   scripts=['Pop-Con','scripts/'],
   description='Visualization of genotype profiles on population genomics data for detection of abnormal genotypes pattern.',
   long_description_content_type="text/markdown",
   long_description=open("README.md").read(),
   author='Yoann Anselmetti',
   author_email='yoann.anselmetti@umontpellier.fr',
   url='https://github.com/YoannAnselmetti/Pop-Con',
   license='CeCILL-2.1',
   install_requires=[
      'matplotlib==2.2.4',
      'numpy==1.16.4',
      'cyvcf2==0.11.4'
   ],
)