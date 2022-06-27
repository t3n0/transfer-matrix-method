from setuptools import setup, find_packages

setup(name='tmm-tool',
      version='0.1',
      description='Tranfer Matrix Method tool (TMM-tool)',
      author='Stefano Dal Forno',
      author_email='tenobaldi@gmail.com',
      url='https://github.com/t3n0/transfer-matrix-method',
      classifiers=[
          'Development Status :: 4 - Beta',

          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Physics',

          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
      ],
      keywords=['transfer', 'matrix'],
      #packages=find_packages(),
      packages = ['tmm'],
      python_requires=">=3.5",
      install_requires=[],
      entry_points={ "console_scripts": [ "tmm-tool=tmm.main:main" ],},
      )
