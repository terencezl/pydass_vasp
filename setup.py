from setuptools import setup
long_description = ''
try:
   import pypandoc
   long_description = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError):
   long_description = open('README.md').read()


setup(name='pydass_vasp',
      version='0.1',
      description='Convenient Python modules and wrapping script executables.',
      long_description=long_description,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      keywords='vasp plotting fitting manipulation',
      url='https://github.com/terencezl/pydass_vasp',
      author='Terence Zhi Liu',
      author_email='zhi.liu@utoledo.edu',
      license='MIT',
      packages=['pydass_vasp', 'pydass_vasp.plotting', 'pydass_vasp.fitting', 'pydass_vasp.manipulation'],
      scripts=['scripts/plot_tdos.py', 'scripts/plot_ldos.py', 'scripts/plot_cohp.py', 'scripts/plot_bs.py'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib >= 1.3'
      ],
      include_package_data=True,
      zip_safe=False)