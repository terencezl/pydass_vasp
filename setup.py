import re
from setuptools import setup


def remove_extra_line(text):
    text = text.split('\n')
    i = 0
    while i < len(text):
        m = re.match(r'\s*:alt:\s*(\w*).*', text[i])
        if m:
            alt_text = m.group(1)
            i += 1
            while i < len(text):
                if alt_text in text[i]:
                    text.pop(i)
                    break
                i += 1
        i += 1
    return '\n'.join(text)

long_description = ''
try:
   import pypandoc
   long_description = pypandoc.convert('README.md', 'rst')
   long_description = remove_extra_line(long_description)
except (IOError, ImportError):
   long_description = open('README.md').read()

setup(name='pydass_vasp',
      version='0.1',
      description='Convenient Python modules and wrapping script executables',
      long_description=long_description,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      keywords='vasp plotting fitting manipulation',
      url='http://terencezl.github.io/pydass_vasp/',
      author='Terence Zhi Liu',
      author_email='zhi.liu@utoledo.edu',
      license='MIT',
      packages=['pydass_vasp', 'pydass_vasp.electronic_structure', 'pydass_vasp.fitting'],
      scripts=['scripts/plot_tdos.py', 'scripts/plot_ldos.py', 'scripts/plot_lobster.py', 'scripts/plot_bs.py'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'matplotlib >= 1.3'
      ],
      include_package_data=True,
      # zip_safe=False
      )