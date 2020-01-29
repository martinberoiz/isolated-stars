from setuptools import setup

# Get the version from astroalign file itself (not imported)
# with open('astroalign.py', 'r') as f:
#     for line in f:
#         if line.startswith('__version__'):
#             _, _, aa_version = line.replace("'", '').split()
#             break

with open('README.md', 'r') as f:
    long_description = f.read()

setup(name='tdagcodechallenge',
      version='1.0',
      description='Time Domain Astronomy Group Code Challenge',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Martin Beroiz',
      author_email='martin.beroiz@utrgv.edu',
      url='https://github.com/ctmobservatory/codechallenge',
      py_modules=['stars', ],
      install_requires=["scipy",],
      test_suite='tests',
      )
