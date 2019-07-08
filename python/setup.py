from setuptools import setup

def readme():
    with open('../README.md') as f:
        return f.read()

setup(name             = 'cahnhilliard_2d',
      version          = '0.1',
      author           = 'Anthony M. DeGennaro',
      author_email     = 'adegennaro@bnl.gov',
      description      = 'Python spectral solver for 2D Cahn-Hilliard',
      long_description = readme(),
      classifiers      = [
        'Topic :: Materials Dynamics :: Cahn-Hilliard',
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: ISC License',
        'Programming Language :: Python :: 3.6',
      ],
      keywords         = 'cahn-hilliard materials dynamics',
      url              = 'http://github.com/adegenna/cahnhilliard_2d',
      license          = 'ISC',
      packages         = ['cahnhilliard_2d','cahnhilliard_2d.src','cahnhilliard_2d.tests'],
      package_dir      = {'cahnhilliard_2d'         : 'cahnhilliard_2d' , \
                          'cahnhilliard_2d.src'     : 'cahnhilliard_2d/src' , \
                          'cahnhilliard_2d.tests'   : 'cahnhilliard_2d/tests' },
      test_suite       = 'cahnhilliard_2d.tests',
      entry_points     = { 'console_scripts': ['Package = cahnhilliard_2d.tests.driver:main' ] },
      install_requires = [ 'numpy', 'scipy', 'matplotlib' ],
      python_requires  = '>=3',
      zip_safe         = False
)
