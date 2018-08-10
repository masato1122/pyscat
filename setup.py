from distutils.core import setup

setup(name='pyscat',
      version='0.1',
      description='T-matrix method for phonon',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=['pyscat'],
      requires=['numpy', 'phonopy'],
      url='https://github.com/masato1122/pyscat',
      license='MIT',
      provides=['pyscat']
      )

