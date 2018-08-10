from distutils.core import setup

packages_pyscat = ['pyscat', 
       'pyscat.calc',
       'pyscat.crystal',
       'pyscat.latdynam',
       'pyscat.utils'
       ]


setup(name='pyscat',
      version='0.1',
      description='calculator of phonon scattering due to an impurity with T-matrix',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=packages_pyscat,
      requires=['numpy', 'phonopy'],
      url='https://github.com/masato1122/pyscat.git',
      license='MIT',
      provides=['pyscat']
      )

