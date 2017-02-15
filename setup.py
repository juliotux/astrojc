from distutils.core import setup

setup(name='astrojc',
      version='1.0',
      description='General astronomy codes used in my day-to-day',
      author='Julio Campagnolo',
      author_email='juliocampagnolo@gmail.com',
      url='https://www.on.br/',
      packages=['astrojc', 'astrojc', 'astrojc.fitting', 'astrojc.io',
                'astrojc.reducing', 'astrojc.photometry', 'astrojc.spectra',
                'astrojc.gui', 'astrojc.math']
     )
