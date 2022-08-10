from setuptools import setup

setup(name='fools',
      version='0.1.0',
      packages=['fools'],
      description='CLI tools for FISH probe generation',
      author='Adrian Tschan',
      author_email='adrian.tschan@uzh.ch',
      url='https://github.com/adrtsc/fools',
      entry_points={'console_scripts': ['generate_probes=fools.generate_probes:main',
                                        'generate_fasta=fools.generate_fasta:main']}
      )