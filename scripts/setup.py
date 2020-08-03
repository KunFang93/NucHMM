from setuptools import setup

setup(
    name='NucHMM',
    version='1.0',
    py_modules=['NucHMM'],
    install_requires=[
        'Click','scipy','numpy','matplotlib','pandas','seaborn','statsmodels'
    ],
    entry_points='''
    [console_scripts]
    NucHMM=NucHMM:cli
    '''
)
