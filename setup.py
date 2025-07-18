from setuptools import setup, find_packages

setup(
    name='igm_model',
    version='0.1',    
    description='Code to calculate the absorption due to foreground H I in the intergalactic medium, based on the models of Kauma et al. (in prep)',
    url='https://github.com/katkauma/igm_model',
    packages=find_packages(),
    install_requires=['numpy',
                      'scipy',                
                      ],
)
