from setuptools import find_packages, setup

app_requirements = [
    'loguru>=0.6',
    'matplotlib>=3.7.1',
    'numpy>=1.24.2',
    'pandas>=1.5.3, <2.0',
    'scipy>=1.10.1'
]

dev_requirements = [
    'pytest>=7.2',
]

setup(name='haem_kinetics',
      description='Haem Speciation Kinetics',
      long_description='Simulates the kinetics of haem speciation in the malaria parasite',
      version='0.1.0',
      url='https://github.com/davidkuter/haem_kinetics',
      license='MIT',
      packages=find_packages(),
      install_requires=app_requirements,
      extras_require={'dev': dev_requirements})
