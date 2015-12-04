from distutils.core import setup
import ensemble_caller

setup(
    name='ensemble_caller',
    version=ensemble_caller.__version__,
    author='Bruno Grande',
    author_email='bgrande@sfu.ca',
    url="https://github.com/brunogrande/ensemble_caller",
    scripts=['ensemble_caller.py'],
    description=ensemble_caller.__desc__,
    long_description=open('README.md').read(),
    install_requires=[
        "pyvcf >= 0.6.7"
    ]
)
