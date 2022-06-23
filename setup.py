from setuptools import setup

setup(
    name='ArachnoFamTox',
    version='1.0.0',
    author = "Fernanda Midori Abukawa",
    author_email = "fernanda.abukawa@gmail.com",
    license = 'GNU',
    packages = ['arachnofamtox'],
    install_requires=[
        'biopython == 1.79',
        'numpy == 1.20.3',
        'pandas == 1.4.1'],
    entry_points={
        'console_scripts': [
            'ArachnoFamTox = arachnofamtox.main:run',
        ]
    }
)
