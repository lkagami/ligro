from setuptools import setup, find_packages, Extension

setup(
    name='ligro',
    packages=['ligro'],
    version='1.1.2',
    description='LiGRO: a graphical user interface for protein-ligand molecular dynamics.',
    author="KAGAMI, Luciano Porto",
    author_email='luciano.kagami@ufrgs.br',
    url='https://github.com/lkagami/ligro',
    keywords=['ligro'],
    include_package_data=True,
    entry_points={'console_scripts': ['ligro = ligro.ligro:main_init']},
    zip_safe=False,
    classifiers=[
        '',
    ]
)
