from setuptools import setup, find_packages
setup(
    name = 'larf',
    version = '0.0',
    package_dir = {'':'python'},
    packages = ['larf'],
    install_requires = [
        'Click',
    ],
    entry_points = dict(
        console_scripts = [
            'larf = larf.cli:main',
        ]
    )
)
