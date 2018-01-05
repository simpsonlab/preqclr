from setuptools import setup

setup(name='preqc-lr',
	version='1.2',
	description='Quality control tool for long read sequencing data',
	url='https://github.com/simpsonlab/preqc-lr',
	author='Simpson Lab',
	author_email='joanna.pineda@oicr.on.ca',
	license='MIT',
	classifiers=[
    	'Development Status :: 3 - Alpha',
		'Topic :: Software Development :: Build Tools',
		'License :: OSI Approved :: MIT License',
    	'Programming Language :: Python :: 2.7',
	],
	python_requires='>=2.7',
	setup_requires=['numpy'],
	install_requires=[
		'Biopython',
        'gevent',
        'greenlet',
		'matplotlib==2.0.0'
	],
	include_package_data=True,
    zip_safe=False)
