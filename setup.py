import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DEL-Analysis",
    version="0.0.1",
    author="Rory Coffey",
    author_email="coffeyrt@gmail.com",
    description="Analysis algorithms that are a compantion to DEL-Decode",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Roco-scientist/DEL-Analysis",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
