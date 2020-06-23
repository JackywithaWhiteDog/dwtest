import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dwtest", # Replace with your own username
    version="0.0.1",
    author="JackywithaWhiteDog",
    author_email="jackyliu0129@gmail.com",
    description="Durbin-Watson Test",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JackywithaWhiteDog/dwtest.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)