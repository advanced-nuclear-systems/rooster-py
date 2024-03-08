#This Python script:
#    1) compiles Fortran files
#    2) generates executable
#    3) creates in /dist directory a timestamped folder and copies there timestamped executable and timestamped zip archive of all source files
#Note:
#    If no Open Multi-Processing is available, set NO_OMP to True.

import datetime
import glob
import os
import subprocess

from zipfile import ZipFile

NO_OMP = False

#timestamp name
name = 'rooster-' + str(datetime.datetime.now())[0:19].replace(" ","-").replace(":","-")

#list of fortran files
filesf90 = glob.glob('./*.f90')
#compile all fortran files
for filef90 in filesf90:
    if NO_OMP:
        subprocess.run(["python3", "-m", "numpy.f2py",           "-c",                        "--build-dir", "_f2py", filef90[2:], "-m", filef90[2:-4]])
    else:
        subprocess.run(["python3", "-m", "numpy.f2py", "-lgomp", "-c", "--f90flags=-fopenmp", "--build-dir", "_f2py", filef90[2:], "-m", filef90[2:-4]])

#list of library files
filesso = glob.glob('./*.so')
#rename library to make name shorter
for fileso in filesso:
    split = fileso[2:].split('.')
    os.rename(fileso[2:], split[0] + ".so")

#generate executable
subprocess.run(["pyinstaller", "-F", "rooster.py"])

#make directory with timestamp in /dist directory
os.mkdir("./dist/" + name)
#add timestamp to the executable name and move it in the timestamp directory
os.rename("./dist/rooster", "./dist/" + name + "/" + name)

#list of python files
filespy = glob.glob('./*.py')

#zip all source files
with ZipFile('./dist/' + name + '/' + name +'.zip', 'w') as zip_object:
    for filepy in filespy:
        zip_object.write(filepy)
        print("add to zip archive " + filepy)
    for filef90 in filesf90:
        zip_object.write(filef90)
        print("add to zip archive " + filef90)
    file = "LICENSE"
    zip_object.write(file)
    print("add to zip archive " + file)
    file = "README.md"
    zip_object.write(file)
    print("add to zip archive " + file)
