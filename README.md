# ROOSTER: Robust Object Oriented Solver of Transport Equations in a Reactor

The ROOSTER projects aims at developing a software for multi-physics modeling of Generation-IV Sodium Fast Reactor.

## How to install and execute on Linux
1. Install dependencies:
- pip3 install CoolProp
- pip3 install pyinstaller

2. Clone the code: `git clone https://github.com/rooster-code/rooster.git or download and unpack ZIP.

3. Compile Fortran source to the `.so` library by running `compile` or `compile_noOMP` batch file. Note that `gfortran` compiler should be installed.

4. Launch ROOSTER by `python rooster.py`.

5. Find the results in the `output` directory.

6. if you want to package this source code into a executable code, please install pyinstaller and uncomment corresponding lines in the `compile`

ROOSTER has not yet been tested for Windows.

More details are at https://rooster-code.github.io/.
