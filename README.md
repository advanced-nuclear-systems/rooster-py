# ROOSTER: Robust Object Oriented Solver of Transport Equations in a Reactor

The ROOSTER projects aims at developing a software for multi-physics modeling of Generation-IV Sodium Fast Reactor.

## How to install and execute on Linux
1. Clone the code: `git clone https://github.com/armstrong-dev/rooster.git` or download and unpack ZIP.

2. Compile Fortran source to the `.so` library by running `compile` or `compile_noOMP` batch file. Note that `gfortran` compiler should be installed.

3. Launch ROOSTER by entering `python3 A_rooster.py`.

4. Find the results in the `output` directory.

5. if you want to package this source code into a executable code, please install pyinstaller and uncomment corresponding lines in the `compile`

ROOSTER has not yet been tested for Windows.

More details are at https://armstrong-dev.github.io/index.html#rooster.

## How to push to armstrong-dev
git push rooster main
