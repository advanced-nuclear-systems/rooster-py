# ROOSTER: Robust Object Oriented Solver of Transport Equations in a Reactor

The ROOSTER projects aims at developing a software for multi-physics modeling of Generation-IV Sodium Fast Reactor.

## How to install on Linux
1. Clone the code:
   `git clone https://github.com/advanced-nuclear-systems/rooster-py.git`
   or download and unpack ZIP.

2. Create virtual environment:
   `python -m  venv venv-rooster`

3. Activate it:
   `source venv-rooster/bin/activate`

4. Update pip
   `pip install --upgrade pip`

5. Install dependencies:
   `pip install CoolProp`
   `pip install pyinstaller`

6. Compile Fortran source to the `.so` library by running 
   `python _compile.py`

Notes: 
- The `gfortran` compiler should be installed.
- If you want to pack the source code into an executable code, uncomment corresponding lines in the `_compile.py`
- ROOSTER has not yet been tested for Windows.
- More details are at https://advanced-nuclear-systems.github.io/.

## How to execute on Linux
1. Launch ROOSTER
   `python rooster.py`.

2. Find the results in the `output` directory.

