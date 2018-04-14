# EWALDOPT - Minimisation of Ewald sums

*EWALDOPT* minimises the [Ewald sum](https://en.wikipedia.org/wiki/Ewald_summation) over a configuration of point charges. It maintains overall charge neutrality and expects integer charges on all atomic sites.

## Algorithm

By default, the program uses a brute-force algorithm for small unit cells (less than ten atoms) and switches to a Monte Carlo algorithm for larger cells as otherwise
computational cost would be prohibitive. The Ewald summation follows the classical paper from Ewald.

## Usage

The program reads input from `stdin` and prints result to `stdout`. Errors are reported to `stderr`.

Use pipes to feed files to the program. For instance, to minimise the structure stored in `tests/test.str`, use

```
cat tests/test.str | ewaldopt
```

We also provide a micro-container. To pipe `stdin` using Docker, use

```
cat tests/test.str | docker run -i kramergroup/ewaldopt
```

*Note*: The `-i` flag is needed for pipes to work even if this is not an interactive use case. See the [Docker docs](https://docs.docker.com/engine/reference/run/#foreground) for details.

### Input

The program reads a structure file from `stdin` and prints a similar structure file to `stdout` with assigned minimal valences. Examples can be found in *tests*.

An example input file for a hypothetical SnF could look like this:

```
5.6402 0.0000 0.0000
0.0000 5.6402 0.0000
0.0000 0.0000 5.6402
8
0.0000 0.0000 0.5000 Sn  1.0 2.0 3.0
0.0000 0.0000 0.0000 F  -1.0
0.5000 0.0000 0.0000 Sn  1.0
0.0000 0.5000 0.0000 Sn  1.0
0.0000 0.5000 0.5000 F  -1.0 -2.0
0.5000 0.0000 0.5000 F  -1.0
0.5000 0.5000 0.5000 Sn  1.0
0.5000 0.5000 0.0000 F  -1.0 -2.0
```

The structure of the input file is as follows:

- The first three lines define the unit cell with one vector per line
- The next line contains the total number *N* of atomic sites in the unit cell
- The following *N* lines describe the position, species, and possible charge of each site

### Output

The returned structure is reported in a similar format with the assigned minimal
valence shown in the last column:

```
0.5640200E+01   0.0000000E+00   0.0000000E+00
0.0000000E+00   0.5640200E+01   0.0000000E+00
0.0000000E+00   0.0000000E+00   0.5640200E+01
8
0.0000000E+00   0.0000000E+00   0.5000000E+00   Sn   0.00
0.0000000E+00   0.0000000E+00   0.0000000E+00    F  -1.00
0.5000000E+00   0.0000000E+00   0.0000000E+00   Sn   3.00
0.0000000E+00   0.5000000E+00   0.0000000E+00   Sn   0.00
0.0000000E+00   0.5000000E+00   0.5000000E+00    F  -1.00
0.5000000E+00   0.0000000E+00   0.5000000E+00    F  -2.00
0.5000000E+00   0.5000000E+00   0.5000000E+00   Sn   3.00
0.5000000E+00   0.5000000E+00   0.0000000E+00    F  -2.00
```

### Parameters

There is a small number of command-line parameters to change the behaviour of the program.

```
ewaldopt - Finds valence distribution in variable valence compounds
           with minimal ewald sum"
Options:
-c    Allow charged cells
-l    Write logging info into file ewaldopt.log
-h    Display this help message
-d    Write debugging info to stdout
-b    Use brute-force regardless of problem size. Beware: this can
      become very expensive very quickly
-v    Display program version
```

## Compile and install

The provided `makefile` uses `gfortran` and it should be as simple as

```
make
```

There are a set recommended of compiler flags for Intel's Fortran compiler `ifort`
as comments in the makefile as well.
