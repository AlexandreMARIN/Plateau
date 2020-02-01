# Plateau
We solve here the discretised Plateau's problem.

We consider only surfaces given by functions in two real variables.


## In order to compile

You can type in your terminal:

```bash
make all
```

## Usage

To run the program, write:

```bash
./solvePlateau inputname tol maxnbiter outputname
```
The arguments are, in that order:
 - the name of a file with extension .mesh
 - the tolerance (a real number)
 - the maximum number of iterations (an integer)
 - the name of a file, which is supposed to be used by `gnuplot` later.

The number of used iterations will be displayed.

If you want to plot the results, run `gnuplot` and type:

```gnuplot
splot "outputname" with lines linewidth 0.1
```

## Structure of the projet

- sources are in `src` and headers are in `hdr`
- some meshes are in `input` : you can test them
- several directories will be created by running the Makefile.

## Meshes
[`FreeFEM++`]: https://freefem.org/ "FreeFEM++"

Meshes have been created with [`FreeFEM++`].

See in the associated documentation to know what .mesh files look like.
