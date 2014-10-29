Casimir
=======

Casimir Effect in the plane-sphere geometry with perfect reflectors.

There is also a website available: http://speicherleck.de/michael/casimir.html
You can cite this software, see: https://zenodo.org/record/12476


Description
===========

The progams in src/ implement the numerics for the Casimir effect in the
plane-sphere geometry with perfect spheres using a scattering approach [1,2].
We use the same approach, but derive slightly different formulas for the matrix
element of the round-trip operator that don't need Wigned-d-symbols.

A sphere of radius R is separated by a distance of L from a plane. The plane is
infinite in the xy-direction and both plane and sphere are assumed to be
perfect reflectors. The programs calculate the free energy F(T,R/L) in scaled
quantities.

This code is part of my master thesis. It costed much time and much work to
write a working and a fast implementation. So, if you find this piece of code
useful and you use it for plots, please consider to cite my work.

The directory talk/ contains a presentation about the Casimir effect in the
plane-sphere geometry in Germand language. The directory report/ contains a
report about the derivation of the matrix elements of the round-trip matrix,
also in German language.


Installation
============

If you use Linux or Unix, you need the gcc and development libraries and header
files for the standard C library. On a Debian-like Linux the command

    $ apt-get install gcc libc6-dev

should install all dependencies. After that, change to the directory src/
and run

    $ make

to compile the sources. Every program prints out a short usage message when
called with the switch -h.


Programs
========

 * casimir

   Calculate the free energy F(T, L/R) in scaled quantities.

 * casimir_hiT

   Calculate the free energy F(T->oo, L/R) for very high temperaturs.

 * casimir_logdetD

  Calculate log(det(D)) for given T, L/R, n and m

 * casimir_tests

   Run a unit test.

 * pfa.py

   Implementation of the PFA (proximity force approximation)

 * Casimir.py

   Python wrapper for casimir


License
=======

The software is licensed as GPLv2, see the file LICENSE. However, if you use
this program for publications, please consider to cite my work.

Bibliography
===========

[1] Michael Hartmann, "Negative Casimir entropies in the plane-sphere geometry", see thesis.pdf

[2] Antoine Canaguier-Durand et al., "Thermal Casimir Effect in the plane-sphere geometry", Phys. Rev. A 82 (1 2010)
