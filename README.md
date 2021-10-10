Image Editing
=============

Table of Contents
-----------------
* [TODO](#todo)
* [Requirement](#requirement)
* [Build](#build)
* [Run](#run)

TODO
----
- [x] Transform	ToGray               (5)
- [x] Quarantined Uniform            (5)
- [x] Quarantined Populosity        (20)
- [x] Dithering	Naive                (3)
- [x] Dithering	Brightness           (7)
- [x] Dithering	Random               (5)
- [x] Dithering	Cluster             (10)
- [x] Dithering	Floyd               (15)
- [x] Dithering	Color Floyd         (10)
- [x] Filter Box                  (15/3)
- [x] Filter Barlette             (15/3)
- [x] Filter Gaussian             (15/3)
- [x] Filter Ab G                   (10)
- [x] Filter Edge Detect          (15/3)
- [x] Filter Edge Enhance         (15/3)
- [x] Resizing Half                  (8)
- [x] Resizing Double               (12)
- [x] Resizing Arbitrary Size    (25/10)
- [x] Resizing Arbitrary Rotate  (25/10)
- [ ] NPR                        (15-50)

Requirement
-----------
Please install the following programs before building. 

* [MinGW]
* [FLTK]
* [CMake]

[MinGW]: https://osdn.net/projects/mingw/
[FLTK]: https://www.fltk.org/
[CMake]: https://cmake.org/

Build
-----
Run `cmake -G "MinGW Makefiles -B build ."` first,

and run `make all` in `build` folder.

Run
---
All executable files can be run both in console or in window form.

You can run `<program> --help` to see more information.
