# HEAD
Hahn-Echo Assisted Deconvolution

Hahn-Echo Assisted Deconvolution (HEAD)

This repository contains the Bruker pulse program for acquire 2D HEAD data in addition to the C++ program used for data processing. Data must be acquired with identical digital resolution in both dimensions. The  resulting 2D spectrum is converted to a ASCII file using the Topspin 'totxt' command. This file can be handled by the HEAD_processing program.

A recent addition is a separate program used to enhance the resolution of a 2D spectrum using a separate 2D HEAD spectrum. Both 2D spectra must have the same F2 spectral widths and number of datapoints, with the 2D HEAD spectrum following the requirements outlined above. This program, including binary and examples, are located in the 2D-HEAD sub-directory.

A binary version of this program is provided for Windows computers. To compile it for other operating systems you must link the GNU scientific library (https://www.gnu.org/software/gsl/doc/html/) and compile as:

g++ HEAD_processing.cpp -o HEAD_processing -lgsl -lgslcblas -lm -O3 -march=native

Copyright 2024. Iowa State University. This material was produced under U.S. Government contract DE-AC02-07CH11358 for the Ames National Laboratory, which is operated by Iowa State University for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software. NEITHER THE GOVERNMENT, AMES NATIONAL LABORATORY, NOR IOWA STATE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from Ames National Laboratory.

Additionally, this program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. (https://www.gnu.org/licenses/gpl-3.0.en.html#license-text)icense as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. (https://www.gnu.org/licenses/gpl-3.0.en.html#license-text)
