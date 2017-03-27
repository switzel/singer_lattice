singer_lattice
==============

This repository contains code relevant to the article [1]. To produce documentation type

    make documentation

The bulk of the code is written in Python 2 [2]. For the stabilizer computations you will also need GAP [3] with GRAPE [4] installed. The code has only been tested on Unix environments and will probably need adjustments on other systems.

The main script is code/hjelmslev.py. Try

    cd code
    python difference_matrices.py 2 --write_lattices
    python hjelmslev.py 2 --moufang --stabilizers --write_lattices
    python summarize.py 2 --properties diff_mat based_diff_mat aut_dm abelianization commutator_abelianization moufang_2 split_1_2 K_1^2 "N/C^2(1)"

or

    python hjelmslev.py --help

for more information.

[1]: Stefan Witzel, On panel-regular ~A_2-lattices, arXiv:1608.07141.  
[2]: Python, http://www.python.org  
[3]: GAP, http://www.gap-system.org  
[4]: Leonard H. Soicher, GRAPE, http://www.maths.qmul.ac.uk/~leonard/grape/
