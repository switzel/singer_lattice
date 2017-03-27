PROJECT=singer_lattice
INSTANCES4=0 4 1 2 7 11 9 16 5 14 15 3 6 8 10 12 13
INSTANCES5=0 22 6 17 4 23 20 9 2 1 21 14 11 15 16 18 19 7 13 12 8 24 10 5 3
LATTICE_ARGS=--number --number_format "\\Gamma_{%d, %d}" --latex
BUILDING_ARGS=--number --number_format "X_{%d, %d}" --latex
LATTICE_PROPS2=diff_mat based_diff_mat aut_dm abelianization commutator_abelianization
LATTICE_PROPS3=diff_mat based_diff_mat aut_dm abelianization commutator_abelianization
LATTICE_PROPS4=diff_mat based_diff_mat aut_dm abelianization commutator_abelianization
LATTICE_PROPS5=num based_diff_mat aut_dm abelianization commutator_abelianization
BUILDING_PROPS2=based_diff_mat moufang_2 moufang_3 split_1_2 K_1^2 "K_1^3(2)" "N/C^2(1)"
BUILDING_PROPS3=based_diff_mat moufang_2 moufang_3 split_1_2 K_1^2 "K_1^3" "N/C^2(1)"
BUILDING_PROPS4=based_diff_mat moufang_2 split_1_2 K_1^2 "K_1^3(2)" "N/C^2(1)"
BUILDING_PROPS5=num based_diff_mat moufang_2 split_1_2 K_1^2 "N/C^2(1)"
DIGEST=--matrices_filename "digest_mats_%d.pickled" --lattice_filename "digest_%d.pickled"

all: article documentation

article: $(PROJECT).pdf

documentation: documentation.pdf

figs: $(PROJECT)-0.mps hjelmslev-0.mps

tabs: lattices_2.tex lattices_3.tex lattices_4.tex lattices_5.tex buildings_2.tex buildings_3.tex buildings_4.tex buildings_5.tex

code: difference_matrices.py hjelmslev.py hjelmslev.gap stabilizers.gap

$(PROJECT).pdf: $(PROJECT).tex $(PROJECT).bbl figs tabs
	pdflatex $(PROJECT).tex
	pdflatex $(PROJECT).tex

$(PROJECT).bbl: $(PROJECT).bib $(PROJECT).aux
	bibtex $(PROJECT)

$(PROJECT).aux: $(PROJECT).tex figs tabs code
	pdflatex $(PROJECT).tex

%-0.mps: %.mp
	mpost $<

lattices_2.tex: code/summarize.py code/lattices_2.pickled
	cd code ;\
	python summarize.py 2 $(LATTICE_ARGS) --properties $(LATTICE_PROPS2) > ../$@ ;\
	cd ..

lattices_3.tex: code/summarize.py code/lattices_3.pickled
	cd code ;\
	python summarize.py 3 $(LATTICE_ARGS) --properties $(LATTICE_PROPS3) > ../$@ ;\
	cd ..

lattices_4.tex: code/summarize.py code/lattices_4.pickled
	cd code ;\
	python summarize.py 4 $(LATTICE_ARGS) --properties $(LATTICE_PROPS4) --instances $(INSTANCES4) > ../$@ ;\
	cd ..

lattices_5.tex: code/summarize.py code/digest_5.pickled
	cd code ;\
	python summarize.py 5 --latex $(DIGEST) --properties $(LATTICE_PROPS5)  --instances $(INSTANCES5) > ../$@ ;\
	cd ..

buildings_2.tex: code/summarize.py code/lattices_2.pickled
	cd code ;\
	python summarize.py 2 $(BUILDING_ARGS) --properties $(BUILDING_PROPS2) > ../$@ ;\
	cd ..

buildings_3.tex: code/summarize.py code/lattices_3.pickled
	cd code ;\
	python summarize.py 3 $(BUILDING_ARGS) --properties $(BUILDING_PROPS3) > ../$@ ;\
	cd ..

buildings_4.tex: code/summarize.py code/lattices_4.pickled
	cd code ;\
	python summarize.py 4 $(BUILDING_ARGS) --properties $(BUILDING_PROPS4) --instances $(INSTANCES4) > ../$@ ;\
	cd ..

buildings_5.tex: code/summarize.py code/digest_5.pickled
	cd code ;\
	python summarize.py 5 --latex $(DIGEST) --properties $(BUILDING_PROPS5) --instances $(INSTANCES5) > ../$@ ;\
	cd ..

code/digest_5.pickled: code/digest.py code/lattices_5.pickled
	cd code ;\
	python digest.py 5

%.py: code/%.py code/%.awk
	gawk -f $(word 2, $^) $< > $@

%.gap: code/%.gap
	cp $< $@

documentation.pdf: documentation.tex documentation.bbl code
	pdflatex documentation.tex
	pdflatex documentation.tex

documentation.bbl: documentation.tex $(PROJECT).bib
	pdflatex documentation.tex

clean:
	rm -f *.aux *.bbl *.log *.blg *.mpx *.out ;\
	rm -f *.mps *.py *.gap ;\
	rm -f lattices_2.tex lattices_3.tex lattices_4.tex lattices_5.tex buildings_2.tex buildings_3.tex buildings_4.tex buildings_5.tex
	rm -f code/digest_*pickled
