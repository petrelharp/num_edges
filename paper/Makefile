writeup.pdf : references.bib

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done
	touch $*.bbl
	touch $@

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux references.bib
	bibtex $<

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.svg
	inkscape $< --export-area-drawing --export-filename=$@
	# chromium --headless --no-pdf-header-footer --print-to-pdf=$@ $<
	# ./svg2pdf.sh $< $@

%.pdf : %.eps
	# inkscape $< --export-filename=$@
	epspdf $<

%.pdf : %.ink.svg
	inkscape $< --export-filename=$@

writeup-diff%.tex : writeup.tex
	latexdiff-git -r $* writeup.tex

diff-to-submitted-version.pdf : writeup-diff654179b2c9eb5788c17bc86aa855cb366835be7c.pdf
	cp $< $@
