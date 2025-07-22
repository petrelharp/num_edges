writeup.pdf : references.bib review-responses2.tex

resubmission_version.pdf : writeup.pdf
	pdfjam --outfile $@ $< 1-20

cover_letter.pdf : writeup.pdf
	pdfjam --outfile $@ $< 21

review_responses.pdf : writeup.pdf
	pdfjam --outfile $@ $< 22-

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

writeup-diff%.tex : writeup.tex review-responses2.tex
	latexdiff-git --force -r $* writeup.tex

diff-to-submitted-version.pdf : writeup-diff800eb7a859e4173f87b569118e6720add966a846.pdf
	cp $< $@
