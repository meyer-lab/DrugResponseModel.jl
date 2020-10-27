
all: time_combin.pdf Bliss.pdf notebookPlots.pdf temporal_combination.pdf replicatesAtOnce.pdf combination.pdf avgRepsAllDrugs.pdf separateDrugsAvg.pdf

%.pdf: %.jmd
	julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.precompile()'
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); ENV["GKSwstype"]="100"; using Weave; weave("$<", doctype = "md2pdf")'
	rm $*.log $*.out $*.tex $*.aux

clean:
	rm -rf *.pdf
