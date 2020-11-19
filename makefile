
all: $(patsubst %.jmd, %.pdf, $(wildcard *.jmd))

%.pdf: %.jmd
	julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.precompile()'
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); ENV["GKSwstype"]="100"; using Weave; weave("$<", doctype = "md2pdf")'
	rm $*.log $*.out $*.tex $*.aux

clean:
	rm -rf *.pdf
