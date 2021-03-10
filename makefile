
all: figure1.svg $(patsubst %.jmd, %.pdf, $(wildcard *.jmd))


%.pdf: %.jmd
	julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.precompile()'
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); ENV["GKSwstype"]="100"; using Weave; weave("$<", doctype = "md2pdf")'
	rm $*.log $*.out $*.tex $*.aux

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using DrugResponseModel; DrugResponseModel.figure$*()'

clean:
	rm -rf *.pdf
