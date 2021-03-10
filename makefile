
all: $(patsubst %.jmd, %.pdf, $(wildcard *.jmd))

%.pdf: %.jmd
	julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.precompile()'
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); ENV["GKSwstype"]="100"; using Weave; weave("$<", doctype = "md2pdf")'
	rm $*.log $*.out $*.tex $*.aux


venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uqr requirements.txt
	touch venv/bin/activate

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using DrugResponseModel; DrugResponseModel.figure$*()'


clean:
	rm -rf *.pdf
