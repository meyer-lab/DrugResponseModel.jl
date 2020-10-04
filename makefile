
all: combination.pdf  replicatesAtOnce.pdf avgRepsAllDrugs.pdf temporal_combination.pdf
# separateDrugsAvg.pdf notebookPlots.pdf
venv: venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uq lcov_cobertura
	. venv/bin/activate && julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.add("Coverage"); Pkg.add("Plots"); Pkg.precompile()'
	touch venv/bin/activate

%.pdf: %.jmd venv
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); ENV["GKSwstype"]="100"; using Weave; weave("$<", doctype = "md2pdf")'

coverage.cob: venv
	julia -e 'using Pkg; using Coverage; Pkg.activate("."); Pkg.test("DrugResponseModel"; coverage=true); coverage = process_folder(); LCOV.writefile("coverage-lcov.info", coverage)'
	. venv/bin/activate && python3 venv/lib/python3.8/site-packages/lcov_cobertura.py coverage-lcov.info -o coverage.cob

clean:
	rm -rf *.html *.log *.log *.aux *.out venv
