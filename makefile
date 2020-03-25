
all: allDrugs.html combination.html notebookPlots.pdf

venv: venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uq jupyter
	. venv/bin/activate && julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.add("Weave"); Pkg.add("Coverage"); Pkg.add("Plots"); Pkg.precompile()'
	touch venv/bin/activate

%.pdf: %.jmd venv
	julia -e 'using Pkg; Pkg.activate("."); using Weave; weave("$<", doctype = "md2pdf")'

coverage.cob: venv
	julia -e 'using Coverage; coverage = process_folder(); LCOV.writefile("coverage-lcov.info", coverage)'
	pip3 install --user lcov_cobertura
	python3 ~/.local/lib/python3.7/site-packages/lcov_cobertura.py coverage-lcov.info -o coverage.cob

%.html: %.ipynb venv
	. venv/bin/activate && jupyter nbconvert --execute --ExecutePreprocessor.timeout=60000 --to html $< --output $@

clean:
	rm -rf *.html *.log *.pdf *.log *.aux *.out
