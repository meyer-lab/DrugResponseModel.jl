
venv: venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uq jupyter
	touch venv/bin/activate

coverage.cob:
	julia -e 'using Pkg; Pkg.add("Coverage"); using Coverage; coverage = process_folder(); LCOV.writefile("coverage-lcov.info", coverage)'
	pip3 install --user lcov_cobertura
	python3 ~/.local/lib/python3.7/site-packages/lcov_cobertura.py coverage-lcov.info -o coverage.cob

%.html: %.ipynb venv
	. venv/bin/activate && julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.precompile()'
	. venv/bin/activate && jupyter nbconvert --execute --ExecutePreprocessor.timeout=60000 --to html $< --output $@

clean:
	rm -rf *.pdf *.aux *.log *.out
