CXX = g++
CXXFLAGS = -Wall -g
SOURCES = B.cpp C.cpp D.cpp E.cpp F.cpp
HEADERS = EquationSolver.hpp Function.hpp
EXECUTABLES = B C D E F
REPORT = report.tex
PDF = report.pdf


all: run report

$(EXECUTABLES): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<


run: $(EXECUTABLES)
	@rm -f run.txt;                      
	@for exe in $(EXECUTABLES); do \
		./$$exe >> run.txt; \
	done


report: $(REPORT)
	xelatex $(REPORT)

clean:
	rm -f $(EXECUTABLES) $(PDF) *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg *-blx.bib run.txt 

.PHONY: all run report clean
