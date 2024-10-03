
OBJ = $(wildcard *.tex)

doc: ${OBJ}
	@xelatex $^
	@bibtex report.aux
	@xelatex $^
	@xelatex $^

.PHONY:clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg *-blx.bib

# Makefile

# 变量定义
CXX = g++
CXXFLAGS = -Wall -g
SOURCES = B.cpp C.cpp D.cpp E.cpp F.cpp
HEADERS = EquationSolver.hpp Function.hpp
EXECUTABLE = program
REPORT = report.tex
PDF = report.pdf

all: $(EXECUTABLE)


$(EXECUTABLE): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(SOURCES)

run: $(EXECUTABLE)
	./$(EXECUTABLE) > run.txt

report: $(REPORT)
	xelatex $(REPORT)

clean:
	rm -f $(EXECUTABLE) $(PDF) *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg *-blx.bib

.PHONY: all run report clean
