CXX = g++
CXXFLAGS = -Wall -g
SOURCES = B.cpp C.cpp D.cpp E.cpp F.cpp
HEADERS = Interpolation.hpp Function.hpp
EXECUTABLES = B C D E F
REPORT = report.tex
PDF = report.pdf

# Doxygen 配置文件
DOXYFILE = Doxyfile

all: run report doxygen

$(EXECUTABLES): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

run: $(EXECUTABLES)
	@rm -f run.txt;                      
	@for exe in $(EXECUTABLES); do \
		./$$exe >> run.txt; \
	done

report: $(REPORT)
	xelatex $(REPORT)
	xelatex $(REPORT)

doxygen:
	@echo "Generating Doxygen documentation..."
	@doxygen $(DOXYFILE)

clean:
	rm -f $(EXECUTABLES) $(PDF) *.aux *.log *.out *.synctex.gz *.toc *.bbl *.run.xml *.blg *-blx.bib run.txt
	rm -rf html latex  # Clean up Doxygen generated documentation

.PHONY: all run report clean doxygen
