
#bdt: TrainBDT.cpp
#	g++ -o likelihood  TrainBDT.cpp `root-config --cflags --libs` -lTMVA
plotsBDT.pdf:
	rm -rf plots
	root -l -b -q "plottingScripts/variables.C(\"QualityNaF.root\")"
	root -l -b -q "plottingScripts/variables.C(\"QualityNaF.root\",\"InputVariables_Deco\", \"Decorellated Variables\")"
	root -l -b -q "plottingScripts/variables.C(\"QualityNaF.root\",\"InputVariables_PCA\", \"Prinipal components\")"
	root -l -b -q "plottingScripts/correlations.C(\"QualityNaF.root\")"
	root -l -b -q "plottingScripts/mvas.C(\"QualityNaF.root\")"
	root -l -b -q "plottingScripts/likelihoodrefs.C(\"QualityNaF.root\")"
	root -l -b -q "plottingScripts/efficiencies.C(\"QualityNaF.root\")"
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=plotsBDT.pdf `ls -hrt plots/*`
		





