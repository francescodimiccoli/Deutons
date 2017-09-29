
#bdt: TrainBDT.cpp
#	g++ -o likelihood  TrainBDT.cpp `root-config --cflags --libs` -lTMVA
plotsBDT.pdf:
	rm -rf plots
	root -l -b -q "plottingScripts/variables.C(\"QualityTOF.root\")"
	root -l -b -q "plottingScripts/variables.C(\"QualityTOF.root\",\"InputVariables_Deco\", \"Decorellated Variables\")"
	root -l -b -q "plottingScripts/variables.C(\"QualityTOF.root\",\"InputVariables_PCA\", \"Prinipal components\")"
	root -l -b -q "plottingScripts/correlations.C(\"QualityTOF.root\")"
	root -l -b -q "plottingScripts/mvas.C(\"QualityTOF.root\")"
	root -l -b -q "plottingScripts/likelihoodrefs.C(\"QualityTOF.root\")"
	root -l -b -q "plottingScripts/efficiencies.C(\"QualityTOF.root\")"
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=plotsBDT.pdf `ls -hrt plots/*`
		





