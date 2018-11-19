float GetW(TH1F * integral_cutoff, float Rmeas, float Rgen){
	
	float Rmin;
	if(Rgen<Rmeas) Rmin = Rgen;
	else Rmin = Rmeas;

	
	return integral_cutoff->GetBinContent(integral_cutoff->FindBin(Rmin));
}
