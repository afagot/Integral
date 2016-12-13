#include "../include/fonctions.h"

using namespace std;

int main(int argc ,char *argv[]){
	stringstream converter;
	converter << argv[0];
	string program;
	converter >> program;
	converter.clear();

	if(argc != 2){
			cout << "[Integrals] expects to have 2 parameters\n";
			cout << "[Offline] USAGE is : " << program << " filebasename\n";
			return -1;
		} else if(argc == 2){
			converter << argv[1];
			string baseName;
			converter >> baseName;
			converter.clear();

			//intput ROOT file
			string ROOTName = baseName + ".root";
			TFile ROOTinput(ROOTName.c_str(),"READ");

			for(int tr=0; tr<5; tr++){
				float RPCthreshold = 2. + tr*0.25;

				//output csv file with integrals results
				string CSVName = baseName.substr(baseName.find_last_of("/")+1) + "-" + floatTostring(RPCthreshold) + ".csv";
				ofstream CSVFile(CSVName.c_str(), ios::out);

				for(int ts=0; ts<=45; ts+=5){

					string FProfName = "prof_f-" + intTostring(ts) + "-" + floatTostring(RPCthreshold);
					TH1I* FProfile = (TH1I*)ROOTinput.Get(FProfName.c_str());

					string BProfName = "prof_b-" + intTostring(ts) + "-" + floatTostring(RPCthreshold);
					TH1I* BProfile = (TH1I*)ROOTinput.Get(BProfName.c_str());

					//*****************************************************************************************************

					TF1* ForwardMuons = new TF1("fit_f","[0]*exp(-0.5*((x-[1])/[2])**2) / (1 + exp(-[3]*(x-[4])))",1,32);
					ForwardMuons->SetParLimits(0,1,10000);
					ForwardMuons->SetParameter(1,26);
					ForwardMuons->SetParLimits(1,1,32);
					ForwardMuons->SetParameter(2,6);
					ForwardMuons->SetParLimits(2,1,32);
					ForwardMuons->SetParameter(3,-1);
					ForwardMuons->SetParLimits(3,-10,-0.01);
					ForwardMuons->SetParameter(4,24);
					ForwardMuons->SetParLimits(4,1,32);
					FProfile->Fit(ForwardMuons,"QR");

					TVirtualFitter* fitterF = TVirtualFitter::GetFitter();
					assert(fitterF != 0);
					double * covMatrixF = fitterF->GetCovarianceMatrix();

					CSVFile << ts << '\t'
							<< FProfile->GetEntries() << '\t'
							<< ForwardMuons->Integral(1.,32.) << '\t'
							<< ForwardMuons->IntegralError(1.,32.) << '\t'
							<< FProfile->Integral(1,16) << '\t'
							<< ForwardMuons->Integral(1.,16.) << '\t'
							<< ForwardMuons->IntegralError(1.,16.) << '\t'
							<< FProfile->GetMaximum() << '\t'
							<< FProfile->GetMaximumBin() << '\t'
							<< ForwardMuons->GetMaximum(19.,23.) << '\t'
							<< ForwardMuons->GetMaximumX(19.,23.) << '\t';

					delete ForwardMuons;
					delete fitterF;

					//*****************************************************************************************************

					TF1* BackwardMuons = new TF1("fit_b","[0]*exp(-0.5*((x-[1])/[2])**2) / (1 + exp(-[3]*(x-[4])))",1,32);
					BackwardMuons->SetParLimits(0,1,10000);
					BackwardMuons->SetParameter(1,27);
					BackwardMuons->SetParLimits(1,1,32);
					BackwardMuons->SetParameter(2,2);
					BackwardMuons->SetParLimits(2,1,32);
					BackwardMuons->SetParameter(3,1);
					BackwardMuons->SetParLimits(3,0.01,10);
					BackwardMuons->SetParameter(4,24);
					BackwardMuons->SetParLimits(4,1,32);
					BProfile->Fit(BackwardMuons,"QR");

					TVirtualFitter* fitterB = TVirtualFitter::GetFitter();
					assert(fitterB != 0);
					double * covMatrixB = fitterB->GetCovarianceMatrix();

					CSVFile << BProfile->GetEntries() << '\t'
							<< BackwardMuons->Integral(1.,32.) << '\t'
							<< BackwardMuons->IntegralError(1.,32.) << '\t'
							<< BProfile->GetMaximum() << '\t'
							<< BProfile->GetMaximumBin() << '\t'
							<< BackwardMuons->GetMaximum(24.,27.) << '\t'
							<< BackwardMuons->GetMaximumX(24.,27.) << '\t';

					delete BackwardMuons;
					delete fitterB;

					//*****************************************************************************************************

					CSVFile << endl;
				}

				CSVFile.close();
			}
			ROOTinput.Close();

			return 0;
		}


}
