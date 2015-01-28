#include <iostream>

#include "TCanvas.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoCone.h"
#include "TGeoPcon.h"
#include "TGeoTube.h"
#include "TROOT.h"
#include "TApplication.h"
//#include "TGeoMatrix.h"
#include "TGeoCompositeShape.h"

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "GeometryRoot.hh"
#include "ViewGeometry.hh"
#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include <omp.h>

using namespace Garfield;
	//create global histograms to use userHandle function in the main function.
	TH1F* energyatend = new TH1F("energyatend","Energy at end",100,0,100);
	TH1F* pofe = new TH1F("pofe","position of electron below -55 um",200,-100,100);
void userHandle(double x, double y, double z, double t,
		double e, double dx, double dy, double dz,
							bool hole) {
  double um = 1e-04;
	if(z< -55*um){energyatend->Fill(e);
			pofe->Fill(z/um);}
	//else std::cout<<"this point : "<<x<<", "<<y<<", "<<z<<std::endl;
}


int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 // plottingEngine.SetDefaultStyle();
 
  const bool debug = true;
  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 293.15;
 
  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();
  const std::string efile = "ELIST.lis";
  const std::string nfile = "NLIST.lis";
  const std::string mfile = "MPLIST.lis";
  const std::string sfile = "PRNSOL.lis";
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();

  double um = 1e-04;
  // Dimensions of the GEM
  const double pitch = 0.014;
  const double length = 100;
  const double kapton_thic = 50;
  const double metal_thic = 5;
  const double outdia = 70;
  const double middia = 30;

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("Ar", 100);
  gas->Initialise();  
  // Set the Penning transfer efficiency.
  const double rPenning = 0.57;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas->LoadIonMobility("IonMobility_Ar+_Ar.txt");
  
  // Associate the gas with the corresponding field map material. 
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-length/2*um, -length/2*um, -60*um,
                   length/2*um,  length/2*um,  60*um);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  //aval->SetElectronTransportCut(0.8);

  //AvalancheMC* drift = new AvalancheMC();
  //drift->SetSensor(sensor);
  //drift->SetDistanceSteps(2.e-4);

  //draw a drift plot(attatched to the geometry)
  const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) {
    driftView->SetArea(-length/2*um, -length/2*um, -60*um,
                        length/2*um,  length/2*um,  60*um);
    // Plot every 10 collisions (in microscopic tracking).
    driftView->EnableDebugging();
    //driftView->SetClusterMarkerSize(0.1);
    //driftView->SetCollisionMarkerSize(10);
    aval->SetCollisionSteps(20); 
    aval->EnablePlotting(driftView);
    //drift->EnablePlotting(driftView);
  }

  //draw a field plot. The field can be a Voltage and an electric field.
  const bool plotField = false;
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    fieldView->SetPlane(0., 0., 1., 0., 0., 30*um);
    fieldView->SetArea(-length*um / 2., -50*um, length*um / 2., 50*um);
    //fieldView->SetVoltageRange(-200., 200.);
    fieldView->SetElectricFieldRange(-40000., 60000.);
    TCanvas* cF = new TCanvas();
    fieldView->SetCanvas(cF);
    fieldView->PlotContour("ez");
  }

  // create histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              nBinsGain, gmin, gmax);
  TH1F* hIons = new TH1F("hIons", "Number of ions",
                         nBinsGain, gmin, gmax);

  int nBinsChrg = 100;
  TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic",
                          nBinsChrg, -0.5 * kapton_thic*um, 0.5 * kapton_thic*um);
  TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic", 
                          nBinsChrg, -0.5 * kapton_thic*um, 0.5 * kapton_thic*um);
  TH1F* initE = new TH1F("initE", "Initial energy of electrons", 1000, 0, 100);
  TH1F* finalE = new TH1F("finalE", "Final energy of electrons", 1000, 0, 100);

  // make variables to sum electron and ion values.
  double sumIonsTotal = 0.;
  double sumIonsDrift = 0.;
  double sumIonsPlastic = 0.;

  double sumElectronsTotal = 0.;
  double sumElectronsPlastic = 0.;
  double sumElectronsUpperMetal = 0.;
  double sumElectronsLowerMetal = 0.;
  double sumElectronsTransfer = 0.;
  double sumElectronsOther = 0.;
  double sumElectronsGain = 0.;

  // input energy from argv[1], total number of events from argv[2], number tag for filename from argv[3]
  const int *nEvents = new int(1);
	
	int* temp = const_cast<int*>(nEvents);
	if(std::atoi(argv[2])>1)*temp = std::atoi(argv[2]);
	std::cout<<"*nEvents : "<<*nEvents<<", argv[2] : "<<argv[2]<<", temp : "<<*temp<<std::endl;
    double e0 = 10;
    if(argv[1])e0 = std::atoi(argv[1]);

    int num=0;
    if(argv[3])num = std::atoi(argv[3]);

  // use this canvas for geometry
  TCanvas* cD = new TCanvas();

  // geometry setting
  const bool plotGeo = true;
  if (plotGeo){// && plotDrift) {
  TGeoManager* mgr = new TGeoManager("GEM","GEM simulation 100*100*60");
  TGeoElementTable* table = gGeoManager->GetElementTable();
  TGeoElement* elH = table->FindElement("H");
  TGeoElement* elC = table->FindElement("C");
  TGeoElement* elN = table->FindElement("N");
  TGeoElement* elO = table->FindElement("O");
  TGeoElement* elCu = table->FindElement("Cu");
  TGeoElement* elAr = table->FindElement("Ar");

  TGeoMixture* elkap = new TGeoMixture("kapton_material", 4); 
  elkap->AddElement(elO,5);
  elkap->AddElement(elC,22);
  elkap->AddElement(elN,2);
  elkap->AddElement(elH,10);
  TGeoMaterial* kapmat = elkap;
  TGeoMaterial* cumat = new TGeoMaterial("copper_material", elCu, 8.96);
  TGeoMaterial* armat = new TGeoMaterial("argon_material", elAr, 1.782e-03);

  TGeoMedium* kapmed = new TGeoMedium("kapton_medium", 1, kapmat);
  TGeoMedium* cumed = new TGeoMedium("copper_medium", 2, cumat);
  TGeoMedium* armed = new TGeoMedium("argon_medium", 3, armat);

  TGeoVolume *world = mgr->MakeBox("world", armed, 50*um, 50*um, 60*um);
  mgr->SetTopVolume(world);

  TGeoBBox* copper = new TGeoBBox("copper", length/2*um, length/2*um, metal_thic/2*um);
  TGeoTube* copper_hole = new TGeoTube("copper_hole", 0.*um, outdia/2.*um, metal_thic/2*um);
  TGeoCompositeShape *copper_solid = new TGeoCompositeShape("copper_solid","copper-copper_hole");
  TGeoVolume *copper_volume = new TGeoVolume("copper_volume", copper_solid, cumed);
  copper_volume->SetLineColor(kBlue);
  copper_volume->SetTransparency(50);

  TGeoTranslation *t1 = new TGeoTranslation("t1",0*um,0*um,(kapton_thic+metal_thic)/2*um);
  t1->RegisterYourself();
  TGeoTranslation *t2 = new TGeoTranslation("t2",0*um,0*um,-(kapton_thic+metal_thic)/2*um);
  t2->RegisterYourself();

  TGeoBBox* kapton = new TGeoBBox("kapton", length/2*um, length/2*um, kapton_thic/2*um);
  TGeoPcon* argon_hole = new TGeoPcon("argon_hole", 0, 360, 3);
  argon_hole->DefineSection(0, -kapton_thic/2*um, 0, outdia/2.*um);
  argon_hole->DefineSection(1, 0*um, 0, middia/2.*um);
  argon_hole->DefineSection(2, kapton_thic/2*um, 0, outdia/2.*um);
  TGeoCompositeShape *kapton_solid = new TGeoCompositeShape("kapton_solid","kapton-argon_hole");
  TGeoVolume *kapton_volume = new TGeoVolume("kapton_volume",kapton_solid, kapmed);
  kapton_volume->SetLineColor(kGreen);
  kapton_volume->SetTransparency(50);

  TGeoVolumeAssembly *gem_volume = new TGeoVolumeAssembly("gem");
  gem_volume->AddNode(kapton_volume,1);
  gem_volume->AddNode(copper_volume,2,t2);
  gem_volume->AddNode(copper_volume,3,t1);

  // set geometry to garfield
  GeometryRoot* geo = new GeometryRoot();
  world->AddNode(gem_volume,1);
  geo->SetGeometry(mgr);
  mgr->CloseGeometry();
  mgr->CheckOverlaps(0.1e-4);
  mgr->SetNmeshPoints(10000);
  cD->cd();
  //mgr->GetTopVolume()->Raytrace();
  mgr->GetTopVolume()->Draw("ogl");
  }

  for (int i = *nEvents; i>0; i-- ) { 
    if (debug && i % 10 == 0) std::cout << i << "/" << *nEvents << "\n";
//    AvalancheMicroscopic avaltemp = *aval;
//    AvalancheMC drifttemp = *drift;
    // Randomize the initial position.
    const double smear = length*um/2; 
    //double x0 = -smear + RndmUniform() * smear;
    //double y0 = -smear + RndmUniform() * smear;
    double x0 = 0;
    double y0 = 0;
    double z0 = 55*um; 
    double t0 = 0.;
	aval->SetUserHandleStep(userHandle);
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 1.);
    //avaltemp.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    //#pragma omp critical
    aval->GetAvalancheSize(ne, ni);
    //avaltemp.GetAvalancheSize(ne, ni);
    //#pragma omp critical
    hElectrons->Fill(ne);
    //#pragma omp critical
    hIons->Fill(ni); 
    const int np = aval->GetNumberOfElectronEndpoints();
    //const int np = avaltemp.GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (int j = np;j>0; j--) {
    //#pragma omp critical
      aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                   xe2, ye2, ze2, te2, e2, status);
	//std::cout<<"end point energy : " << e2 << std::endl;
	initE->Fill(e1);
	finalE->Fill(e2);
      //avaltemp.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
      //                             xe2, ye2, ze2, te2, e2, status);
      sumElectronsTotal += 1.;
      if (ze2 > -kapton_thic*um / 2. && ze2 < kapton_thic*um / 2.) {
    //#pragma omp critical
        hChrgE->Fill(ze2);
        sumElectronsPlastic += 1.;
      } else if (ze2 >= kapton_thic*um / 2. && ze2 <= kapton_thic*um  / 2. + metal_thic*um) {
        sumElectronsUpperMetal += 1.;
      } else if (ze2 <= -kapton_thic*um / 2. && ze2 >= -kapton_thic*um / 2. - metal_thic*um) {
        sumElectronsLowerMetal += 1.;
      } else if (ze2 < -55*um) {
        sumElectronsGain += 1.;
      } else if (ze2 < -kapton_thic*um / 2. - metal_thic*um) {
        sumElectronsTransfer += 1.;
      } else {
        sumElectronsOther += 1.;
      }
      //drift->DriftIon(xe1, ye1, ze1, te1);
      //drifttemp.DriftIon(xe1, ye1, ze1, te1);
    //#pragma omp critical
      //drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
      //                         xi2, yi2, zi2, ti2, status);
      //drifttemp.GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
      //                         xi2, yi2, zi2, ti2, status);
      if (zi1 < 0.01) {
        sumIonsTotal += 1.;
        if (zi2 > 0.01) sumIonsDrift += 1.;
      }
      if (zi2 > -kapton_thic / 2. && zi2 < kapton_thic / 2.) {
    //#pragma omp critical
        hChrgI->Fill(zi2);
        sumIonsPlastic += 1.;
      }
    }
  }
  double fFeedback = 0.;
  if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
  std::cout << "Fraction of ions drifting back: " << fFeedback << "\n"; 

  const double neMean = hElectrons->GetMean();
  std::cout << "Mean number of electrons: " << neMean << "\n";
  const double niMean = hIons->GetMean();
  std::cout << "Mean number of ions: " << niMean << "\n";

  std::cout << "Mean number of electrons on plastic: "
            << sumElectronsPlastic / *nEvents << "\n";
  std::cout << "Mean number of ions on plastic: "
            << sumIonsPlastic / *nEvents << "\n";
 
  std::cout << "Electron endpoints:\n";
  const double fUpperMetal = sumElectronsUpperMetal / sumElectronsTotal;
  const double fPlastic = sumElectronsPlastic / sumElectronsTotal;
  const double fLowerMetal = sumElectronsLowerMetal / sumElectronsTotal;
  const double fTransfer = sumElectronsTransfer / sumElectronsTotal;
  const double fOther = sumElectronsOther / sumElectronsTotal;
  std::cout << "    upper metal: " << fUpperMetal * 100. << "%\n";
  std::cout << "    plastic:     " << fPlastic * 100. << "%\n";
  std::cout << "    lower metal: " << fLowerMetal * 100. << "%\n";
  std::cout << "    transfer:    " << fTransfer * 100. << "%\n";
  std::cout << "    other:       " << fOther * 100. << "%\n";
  std::cout << "    gain:        " << sumElectronsTransfer / *nEvents << std::endl;
  std::cout << "    gain:        " << sumElectronsGain / *nEvents << std::endl;

  // draw the geometry with drifts
  if (plotDrift) {
    driftView->SetCanvas(cD);
    driftView->Plot(false, false);
  char hname[3];
  sprintf(hname, "geo%f_%d.root",e0,num);
  cD->SaveAs(hname);
  }

  // draw histogram with summation values
  const bool plotHistogram = false;
  if (plotHistogram) {
    TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
    cH->Divide(2, 2);
    cH->cd(1);
    hElectrons->Draw();
    cH->cd(2);
    hIons->Draw();
    cH->cd(3);
    hChrgE->Draw();
    cH->cd(4);
    hChrgI->Draw();
  }


  // view geometry
//  ViewGeometry* view = new ViewGeometry();

  // Set the field range to be covered by the gas table. 
//  const int nFields = 20;
//  const double emin =    100.;
//  const double emax = 100000.;
  // Flag to request logarithmic spacing.
//  const bool useLog = false;
 // gas->SetFieldGrid(emin, emax, nFields, useLog); 

//  const int ncoll = 10;
  // Switch on debugging to print the Magboltz output.
//  gas->EnableDebugging();
  // Run Magboltz to generate the gas table.
// gas->GenerateGasTable(ncoll);
  //gas->DisableDebugging();
  // Save the table. 
//  gas->WriteGasFile("ar_93_co2_7.gas");
 // file save as root
 char hname[3];
 sprintf(hname, "energy%f_%d.root",e0,num);
 energyatend->SaveAs(hname);
 sprintf(hname, "pofe%f_%d.root",e0,num);
 pofe->SaveAs(hname);
 sprintf(hname, "initenergy%f_%d.root",e0,num);
 initE->SaveAs(hname);
 sprintf(hname, "finalenergy%f_%d.root",e0,num);
 finalE->SaveAs(hname);

  app.Run(kTRUE);
//	std::exit(1);

}
