#include <iostream>
#include <fstream>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const bool debug = true;

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

  // Dimensions of the GEM
  const double pitch = 0.014;
  const double kapton = 50.e-4;
  const double metal = 5.e-4;
  const double outdia = 70.e-4;
  const double middia = 50.e-4;

  const bool plotField = true;
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
    fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
    fieldView->SetVoltageRange(-160., 160.);
    TCanvas* cF = new TCanvas();
    fieldView->SetCanvas(cF);
    fieldView->PlotContour();
  }

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 70., "co2", 30.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();
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
  sensor->SetArea(-5 * pitch, -5 * pitch, -0.03,
                   5 * pitch,  5 * pitch,  0.03);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) {
    driftView->SetArea(-2 * pitch, -2 * pitch, -0.02,
                        2 * pitch,  2 * pitch,  0.02);
    // Plot every 10 collisions (in microscopic tracking).
    aval->SetCollisionSteps(10); 
    aval->EnablePlotting(driftView);
    drift->EnablePlotting(driftView);
  }

  // Histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              nBinsGain, gmin, gmax);
  TH1F* hIons = new TH1F("hIons", "Number of ions",
                         nBinsGain, gmin, gmax);

  int nBinsChrg = 100;
  TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic",
                          nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
  TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic", 
                          nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);

  double sumIonsTotal = 0.;
  double sumIonsDrift = 0.;
  double sumIonsPlastic = 0.;

  double sumElectronsTotal = 0.;
  double sumElectronsPlastic = 0.;
  double sumElectronsUpperMetal = 0.;
  double sumElectronsLowerMetal = 0.;
  double sumElectronsTransfer = 0.;
  double sumElectronsOther = 0.;

  const int nEvents = 10;
  for (int i = nEvents; i--;) { 
    if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position.
    const double smear = pitch / 2.; 
    double x0 = -smear + RndmUniform() * smear;
    double y0 = -smear + RndmUniform() * smear;
    double z0 = 0.025; 
    double t0 = 0.;
    double e0 = 0.1;
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    hElectrons->Fill(ne);
    hIons->Fill(ni); 
    const int np = aval->GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (int j = np; j--;) {
      aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                   xe2, ye2, ze2, te2, e2, status);
      sumElectronsTotal += 1.;
      if (ze2 > -kapton / 2. && ze2 < kapton / 2.) {
        hChrgE->Fill(ze2 * 1.e4);
        sumElectronsPlastic += 1.;
      } else if (ze2 >= kapton / 2. && ze2 <= kapton  / 2. + metal) {
        sumElectronsUpperMetal += 1.;
      } else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
        sumElectronsLowerMetal += 1.;
      } else if (ze2 < -kapton / 2. - metal) {
        sumElectronsTransfer += 1.;
      } else {
        sumElectronsOther += 1.;
      }
      drift->DriftIon(xe1, ye1, ze1, te1);
      drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
                               xi2, yi2, zi2, ti2, status);
      if (zi1 < 0.01) {
        sumIonsTotal += 1.;
        if (zi2 > 0.01) sumIonsDrift += 1.;
      }
      if (zi2 > -kapton / 2. && zi2 < kapton / 2.) {
        hChrgI->Fill(zi2 * 1.e4);
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
            << sumElectronsPlastic / nEvents << "\n";
  std::cout << "Mean number of ions on plastic: "
            << sumIonsPlastic / nEvents << "\n";
 
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

  TCanvas* cD = new TCanvas();
  const bool plotGeo = true;
  if (plotGeo && plotDrift) {
    // Build the geometry in Root.
    TGeoManager* geoman = new TGeoManager("world", "geometry");
    TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    TGeoMaterial* matKapton = new TGeoMaterial("Kapton", 12, 6, 1.42);
    TGeoMedium* medKapton = new TGeoMedium("Kapton", 2, matKapton);
    TGeoMaterial* matCopper = new TGeoMaterial("Copper", 63, 29, 8.94);
    TGeoMedium* medCopper = new TGeoMedium("Copper", 3, matCopper);
    TGeoVolume* volTop = geoman->MakeBox("TOP", 
                                         medVacuum, pitch, pitch, 0.02);
    volTop->SetVisibility(0);
    TGeoBBox* shpKapton = new TGeoBBox("K", pitch / 2., 
                                            pitch / 2., 
                                            kapton / 2.);
    TGeoPcon* shpHole = new TGeoPcon("H", 0., 360., 3);
    shpHole->DefineSection(0, -kapton / 2., 0., outdia / 2.);
    shpHole->DefineSection(1,           0., 0., middia / 2.);
    shpHole->DefineSection(2,  kapton / 2., 0., outdia / 2.);

    TGeoCompositeShape* shpGem = new TGeoCompositeShape("G", "K - H");
    TGeoVolume* volKapton = new TGeoVolume("Kapton", shpGem, medKapton);
    volKapton->SetLineColor(kGreen);
    volKapton->SetTransparency(50);

    TGeoBBox* shpMetal = new TGeoBBox("M", pitch / 2., 
                                           pitch / 2., 
                                           metal / 2.);
    TGeoTube* shpTube = new TGeoTube("T", 0., outdia / 2., metal / 2.);
    TGeoCompositeShape* shpElectrode = new TGeoCompositeShape("E", "M - T");
    TGeoVolume* volElectrode = new TGeoVolume("Electrode", 
                                              shpElectrode, medCopper);
    volElectrode->SetLineColor(kBlue);
    volElectrode->SetTransparency(50);

    TGeoVolumeAssembly* volGem = new TGeoVolumeAssembly("Gem");
    const double shift =  0.5 * (metal + kapton);
    volGem->AddNode(volKapton, 1);
    volGem->AddNode(volElectrode, 2, new TGeoTranslation(0., 0.,  shift));
    volGem->AddNode(volElectrode, 3, new TGeoTranslation(0., 0., -shift));

    volTop->AddNode(volGem, 1);
    volTop->AddNode(volGem, 2, new TGeoTranslation(-pitch, 0., 0.));
    volTop->AddNode(volGem, 3, new TGeoTranslation(+pitch, 0., 0.));
    volTop->AddNode(volGem, 4, 
               new TGeoTranslation(-pitch / 2., sqrt(3) * pitch / 2., 0.));
    volTop->AddNode(volGem, 5, 
               new TGeoTranslation(+pitch / 2., sqrt(3) * pitch / 2., 0.));
    volTop->AddNode(volGem, 6,
               new TGeoTranslation(-pitch / 2., -sqrt(3) * pitch / 2., 0.));
    volTop->AddNode(volGem, 7,
               new TGeoTranslation(+pitch / 2., -sqrt(3) * pitch / 2., 0.));
    geoman->SetVerboseLevel(0);
    geoman->SetTopVolume(volTop);
    geoman->CloseGeometry();
    geoman->CheckOverlaps(0.1e-4);
    geoman->SetNmeshPoints(100000);
    cD->cd();
    geoman->GetTopVolume()->Draw("ogl");
  }

  if (plotDrift) {
    driftView->SetCanvas(cD);
    driftView->Plot(false, true);
  }

  const bool plotHistogram = true;
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

  app.Run(kTRUE);

}
