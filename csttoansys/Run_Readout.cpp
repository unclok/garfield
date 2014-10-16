#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <float.h>
#include "CST_Readout.h"
#include "CSTResultReaderInterf.h"

using namespace std;

int main() {
  string projectName;
  cout << "Give name and the path of the project to read (...Project.cst):\n"
      << endl;
  cin >> projectName;
  //** Here it is assumed that the potential (Potential [Es]) was copied and the new name is Potential
  string TreeFieldName("2D/3D Results\\Potential");
  CST_Readout object(projectName.c_str());

  int *meshnodes = object.Getnxyz();
  cout << "Mesh nodes x: " << meshnodes[0] << "\n Mesh nodes y: "
      << meshnodes[1] << "\n Mesh nodes z: " << meshnodes[2]
      << "\n Mesh nodes all: " << object.GetNumberOfNodes()
      << "\n Meshelements: " << object.GetNumberOfElements() << endl;
  //**Here all information needed by Garfield is exported except the potentials (due to a bug this needs to be done within CST using  Get_Field.bas)
  string target;
  cout << "Give the name of the directory where to write the resulting file:\n"
      << endl;
  cin >> target;
  object.WriteElementComp(target.c_str());
  object.WritePosLines(target.c_str());
  object.WritePotentials(TreeFieldName,target.c_str());
  cout << "Everything is done." << endl;
  return 1;
}
