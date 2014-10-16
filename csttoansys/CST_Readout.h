/**
 * Author: Klaus Zenker (DESY)
 * Version: 15.10.2013
 * - includes bugfix related to the handling of elements having PEC nodes
 * - includes potential export
 */

#pragma once

#ifndef CST_Readout_h
#define CST_Readout_h
#include "CSTResultReaderInterf.h"
#include <vector>

using namespace std;

#define eps0 8.85418782e-12
#define mue0 M_PI*4e-7

class CST_Readout {
private:
  CSTProjHandle m_pHandle;
  int uReturnVal;
  int m_fieldType;
  float* m_fieldData;
  float* m_meshMatrix;
  //mesh information
  int* m_nxyz; //array node numbers in x,y,z (nxyz[0...3])
  double* m_xyzlines;
  vector<double> m_xlines;
  vector<double> m_ylines;
  vector<double> m_zlines;
  void GetInv(int i, int j, int k, double &inv_x, double &inv_y, double &inv_z);

public:
  CST_Readout(string cst_file);
  virtual ~CST_Readout();
  static void PrintStatus(int event, int all_events, double percent = 0.1);
  static double Round(double number, int digits);
  void Node2Index(int node, int &i, int &j, int &k);
  void Element2Index(int element, int &i, int &j, int &k);
  void GetPositionFromNode(int node, double* &xyz);
  void GetPositionFromIndex(int i, int j, int k, double* &xyz);
  int GetNumberOfNodes() {
    return m_nxyz[0] * m_nxyz[1] * m_nxyz[2];
  }
  int GetNumberOfElements() {
    return (m_nxyz[0] - 1) * (m_nxyz[1] - 1) * (m_nxyz[2] - 1);
  }
  double GetMaterialPropertyFromIndex(int i, int j, int k, int MaterialType);
  double GetMaterialPropertyFromIndex(float* matMatrix_cst, int i, int j, int k,
      int MaterialType);
  void GetNodesForElement(int element, vector<int> &nodes);
  double GetElementMaterialProperty(int element, int MaterialType,
      int testnode = 3);
  int GetElementMaterial(int element, int node = 3);
  int* Getnxyz() {
    return m_nxyz;
  }
  double* Getxyzlines() {
    return m_xyzlines;
  }
  float* GetFieldData() {
    return m_fieldData;
  }
  void WriteElementComp(string path, string prefix = "");
  void WritePosLines(string path, string prefix = "");
  void WritePotentials(string TreeFieldName, string path, string prefix = "");
  int  RetrieveFieldData(string TreeFieldName);
  void AnalyseField(double scale = 1.);
  void AnalyseMaterial(double scale = 1.);
};

struct IntCmp {
  bool operator()(const pair<double, int> &lhs, const pair<double, int> &rhs) {
    return lhs.second > rhs.second;
  }
};

struct MaterialObject {
  vector<int*> position;
  double eps;
  int neps;
  bool MaterialObject::operator==(const MaterialObject &other) {
    return eps == other.eps;
  }
  friend ostream& operator<<(ostream &o, const MaterialObject & obj) {
    if (obj.neps != obj.position.size()) {
      o << obj.eps << "\t:" << obj.neps << "(but vec_size is: "
          << obj.position.size() << ")\t at (" << obj.position.at(0)[0] << ","
          << obj.position.at(0)[1] << "," << obj.position.at(0)[2] << ")"
          << endl;
    } else if (obj.position.size() > 1) {
      o << obj.eps << "\t:" << obj.neps << "\t at e.g. ("
          << obj.position.at(0)[0] << "," << obj.position.at(0)[1] << ","
          << obj.position.at(0)[2] << ")";
    } else {
      o << obj.eps << "\t:" << obj.neps << "\t at\t (" << obj.position.at(0)[0]
          << "," << obj.position.at(0)[1] << "," << obj.position.at(0)[2]
          << ")";
    }
    return o;
  }
};

struct Sorteps {
  bool operator()(const MaterialObject &left, const MaterialObject &right) {
    return left.eps < right.eps;
  }
};

struct Sortneps {
  bool operator()(const MaterialObject &left, const MaterialObject &right) {
    return left.neps > right.neps;
  }
};

#endif
