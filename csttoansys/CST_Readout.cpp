/**
 * Author: Klaus Zenker (DESY)
 * Version: 15.10.2013
 * - includes bugfix related to the handling of elements having PEC nodes
 * - includes potential export
 */

#define _USE_MATH_DEFINES
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

#include "CST_Readout.h"
#include "CSTResultReaderInterf.h"

void CST_Readout::GetInv(int i, int j, int k, double &inv_x, double &inv_y,
    double &inv_z) {
  /*Here the scaling factor for the material calculation is computed*/
  // Cell lengths in x/y/z direction at index i/j/k and (i+1)/(j+1)/(k+1).
  double pLxi1 = 0;
  //for the last node choose also half of the element size (like for the first cell in the dual grid
  if (i == (m_xlines.size() - 1)) {
    pLxi1 = m_xlines.at(i) - m_xlines.at(i - 1);
  } else {
    pLxi1 = m_xlines.at(i + 1) - m_xlines.at(i);
  }
  double pLyj1 = 0;
  if (j == (m_ylines.size() - 1)) {
    pLyj1 = m_ylines.at(j) - m_ylines.at(j - 1);
  } else {
    pLyj1 = m_ylines.at(j + 1) - m_ylines.at(j);
  }
  double pLzk1 = 0;
  if (k == (m_zlines.size() - 1)) {
    pLzk1 = m_zlines.at(k) - m_zlines.at(k - 1);
  } else {
    pLzk1 = m_zlines.at(k + 1) - m_zlines.at(k);
  }
  // Cell face areas with face normal x/y/z at index i/j/k -> not used at all!
  double pAxjk = pLyj1 * pLzk1;
  double pAyik = pLxi1 * pLzk1;
  double pAzij = pLxi1 * pLyj1;
  /* Now dual (d) grid: Lengths L, Areas A

   Cell lengths in x/y/z direction at index i/j/k.
   Dual grid is a bit more complex since first and last cell are different
   from the rest.*/
  double dLxi, dLyj, dLzk;
  if (i == 0 || i == (m_xlines.size() - 1)) {
    dLxi = pLxi1 / 2; // first dual cell length is half the length of first primary cell
  } else {
    double pLxi2 = m_xlines.at(i) - m_xlines.at(i - 1);
    dLxi = pLxi1 / 2 + pLxi2 / 2;
  }

  if (j == 0 || j == (m_ylines.size() - 1)) {
    dLyj = pLyj1 / 2; // first dual cell length is half the length of first primary cell
  } else {
    double pLyj2 = m_ylines.at(j) - m_ylines.at(j - 1);
    dLyj = pLyj1 / 2 + pLyj2 / 2;
  }

  if (k == 0 || k == (m_zlines.size() - 1)) {
    dLzk = pLzk1 / 2; // first dual cell length is half the length of first primary cell
  } else {
    double pLzk2 = m_zlines.at(k) - m_zlines.at(k - 1);
    dLzk = pLzk1 / 2 + pLzk2 / 2;
  }
  /*Cell face areas with face normal x/y/z at index i/j/k
   One more cell per direction to cover the same total area as the primary grid*/
  double dAxjk = dLyj * dLzk;
  double dAyik = dLxi * dLzk;
  double dAzij = dLxi * dLyj;

  // x components
  inv_x = dAxjk / pLxi1;
  // y components
  inv_y = dAyik / pLyj1;
  // z components
  inv_z = dAzij / pLzk1;
}

CST_Readout::CST_Readout(string cst_file) {
  m_meshMatrix = 0;
  m_fieldType = 0;
  m_fieldData = 0;
  m_pHandle.m_pProj = NULL;
  uReturnVal = CST_OpenProject(cst_file.c_str(), &m_pHandle);
  assert(!uReturnVal);
  m_nxyz = new int[3];
  uReturnVal = CST_GetHexMeshInfo(&m_pHandle, m_nxyz);
  assert(!uReturnVal);
  m_xyzlines = new double[m_nxyz[0] + m_nxyz[1] + m_nxyz[2]];
  uReturnVal = CST_GetHexMesh(&m_pHandle, m_xyzlines);
  assert(!uReturnVal);
  m_xlines.resize(m_nxyz[0]);
  m_ylines.resize(m_nxyz[1]);
  m_zlines.resize(m_nxyz[2]);
  copy(m_xyzlines, m_xyzlines + m_nxyz[0], m_xlines.begin());
  copy(m_xyzlines + m_nxyz[0], m_xyzlines + m_nxyz[0] + m_nxyz[1],
      m_ylines.begin());
  copy(m_xyzlines + m_nxyz[0] + m_nxyz[1],
      m_xyzlines + m_nxyz[0] + m_nxyz[1] + m_nxyz[2], m_zlines.begin());
  cout << "xlines: " << m_xlines.size() << "\t ylines: " << m_ylines.size()
      << "\t zlines: " << m_zlines.size() << endl;
}

CST_Readout::~CST_Readout() {
  uReturnVal = CST_CloseProject(&m_pHandle);
  assert(!uReturnVal);
  cout << "Prejct was closed without problems" << endl;
  delete[] m_fieldData;
  delete[] m_meshMatrix;
  delete[] m_nxyz;
  delete[] m_xyzlines;
}

void CST_Readout::PrintStatus(int event, int all_events, double percent) {
  int tmp = (int) (percent * all_events);
  if (event % tmp == 0) {
    cout << "--> " << (int) (event * 100. / all_events + .5)
        << "% \t is done so far." << endl;
  }
}

double CST_Readout::Round(double number, int digits) {
  double v[] = { 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 }; // to be enlarged for more digits
  return floor(number * v[digits] + 0.5) / v[digits];
}

void CST_Readout::Node2Index(int node, int &i, int &j, int &k) {
  j = 0;
  k = 0;
  int tmp = node;
  while (1) {
    k++;
    tmp -= m_nxyz[0] * m_nxyz[1];
    if (tmp < 0)
      break;
  }
  k--;
  tmp = node - k * m_nxyz[0] * m_nxyz[1];
  while (1) {
    j++;
    tmp -= m_nxyz[0];
    if (tmp < 0)
      break;
  }
  j--;
  i = node - j * m_nxyz[0] - k * m_nxyz[0] * m_nxyz[1];
}

void CST_Readout::Element2Index(int element, int &i, int &j, int &k) {
  j = 0;
  k = 0;
  int tmp = element;
  while (1) {
    k++;
    tmp -= (m_nxyz[0] - 1) * (m_nxyz[1] - 1);
    if (tmp < 0)
      break;
  }
  k--;
  tmp = element - k * (m_nxyz[0] - 1) * (m_nxyz[1] - 1);
  while (1) {
    j++;
    tmp -= (m_nxyz[0] - 1);
    if (tmp < 0)
      break;
  }
  j--;
  i = element - j * (m_nxyz[0] - 1) - k * (m_nxyz[0] - 1) * (m_nxyz[1] - 1);
}

void CST_Readout::GetPositionFromNode(int node, double* &xyz) {
  int i, j, k;
  CST_Readout::Node2Index(node, i, j, k);
  CST_Readout::GetPositionFromIndex(i, j, k, xyz);
}

void CST_Readout::GetPositionFromIndex(int i, int j, int k, double* &xyz) {
  xyz[0] = m_xlines.at(i);
  xyz[1] = m_ylines.at(j);
  xyz[2] = m_zlines.at(k);
}

double CST_Readout::GetMaterialPropertyFromIndex(int i, int j, int k,
    int MaterialType) {
  /*See CST_Readout::GetMaterialFromNode(float* matMatrix_cst,int i, int j, int k, int MaterialType) for details.     Here just the matMatrix_cst is retrieved before.*/
  if (m_meshMatrix == 0) {
    m_meshMatrix = new float[CST_Readout::GetNumberOfNodes() * 3];
    uReturnVal = CST_GetMaterialMatrixHexMesh(&m_pHandle, MaterialType,
        m_meshMatrix);
  }

  if (uReturnVal) {
    cout << "Didn't found quantity " << MaterialType << " due to error "
        << uReturnVal << " and skip it." << endl;
    assert(!uReturnVal);
    return 0;
  }
  return CST_Readout::GetMaterialPropertyFromIndex(m_meshMatrix, i, j, k,
      MaterialType);;
}

double CST_Readout::GetMaterialPropertyFromIndex(float* matMatrix_cst, int i,
    int j, int k, int MaterialType) {
  /*Here for a given position (node with indeces i,j,k) eps (matType=0), mu (matType=1), kappa (matType=2)
   or rho (matType	=3) are calculated. In the CST macro there is no problem with the last node, because the last node will never be used (for a given (x,y,z) the node before is choosen).*/
  // Define super-index and its maximum value
  int n = 1 + i + j * m_nxyz[0] + k * m_nxyz[0] * m_nxyz[1];
  if (n < 0 || n > CST_Readout::GetNumberOfNodes()) {
    cout
        << "In CST_Readout::GetMaterialFromIndex the superindex is greater #nodes or negative!"
        << endl;
    cout << "(i = " << i << ",j = " << j << ",k =  " << k << ")" << endl;
    return -1;
  }
  double inv_x, inv_y, inv_z;
  CST_Readout::GetInv(i, j, k, inv_x, inv_y, inv_z);

  double matValue[3];
  switch (MaterialType) {
  case 0: //eps
    matValue[0] = 1. / (matMatrix_cst[n - 1] * inv_x * eps0);
    matValue[1] =
        1.
            / (matMatrix_cst[CST_Readout::GetNumberOfNodes() + n - 1] * inv_y
                * eps0);
    matValue[2] = 1.
        / (matMatrix_cst[2 * CST_Readout::GetNumberOfNodes() + n - 1] * inv_z
            * eps0);
    break;
  case 1: //mue matrix
    matValue[0] = 1. / (matMatrix_cst[n - 1] * inv_x * mue0);
    matValue[1] =
        1.
            / (matMatrix_cst[CST_Readout::GetNumberOfNodes() + n - 1] * inv_y
                * mue0);
    matValue[2] = 1.
        / (matMatrix_cst[2 * CST_Readout::GetNumberOfNodes() + n - 1] * inv_z
            * mue0);
    break;
  case 2: // kap
    matValue[0] = matMatrix_cst[n - 1] / inv_x;
    matValue[1] = matMatrix_cst[CST_Readout::GetNumberOfNodes() + n - 1]
        / inv_y;
    matValue[2] = matMatrix_cst[2 * CST_Readout::GetNumberOfNodes() + n - 1]
        / inv_z;
    break;
  case 3: // rho
    matValue[0] = matMatrix_cst[n - 1];
    matValue[1] = matMatrix_cst[CST_Readout::GetNumberOfNodes() + n - 1];
    matValue[2] = matMatrix_cst[2 * CST_Readout::GetNumberOfNodes() + n - 1];
    break;
  default: // unsupported material MaterialType
    cout << "unsupported maeterial!" << endl;
    matValue[0] = -99999;
    matValue[1] = -99999;
    matValue[2] = -99999;
    break;
  }
  if (CST_Readout::Round(matValue[0], 5) != CST_Readout::Round(matValue[1], 5)
      || CST_Readout::Round(matValue[0], 5)
          != CST_Readout::Round(matValue[2], 5)) {
    if (matValue[0] <= DBL_MAX && matValue[0] >= -DBL_MAX
        && matValue[1] <= DBL_MAX && matValue[1] >= -DBL_MAX
        && matValue[2] <= DBL_MAX && matValue[2] >= -DBL_MAX) {
    }
  }
  return matValue[0];

}

void CST_Readout::GetNodesForElement(int element, vector<int> &nodes) {
  /*
   global coordinates   8 _ _ _ _7    t3    t2
                       /       /|     ^   /|
     ^ z              /       / |     |   /
     |             5 /_______/6 |     |  /
     |              |        |  |     | /
     |              |  4     |  /3    |/     t1
     ------->       |        | /       ------->
    /      y        |        |/       local coordinates
   /                1--------2
  /
 v x
 */

  int i, j, k;
  CST_Readout::Element2Index(element, i, j, k);
  nodes.clear();
  nodes.push_back((i + 1) + j * m_nxyz[0] + k * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back((i + 1) + (j + 1) * m_nxyz[0] + k * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back(i + (j + 1) * m_nxyz[0] + k * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back(i + j * m_nxyz[0] + k * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back((i + 1) + j * m_nxyz[0] + (k + 1) * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back(
      (i + 1) + (j + 1) * m_nxyz[0] + (k + 1) * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back(i + (j + 1) * m_nxyz[0] + (k + 1) * m_nxyz[0] * m_nxyz[1]);
  nodes.push_back(i + j * m_nxyz[0] + (k + 1) * m_nxyz[0] * m_nxyz[1]);
}

double CST_Readout::GetElementMaterialProperty(int element, int MaterialType,
    int testnode) {
  /*Here the material is given for the element like it is done GetPhzsicalMatMatrix.m.
   This  means the material of the element is given by node number 1 of the element!
   Thus here only the material for this node is read out.*/

  vector<int> nodes;
  GetNodesForElement(element, nodes);
  int i_2, j_2, k_2;
  CST_Readout::Node2Index(nodes.at(testnode), i_2, j_2, k_2);
  return CST_Readout::GetMaterialPropertyFromIndex(i_2, j_2, k_2, MaterialType);
}

int CST_Readout::GetElementMaterial(int element, int node) {
  double eps = CST_Readout::Round(
      CST_Readout::GetElementMaterialProperty(element, 0, node), 2);
  int mat;
  if (eps < 3.6 && eps > 3.4) {
    // kapton: epsilon = 3.5
    mat = 3;
  } else if (eps < 4.9 && eps > 4.7) {
    // G10: epsilon = 4.8
    mat = 4;
  } else if (eps < 9.5 && eps > 9.3) {
    // Alu: epsilon = 9.4
    mat = 5;
  } else if (eps < 3.4 && eps > 3.2) {
    // PEEK: epsilon = 3.3
    mat = 6;
  } else if (!(eps <= DBL_MAX && eps >= -DBL_MAX)) {
    mat = 2;
  } else {
    mat = 1;
  }
  return mat;
}

int CST_Readout::RetrieveFieldData(string TreeFieldName) {
  /*Here the field specified by Tree3DField is read.*/
  int n3DDataSize = 0;
  uReturnVal = CST_Get3DHexResultSize(&m_pHandle, TreeFieldName.c_str(), 0,
      &n3DDataSize);
  assert(!uReturnVal);
  int charBufSize = 0;
  int intInfo = 0;

  //****not necessary - just info****
  /*get the information of field in  Tree3DField (info is given for intSize fields,e.g. if there are
   3 fields in Tree3DField than intSize should be set to 3*/
  int infoSize = 1;
  uReturnVal = CST_Get3DHexResultInfo(&m_pHandle, TreeFieldName.c_str(), 0,
      infoSize, charBufSize, NULL, &intInfo, NULL);
  assert(!uReturnVal);
  m_fieldType = intInfo;
  cout << "Field Type: " << intInfo << endl;
  switch (intInfo) {
  case 1:
    cout << "Number of data found (Re,Im for E_x,E_y,E_z -> 6x#nodes): "
        << n3DDataSize << endl;
    break;
  case 2:
    cout << "This is a real vector field (Re for E_x,E_y,E_z -> 3x#nodes): "
        << n3DDataSize << endl;
    break;
  case 3:
    cout << "This is a complex scalar field (RE,Im -> 2x#nodes): "
        << n3DDataSize << endl;
    break;
  case 4:
    cout << "Number of data found (real scalar field -> #nodes): "
        << n3DDataSize << endl;
    break;
  }
  //****not necessary - just info****

  m_fieldData = new float[n3DDataSize * 6];
  // get the field vector
  uReturnVal = CST_Get3DHexResult(&m_pHandle, TreeFieldName.c_str(), 0,
      m_fieldData);
  if (uReturnVal != 0)
    cout << "Field readout error code: " << uReturnVal << endl;
  assert(!uReturnVal);
  return n3DDataSize;
}

void CST_Readout::WriteElementComp(string path, string prefix) {
  /*Here the file ELIST.lis is created in path. It contains the element number, the type of material (
   mat == 3: kapton [eps +- 0.4 = 3. ~ 4.]
   (check CST_Readout::AnalyseMaterial to see how much nodes are not equal 3.5),
   mat == 5: Al
   mat == 4: G10
   mat == 3: kapton
   mat == 2: PEC    [eps = #inf],
   mat == 1: vacuum [else] (should be mostly 1.),
   and the nodes corresponding to the element (element,node1,node2...).*/
  cout << "--->Working on CST_Readout::WriteElementComp" << endl;
  ofstream file;
  path += "\\";
  path += prefix;
  path += "\\ELIST1.lis";
  file.open(path.c_str());
  int element = 0;
  int mat_tmp;
  int vacuum = 0;
  int pec = 0;
  int other = 0;
  vector<int> nodes;
  for (int z = 0; z < m_nxyz[2] - 1; z++) {
    for (int y = 0; y < m_nxyz[1] - 1; y++) {
      for (int x = 0; x < m_nxyz[0] - 1; x++) {
        CST_Readout::PrintStatus(element, CST_Readout::GetNumberOfElements());
        //get eps and decide what is the coresponding material
        mat_tmp = CST_Readout::GetElementMaterial(element);
        // check also the 7th node to see if the whole element is made of PEC
        // for other material it seem to work to only consider the 1st node but not for PEC
        if (mat_tmp == 2)
          mat_tmp = CST_Readout::GetElementMaterial(element, 5);

        if (mat_tmp == 1) {
          vacuum++;
        } else if (mat_tmp == 2) {
          pec++;
        } else {
          other++;
        }
        file << " " << element << "    " << mat_tmp << "    \n";
        element++;
      }
    }
  }
  file.close();
  cout << "Wrote file " << path.c_str() << "\n and found \n" << vacuum
      << " vacuum elements, \n" << pec << " pec elements and \n" << other
      << " other elements." << endl;
  cout << "--->Leaving CST_Readout::WriteElementComp" << endl;
}

void CST_Readout::WritePosLines(string path, string prefix) {
  /*Here the file NODELINES.lis is created in path. It crontains the positions
   of the nodes for all possible positions (but not for each node)*/
  cout << "--->Working on CST_Readout::WritePosLines" << endl;
  ofstream file;
  path += "\\";
  path += prefix;
  path += "\\NODELINES.lis";
  file.open(path.c_str());
  file << " xmax " << m_xlines.size() << " ymax " << m_ylines.size() << " zmax "
      << m_zlines.size() << "\n";
  file << " x-lines \n";
  vector<double>::iterator it_lines = m_xlines.begin();
  int nall_lines = m_xlines.size() + m_ylines.size() + m_zlines.size();
  int line = 0;
  double scale_x, scale_y, scale_z;
  while (it_lines != m_xlines.end()) {
    CST_Readout::PrintStatus(line, nall_lines);
    CST_Readout::GetInv(line, 0, 0, scale_x, scale_y, scale_z);
    file << " " << (*it_lines) << "\n";
    line++;
    it_lines++;
  }
  file << " y-lines \n";
  it_lines = m_ylines.begin();
  while (it_lines != m_ylines.end()) {
    CST_Readout::PrintStatus(line, nall_lines);
    file << " " << (*it_lines) << "\n";
    line++;
    it_lines++;
  }
  file << " z-lines \n";
  it_lines = m_zlines.begin();
  while (it_lines != m_zlines.end()) {
    CST_Readout::PrintStatus(line, nall_lines);
    file << " " << *it_lines << "\n";
    line++;
    it_lines++;
  }
  file.close();
  cout << "Wrote file " << path.c_str() << endl;
  cout << "--->Leaving CST_Readout::WritePosLines" << endl;
}


void
CST_Readout::WritePotentials(string TreeFieldName, string path, string prefix){
/*Here the file PRNSOL.lis is created in path. It crontains the element
  number and the corresponding potential*/
  cout << "--->Working on CST_Readout::WritePotentials" << endl;
  ofstream file;
  path += "\\";
  path += prefix;
  path += "\\PRNSOL.lis";
  file.open(path.c_str());
  int n3dData = RetrieveFieldData(TreeFieldName);
  for(int data = 0; data < n3dData; data++){
    CST_Readout::PrintStatus(data,n3dData);
    file << data+1 << " " << setprecision(15) << m_fieldData[data] << "\n";
  }
  file.close();
  cout << "Wrote file " << path.c_str() << endl;
  cout << "Leaving CST_Readout::WritePotentials" << endl;
}

void CST_Readout::AnalyseField(double scale) {
  /*Here the field retrieved by CST_Readout::RetrieveFieldData are analysed (max,min,...). Scale is used to   get coordinates in project units instead of m (default after read out).*/
  cout << "--->Working on CST_Readout::AnalyseField" << endl;
  if (m_fieldData == 0) {
    cout
        << "You have to call CST_Readout::RetrieveFieldData before you can analyse the field!"
        << endl;
    return;
  }
  vector<float> RE_x, RE_y, RE_z;
  vector<float> IM_x, IM_y, IM_z;
  vector<float>::iterator it_RE_x;
  vector<float>::iterator it_RE_y;
  vector<float>::iterator it_RE_z;
  vector<float>::iterator it_IM_x;
  vector<float>::iterator it_IM_y;
  vector<float>::iterator it_IM_z;
  vector<float> Abs_field;
  vector<float>::iterator maximum;
  double* tmp_xyz;
  int nodes = CST_Readout::GetNumberOfNodes();
  float tmp_max = 0;
  switch (m_fieldType) {
  case 1:
    for (int i = 0; i < 2 * nodes; i++) {
      if (i % 2 == 0) {
        RE_x.push_back(m_fieldData[i]);
        RE_y.push_back(m_fieldData[i + 2 * nodes]);
        RE_z.push_back(m_fieldData[i + 2 * 2 * nodes]);
      } else {
        IM_x.push_back(m_fieldData[i]);
        IM_y.push_back(m_fieldData[i + 2 * nodes]);
        IM_z.push_back(m_fieldData[i + 2 * 2 * nodes]);
      }
    }
    it_RE_x = RE_x.begin();
    it_RE_y = RE_x.begin();
    it_RE_z = RE_x.begin();
    it_IM_x = IM_x.begin();
    it_IM_y = IM_y.begin();
    it_IM_z = IM_z.begin();
    while (it_RE_x != RE_x.end()) {
      Abs_field.push_back(
          sqrt(
              (*it_RE_x) * (*it_RE_x) + (*it_RE_y) * (*it_RE_y)
                  + (*it_RE_z) * (*it_RE_z) + (*it_IM_x) * (*it_IM_x)
                  + (*it_IM_y) * (*it_IM_y) + (*it_IM_z) * (*it_IM_z)));
      it_RE_x++;
      it_RE_y++;
      it_RE_z++;
      it_IM_x++;
      it_IM_y++;
      it_IM_z++;
    }
    maximum = max_element(RE_x.begin(), RE_x.end());
    tmp_xyz = new double[3];
    CST_Readout::GetPositionFromNode(distance(RE_x.begin(), maximum), tmp_xyz);
    cout << "Max. Field RE(x):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(RE_x.begin(), maximum)
        << " of " << RE_x.size() << ")\n" << endl;
    maximum = max_element(RE_y.begin(), RE_y.end());
    CST_Readout::GetPositionFromNode(distance(RE_y.begin(), maximum), tmp_xyz);
    cout << "Max. Field RE(y):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(RE_y.begin(), maximum)
        << " of " << RE_y.size() << ")\n" << endl;
    maximum = max_element(RE_z.begin(), RE_z.end());
    CST_Readout::GetPositionFromNode(distance(RE_z.begin(), maximum), tmp_xyz);
    cout << "Max. Field RE(z):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << " (node " << distance(RE_z.begin(), maximum)
        << " of " << RE_z.size() << ")\n" << endl;
    maximum = max_element(IM_x.begin(), IM_x.end());
    CST_Readout::GetPositionFromNode(distance(IM_x.begin(), maximum), tmp_xyz);
    cout << "Max. Field IM(x):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(IM_x.begin(), maximum)
        << " of " << IM_x.size() << ")\n" << endl;
    maximum = max_element(IM_y.begin(), IM_y.end());
    CST_Readout::GetPositionFromNode(distance(IM_y.begin(), maximum), tmp_xyz);
    cout << "Max. Field IM(y):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(IM_y.begin(), maximum)
        << " of " << IM_y.size() << ")\n" << endl;
    maximum = max_element(IM_z.begin(), IM_z.end());
    CST_Readout::GetPositionFromNode(distance(IM_z.begin(), maximum), tmp_xyz);
    cout << "Max. Field IM(z):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(IM_z.begin(), maximum)
        << " of " << IM_z.size() << ")\n" << endl;
    maximum = max_element(Abs_field.begin(), Abs_field.end());
    CST_Readout::GetPositionFromNode(distance(Abs_field.begin(), maximum),
        tmp_xyz);
    cout << "Max. Field Absolute:\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node "
        << distance(Abs_field.begin(), maximum) << " of " << Abs_field.size()
        << ")\n" << endl;
    delete tmp_xyz;
    break;

  case 4:
    for (int i = 0; i < nodes; i++) {
      RE_x.push_back(m_fieldData[i]);
    }
    maximum = max_element(RE_x.begin(), RE_x.end());
    tmp_xyz = new double[3];
    CST_Readout::GetPositionFromNode(distance(RE_x.begin(), maximum), tmp_xyz);
    cout << "Max. Field RE(x):\t" << (*maximum) << "\n at position: "
        << tmp_xyz[0] * scale << "\t" << tmp_xyz[1] * scale << "\t"
        << tmp_xyz[2] * scale << "\t(node " << distance(RE_x.begin(), maximum)
        << " of " << RE_x.size() << ")\n" << endl;
    delete tmp_xyz;
    break;

  default:
    cout << "The information for this kind of field (" << m_fieldType
        << ") is not yet implemented." << endl;
    break;
  }
  cout << "--->Leaving CST_Readout::AnalyseField" << endl;

}

void CST_Readout::AnalyseMaterial(double scale) {
  /*Here eps is read out (rounded to 2 digits) and then all eps are sorted by number
   of nodes corresponding to this eps. Afterwards eps for all elements (eps(element) = eps(first node))
   is calculated and it can be compared to the eps of the other nodes of the element. In the end
   the number of all kins of eps is given.*/
  cout << "--->Working on CST_Readout::AnalyseMaterial" << endl;
  bool test = true;
  int *indices;
  double eps;
  int is_event = 0;
  int nnodes = CST_Readout::GetNumberOfNodes();
  int nelements = CST_Readout::GetNumberOfElements();
  vector<MaterialObject> mat_vec;
  for (int i_x = 0; i_x < m_nxyz[0]; i_x++) {
    for (int i_y = 0; i_y < m_nxyz[1]; i_y++) {
      for (int i_z = 0; i_z < m_nxyz[2]; i_z++) {
        CST_Readout::PrintStatus(is_event, nnodes + nelements);
        eps = CST_Readout::Round(
            CST_Readout::GetMaterialPropertyFromIndex(i_x, i_y, i_z, 0), 2);
        MaterialObject tmp_obj;
        tmp_obj.neps = 1;
        tmp_obj.eps = eps;
        indices = new int[3];
        indices[0] = i_x;
        indices[1] = i_y;
        indices[2] = i_z;
        vector<MaterialObject>::iterator it_mat_vec = find(mat_vec.begin(),
            mat_vec.end(), tmp_obj);
        if (it_mat_vec == mat_vec.end()) {
          (tmp_obj.position).push_back(indices);
          mat_vec.push_back(tmp_obj);
        } else {
          (*it_mat_vec).neps += 1;
          ((*it_mat_vec).position).push_back(indices);
        }
        is_event++;
      }
    }

  }

  sort(mat_vec.begin(), mat_vec.end(), Sortneps());

  int check_sum = 0;
  vector<MaterialObject>::iterator it_mat_vec = mat_vec.begin();
  while (it_mat_vec != mat_vec.end()) {
    cout << (*it_mat_vec) << "\t node: "
        << (*it_mat_vec).position.at(0)[0]
            + (*it_mat_vec).position.at(0)[1] * m_nxyz[0]
            + (*it_mat_vec).position.at(0)[2] * m_nxyz[0] * m_nxyz[1] << endl;
    vector<int*>::iterator it = (*it_mat_vec).position.begin();
    check_sum += (*it_mat_vec).neps;
    while (it != (*it_mat_vec).position.end()) {
      delete (*it);
      it++;
    }
    it_mat_vec++;
  }
  printf("Analysed nodes: %d\n", is_event);
  printf(
      "Size of Material vecotor is (corresponds to different kinds of eps): %d\n",
      mat_vec.size());
  printf("Consitency check: nNodes = %d; sum_all_eps = %d\n",
      CST_Readout::GetNumberOfNodes(), check_sum);
  //Compare node material for one Element
  int element = 0;
  bool print_debug = false;
  vector<int> element_nodes;
  int count_elements = 0;
  vector<int>::iterator it_element_nodes;
  vector<double> all_eps;
  vector<double>::iterator it_all_eps;
  double *xyz = new double[3];
  for (int element = 0; element < nelements; element++) {
    CST_Readout::PrintStatus(is_event, nnodes + nelements);
    vector<double> eps_vec;
    vector<double>::iterator it_eps_vec;
    CST_Readout::GetNodesForElement(element, element_nodes);
    it_element_nodes = element_nodes.begin();
    while (it_element_nodes != element_nodes.end()) {
      int i, j, k;
      CST_Readout::Node2Index(*it_element_nodes, i, j, k);
      CST_Readout::GetPositionFromIndex(i, j, k, xyz);
      double eps = CST_Readout::Round(
          CST_Readout::GetMaterialPropertyFromIndex(i, j, k, 0), 2);
      it_eps_vec = find(eps_vec.begin(), eps_vec.end(), eps);
      if (it_eps_vec == eps_vec.end())
        eps_vec.push_back(eps);
      it_element_nodes++;
    }
    it_all_eps = find(all_eps.begin(), all_eps.end(), eps_vec.at(0));
    if (it_all_eps == all_eps.end())
      all_eps.push_back(eps_vec.at(0));
    if (eps_vec.size() != 1) {
      count_elements++;
      if (print_debug == true) {
        cout << "Vec size is: " << eps_vec.size() << endl;
        it_element_nodes = element_nodes.begin();
        while (it_element_nodes != element_nodes.end()) {
          int i, j, k;
          CST_Readout::Node2Index(*it_element_nodes, i, j, k);
          CST_Readout::GetPositionFromIndex(i, j, k, xyz);
          cout << "Element: " << element << "\t node: " << *it_element_nodes
              << "\t eps: "
              << CST_Readout::GetMaterialPropertyFromIndex(i, j, k, 0)
              << "\t at: (" << xyz[0] * scale << "," << xyz[1] * scale << ","
              << xyz[2] * scale << ")" << endl;
          it_element_nodes++;
        }
      }
    }
    is_event++;
  }
  cout << "I found " << all_eps.size() << " different eps and "
      << count_elements << " elements with different eps for their nodes!"
      << endl;
  delete xyz;
  xyz = 0;
  cout << "--->Leaving CST_Readout::AnalyseMaterial" << endl;
}
