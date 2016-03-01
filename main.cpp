#include "MyMolecules.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <iomanip> //control the precision
#include <cmath>
#include <time.h>
using namespace std;
/*
static void outputWelcome()
{
cout << "##################################################" << endl;
cout << "#      Welcome to use Molecule Recognition       #" << endl;
cout << "# Many stable configurations will be produced    #" << endl;
cout << "##################################################" << endl;
}
int main()
{
outputWelcome();
DoubleMolecule geoinfo;
ifstream thisfile("task.txt",ios::in);
string ifile; int a1,b1;
thisfile>>ifile>>a1>>b1;
double StepLength = 4;
double StepPrecision = 0.15;
double MaxRotTime = 20000;
double RotPrecision = 20;
thisfile>>StepLength>>StepPrecision>>MaxRotTime>>RotPrecision;
thisfile.close();
geoinfo.ReadFromTinkerXYZ("InitiConfig/"+ifile,a1,b1);
Molecule a, b;
geoinfo.GetABInfo(a, b);
a.PerformRandomRotEuler(10);
a.output();
b.PerformRandomRotEuler(10);
b.output();
string forcefieldORbasis = "oplsaa.prm";
//This part is suitable for single test
string MinConfigName = "SaveConfigs/Min_";
string RelativeMinConfigName = "SaveConfigs/RelaMin_";
MonteCarlo(a, b, StepLength, StepPrecision, MaxRotTime, RotPrecision, forcefieldORbasis, MinConfigName, RelativeMinConfigName);

return 0;
}

*/
void GenerateFunction(string x_filename, const int x_number, int x_matrix[][2], string y_filename, const int y_number, int y_matrix[][2], int OutputNumbers, string Basis);

int main()
{
	int matrix[8][2] = { {0,2},{3,4},{5,7},{8,10},{11,12},{13,14},{15,17},{18,20} };
	int index = 8;
	int matrix2[2][2] = { {0,0},{1,1} };
	int index2 = 2;
	Fragments A,B;
	A.ReadFromXYZfile("example.xyz", index, matrix);
	B.ReadFromXYZfile("HCl.xyz", index2, matrix2);
    Molecule A1 = A.ThisFragment();
	Molecule A2 = A.OtherFragments();
	Molecule B1 = B.ThisFragment();
	Molecule B2 = B.OtherFragments();
	Eigen::Vector3d MC_A1 = A1.MassCenter();
	A1.PerformTrans(-1 * MC_A1);  A2.PerformTrans(-1 * MC_A1);
	Eigen::Vector3d MC_A2 = A2.MassCenter();
	A1.PerformPointRotToXMinus(MC_A2);  A2.PerformPointRotToXMinus(MC_A2);//这里有问题
	cout << A1.MassCenter() << endl;
	cout<< A2.MassCenter() << endl;
	


	// GenerateFunction("example.xyz",index,matrix,"HCl.xyz",index2,matrix2,100,"6-31g");
	system("PAUSE");
	return 0;
}