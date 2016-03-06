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
int CallTimes = 0;

double ReadFile(string Tempfilename)
{
	ifstream readfile(Tempfilename.c_str());
	if (!readfile)
	{
		cerr << "Error to read " << Tempfilename << endl;
		exit(1);
	}
	double num;
	readfile >> num;
	readfile.close();
	return num;
}
double G09energy(Molecule &a, string basis = "6-31g", string functional = "b3lyp")
{
#ifdef _WIN32
	clock_t now = clock();
	srand(now);
	double x=(rand() % 1000 + 1);
	return x;
#else
	string filename("DATA/system");
	a.ToG09FileDFT(filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
#endif

}
double G09energy(Molecule &a, Molecule &b, string basis = "6-31g", string functional = "b3lyp")
{
#ifdef _WIN32
	clock_t now = clock();
	srand(now);
	double x= (rand() % 1000 + 1);
	return x;
#else
	string filename("DATA/system");
	ToG09FileDFT(a, b, filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
#endif
}


/*
static int MonteCarlo01Distribution(double delta_potential, double T)
{
double probability = exp(-abs(delta_potential) * 315772 / T);
int MaxNums = (int)(1 / probability);
clock_t now = clock();

//std::default_random_engine generator(now);
//std::uniform_int_distribution<int> dis(1, MaxNums);
//if (dis(generator) == 1)

srand(now);
if ((rand() % MaxNums + 1) == 1)
return 1;
else
return 0;
}
static void outputWelcome()
{
	cout << "##################################################" << endl;
	cout << "#      Welcome to use Monte Carlo Gaussian09         #" << endl;
	cout << "# Many first-opt configurations will be produced  #" << endl;
	cout << "##################################################" << endl;
}


static void StepOpt(Molecule a, Molecule b, double &relativeMinPotential, DoubleMolecule &relativeMinConfig, double potential, double ZeroEnergy, string forcefieldORbasis, double StepLength, double StepPrecision)
{
	cout << "Enter Step-opt step" << endl;
	b.PerformXTrans(StepLength / 2);
	relativeMinPotential = potential;//relative** is always in step-opt
	relativeMinConfig.Set(a, b, relativeMinPotential);
	const int stepCountTimes = (int)(StepLength / StepPrecision) + 1;
	for (int iStep = 0; iStep != stepCountTimes; iStep++)
	{
		b.PerformXTrans(-1 * StepPrecision);
		potential = G09energy(a, b, forcefieldORbasis) - ZeroEnergy;
		if (potential < relativeMinPotential)
		{
			relativeMinPotential = potential;
			relativeMinConfig.Set(a, b, relativeMinPotential);
		}
	}
}
static string CombineFileName(string fileName, double num)
{
	stringstream is;
	string IS;
	is << num;
	is >> IS;
	fileName = fileName + IS + ".xyz";
	return fileName;
}
void MonteCarlo(Molecule &a, Molecule &b, vector<DoubleMolecule> &SaveMinConfig, Eigen::Vector3d bposition, string forcefieldORbasis, string MinConfigName, string RelativeMinConfigName)
{
//initi parameter 预设A和B都已经在合适的位置
const double StepPrecision = 0.15;
const double StepLength = 3.5;
const int stepCountTimes = (int)(StepLength / StepPrecision) + 1;
const double RotPrecision = 90;
const double MaxRotTime = 100;
const double OriginTemp = 200;
double Temp = OriginTemp;
double dT = 10;
int Tcount = 0;
const double FirstDistance = 8;

int MinConfigNum = 0;
int RelativeMinConfigNum = 0;
//const string forcefieldORbasis("6-31g");
//string MinConfigName = "SaveConfigs/Min_";
//string RelativeMinConfigName = "SaveConfigs/RelaMin_";

outputWelcome();
//Calculate zeroEnergy
double zeroEnergy = 0;
zeroEnergy += G09energy(a, forcefieldORbasis);
zeroEnergy += G09energy(b, forcefieldORbasis);
const double ZeroEnergy = zeroEnergy;


//initilize min value with first step opt
int firstStepOptTimes = (int)((FirstDistance - 1) / StepPrecision);
//Calculate initi potential
double potential = 0;
potential = G09energy(a, b, forcefieldORbasis) - ZeroEnergy;
double relativeMinPotential = potential;//relative*** is in the step-opt procedure
DoubleMolecule relativeMinConfig(a, b, relativeMinPotential);
for (int i = 0; i != firstStepOptTimes; i++)
{
b.PerformXTrans(-1 * StepPrecision);
potential = G09energy(a, b, forcefieldORbasis) - ZeroEnergy;
if (potential < relativeMinPotential)
{
relativeMinPotential = potential;
relativeMinConfig.Set(a, b, relativeMinPotential);
}
}

//Set global Min config & potential
double MinPotential = relativeMinPotential;
DoubleMolecule MinConfig(relativeMinConfig);
MinConfig.ToXYZ("init.txt", CallTimes); MinConfigNum += 1;//Save this one
cout << "No." << MinConfigNum << " is obtained " << endl;
MinConfig.output();



//Random-Rot

for (int i = 0; i != MaxRotTime; i++)
{
Temp += dT;
//This is important for each loop
MinConfig.GetInfo(a, b, potential);//Form global min config, to perform rot
b.PerformRandomRotEuler(bposition, RotPrecision);
potential = G09energy(a, b, forcefieldORbasis) - ZeroEnergy;
if (MonteCarlo01Distribution(potential - MinPotential, Temp) != 0)
{
cout << "Out rot permitted!, T is " << Temp << endl;
//step-opt
StepOpt(a, b, relativeMinPotential, relativeMinConfig, potential, ZeroEnergy, forcefieldORbasis, StepLength, StepPrecision);
relativeMinConfig.ToXYZ(CombineFileName(RelativeMinConfigName, RelativeMinConfigNum), CallTimes); RelativeMinConfigNum += 1;//Save this relative one
cout << "No." << RelativeMinConfigNum << " Relative MinConfig in outer loop and potential is " << relativeMinPotential << endl;
//If after-rot step-opt config has lower energy, save as global min
if (relativeMinPotential < MinPotential)
{
MinPotential = relativeMinPotential;
MinConfig = relativeMinConfig;
//save MIn
Temp = OriginTemp;
MinConfig.ToXYZ(CombineFileName(MinConfigName, MinConfigNum), CallTimes); MinConfigNum += 1;//Save this one
SaveMinConfig.push_back(MinConfig);
cout << "No." << MinConfigNum << " MinConfig in outer loop " << endl;
MinConfig.output();
}
//If not, keep rotation until forbidden!
else
{
double BranchTemp = Temp;
for (int iBranch = 0; iBranch != MaxRotTime; iBranch++)
{
//This is Important for each loop
relativeMinConfig.GetInfo(a, b, potential);//set a,b to a step-opt state
relativeMinPotential = potential;
//perform random rot
if (BranchTemp>(0.8*OriginTemp))
(BranchTemp -= 2 * dT);
b.PerformRandomRotEuler(bposition, RotPrecision);
potential = G09energy(a, b, forcefieldORbasis) - ZeroEnergy;

if (MonteCarlo01Distribution(potential - relativeMinPotential, BranchTemp) != 0)
{
//do step-opt
cout << "Inner rot permitted, Branch T is " << BranchTemp << endl;
StepOpt(a, b, relativeMinPotential, relativeMinConfig, potential, ZeroEnergy, forcefieldORbasis, StepLength, StepPrecision);
relativeMinConfig.ToXYZ(CombineFileName(RelativeMinConfigName, RelativeMinConfigNum), CallTimes); RelativeMinConfigNum += 1;//Save this relative one
cout << "No." << RelativeMinConfigNum << " Relative MinConfigNum in inner loop and potential is " << relativeMinPotential << endl;
//If relativeMin is global min, save it.
if (relativeMinPotential <MinPotential)
{
Temp = OriginTemp;
MinPotential = relativeMinPotential;
MinConfig = relativeMinConfig;
MinConfig.ToXYZ(CombineFileName(MinConfigName, MinConfigNum), CallTimes); MinConfigNum += 1;//Save this one
SaveMinConfig.push_back(MinConfig);
cout << "No." << MinConfigNum << " MinConfigNum in inner loop " << endl;
MinConfig.output();
break;
}
}
else
{
cout << "Inner rot forbidden!!, T is " << Temp << endl;
break;
}
}
}
}
else
{
cout << "Out random rot forbidden!! T is " << Temp << endl;
}
}
}

void OutputSomeConfig(vector<DoubleMolecule> &a, int numbers, string dir_name = "Result")
{
int total_num = a.size();
int MinIndex = 0, Index = 0;
DoubleMolecule temp;
for (int i = 0; i != total_num; i++)
for (int j = i; j != total_num; j++)
{
if (a[j].EnergyValue() < a[i].EnergyValue())
{
temp = a[j]; a[j] = a[i]; a[i] = temp;
}
}
for (int i = 0; i < numbers || i < total_num; i++)
{
a[i].ToXYZ(dir_name+"/PossibleConfig_"+X_ToStr<int>(i)+".xyz");
}
}

*/

void GenerateFunction()
{
	int matrix[3][2] = { { 0,1 },{ 2,4 },{ 5,6 } };
	int index = 3;
	int matrix2[2][2] = { { 0,2 },{ 3,3 } };
	int index2 = 2;
	const int MaxRotTimes = 50;
	const double RotPrecision = 60.0;
	const double B1_default_value = 3.00;
	const int EachPairSaveNumber = 4;
	const int OutPutNumber = 50;
	Fragments FA, FB;
	FA.ReadFromXYZfile("ch2choh_1.xyz", index, matrix);
	FB.ReadFromXYZfile("ch2o.xyz", index2, matrix2);

	double RestEnergies = 0;
	RestEnergies = G09energy(FA.TotalFragments()) + G09energy(FB.TotalFragments());
	vector<DoubleMolecule> SaveSuitableCofigs;

	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			cout << "Do calculate of " << FA.ThisFragment().MoleculeName() << " No." << i << ", and " << FB.ThisFragment().MoleculeName() << " No." << j << endl;
			Molecule A1 = FA.ThisFragment();
			Molecule A2 = FA.OtherFragments();
			Molecule B1 = FB.ThisFragment();
			Molecule B2 = FB.OtherFragments();
			//translate A1-MC to origin
			Eigen::Vector3d MC_A1 = A1.MassCenter();
			A1.PerformTrans(-1 * MC_A1);  A2.PerformTrans(-1 * MC_A1);
			//keep A1-MC to origin, put A2-MC to x- axis
			Eigen::Vector3d MC_A2 = A2.MassCenter();
			A1.PerformOnePointRotToXMinus(MC_A2);  A2.PerformOnePointRotToXMinus(MC_A2);
			//translate B1-MC to x+ axis(default distance is 3.00), B2 at random position
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;//3.0, 0, 0
			Eigen::Vector3d MC_B1 = B1.MassCenter();
			B1.PerformTrans(-1 * MC_B1 + Default_B1); B2.PerformTrans(-1 * MC_B1 + Default_B1);
			//rot B at point B1_MC randomly
			Molecule A = A1 + A2;
			Molecule B = B1 + B2;
			MC_B1 = B1.MassCenter();
			DoubleMolecule TempConfigs[MaxRotTimes];
			for (int k = 0; k != MaxRotTimes; k++)
			{
				B.PerformRandomRotEuler(MC_B1, RotPrecision);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsShortestDistanceMoveB(A, B, B1_default_value);
				double  potential = G09energy(A, B) - RestEnergies;
				TempConfigs[k].Set(A, B, potential);
			}
			//Sort configurations 
			DoubleMolecule temp(TempConfigs[0]);
			for (int ii = 0; ii != MaxRotTimes - 1; ii++)
				for (int jj = i + 1; jj != MaxRotTimes; jj++)
				{
					if (TempConfigs[ii].Energy() > TempConfigs[jj].Energy())
					{
						temp = TempConfigs[ii];
						TempConfigs[ii] = TempConfigs[jj];
						TempConfigs[jj] = temp;
					}
				}
			//Save Least Energy 4 configs
			for (int ii = 0; ii != EachPairSaveNumber; ii++)
				SaveSuitableCofigs.push_back(TempConfigs[ii]);
		}
	}

	//sort SaveSuitableConfigs
	DoubleMolecule temp;
	int total = SaveSuitableCofigs.size();
	for (int i = 0; i < total-1; i++)
		for (int j = i + 1; j < total; j++)
		{
			if (SaveSuitableCofigs[i] > SaveSuitableCofigs[j])
			{
				temp = SaveSuitableCofigs[i];
				SaveSuitableCofigs[i] = SaveSuitableCofigs[j];
				SaveSuitableCofigs[j] = temp;
			}
		}
	//output some configs
	cout << "Here output " << ((total > OutPutNumber) ? OutPutNumber : total) << " .xyz files to ./SaveConfigs/ as the final result" << endl;
	for (int i = 0; i < SaveSuitableCofigs.size() && i < OutPutNumber; i++)
		SaveSuitableCofigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
}


