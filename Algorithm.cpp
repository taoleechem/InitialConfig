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

static double ReadFile(string Tempfilename)
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
double G09energy(Molecule a, string basis = "6-31g", string functional = "b3lyp")
{
	string filename("DATA/system");
	a.ToG09FileDFT(filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
}
double G09energy(Molecule a, Molecule b, string basis = "6-31g", string functional = "b3lyp")
{
	string filename("DATA/system");
	ToG09FileDFT(a, b, filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
}
static int MonteCarlo01Distribution(double delta_potential, double T)
{
	double probability = exp(-abs(delta_potential) * 315772 / T);
	int MaxNums = (int)(1 / probability);
	clock_t now = clock();
	/*
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(1, MaxNums);
	if (dis(generator) == 1)
	*/
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
void GenerateFunction(string x_filename,const int x_number, int x_matrix[][2], string y_filename, const int y_number, int y_matrix[][2], int OutputNumbers,string Basis)
{
	const double B1MC_XValue = 5.0;
	const double PrecisionDegree = 90;
	vector<DoubleMolecule> MinConfig;
	Fragments A,B;
	A.ReadFromXYZfile(x_filename, x_number, x_matrix);
	B.ReadFromXYZfile(y_filename, y_number, y_matrix);
	for (int i = 0; i != A.FragNumbers(); i++, ++A)
	{
		for (int j = 0; j != B.FragNumbers(); j++, ++B)
		{
			Molecule A1 = A.ThisFragment();
			Molecule A2 = A.OtherFragments();
			Molecule B1 = B.ThisFragment();
			Molecule B2 = B.OtherFragments();
			Eigen::Vector3d MC_A1 = A1.MassCenter();
			//将A1质心在原点，A2质心转到x-
			A1.PerformTrans(-1 * MC_A1);  A2.PerformTrans(-1 * MC_A1);
			Eigen::Vector3d MC_A2 = A2.MassCenter();
			A1.PerformPointRotToXMinus(MC_A2);  A2.PerformPointRotToXMinus(MC_A2);
			//将B1质心在x+
			Eigen::Vector3d TransB1 = -1 * B1.MassCenter();
			TransB1(0) += B1MC_XValue;
			B1.PerformTrans(TransB1);  B2.PerformTrans(TransB1);
			//保持A不动 随机转动B 绕B1质心随机转动, 优化后，输出最佳构型；
			Molecule a = A1 + A2;
			Molecule b = B1 + B2;
			Eigen::Vector3d B_RotCenterPoint = B1.MassCenter();
			MonteCarlo(a, b, MinConfig, B_RotCenterPoint, Basis, "SaveConfigs/" + X_ToStr<int>(i) + "_" + X_ToStr<int>(j) + "_Min", "SaveConfigs/" + X_ToStr<int>(i) + "_" + X_ToStr<int>(j) + "_RelativeMin");
			OutputSomeConfig(MinConfig, OutputNumbers);

		}
		B.IndexToZero();
	}
}


