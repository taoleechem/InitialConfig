/*
 * =====================================================================================
 *
 *       Filename:  molecule.h
 *
 *    Description:  Core design headfile for chemistry study
 *
 *        Version:  1.0
 *        Created:  2016年03月21日 14时39分14秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tao Li (), taoleechem@outlook.com
 *   Organization:  Nanjing University
 *
 * =====================================================================================
 */

#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <iomanip> //control the precision
#include <cmath>
#include <time.h>
using namespace std;
const int MaxAtom = 200;
const int MAXFRAGMENT = 8;
const double PI = 3.1415926;
const double HARTREE = 4.359744e-18;
const double R_BOLTZMAN = 8.314;
const double AVOGADRO_CONST = 6.023e23;
const double K_B_BOLTZMAN = 1.38037e-23;
const double ROOM_TEMPERATURE = 300;
const int MaxElementOrder=109;
const string ELEMENT[MaxElementOrder]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K",
                  "Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
                  "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs",
                  "Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
                  "W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa",
                  "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt"};
                  //109 ge yuan su
const double ELEMENT_RADIUS[MaxElementOrder]={0.30, 1.16, 1.23, 0.89, 0.88,0.77, 0.70, 0.66, 0.58, 0.55,1.40, 1.36, 1.25, 1.17, 1.10,
                    1.11, 0.99, 1.58, 2.03, 1.74,1.44, 1.32, 1.20, 1.13, 1.17,1.16, 1.16, 1.15, 1.17, 1.25,
                    1.25, 1.22, 1.21, 1.17, 1.14,1.89, 2.25, 1.92, 1.62, 1.45,1.34, 1.29, 1.23, 1.24, 1.25,
                    1.28, 1.34, 1.41, 1.50, 1.40,1.41, 1.37, 1.33, 2.09, 2.35,1.98, 1.69, 1.65, 1.65, 1.64,
                    1.64, 1.66, 1.85, 1.61, 1.59,1.59, 1.58, 1.57, 1.56, 1.70,1.56, 1.44, 1.34, 1.30, 1.28,
                    1.26, 1.26, 1.29, 1.34, 1.44,1.55, 1.54, 1.52, 1.53, 1.52,1.53, 2.45, 2.02, 1.70, 1.63,
                    1.46, 1.40, 1.36, 1.25, 1.57,1.58, 1.54, 1.53, 1.84, 1.61,1.50, 1.49, 1.38, 1.36, 1.26,
                    1.20, 1.16, 1.14, 1.06 };

template <class T> string X_ToStr(T tmp)
{
	stringstream ss;
	ss << tmp;
	return ss.str();
}
double RandomNumber(double MaxValue);


class Bond
{
public:
	int i, j;
	int order;
	Bond() :i(0), j(0), order(0) {}
	Bond(int ii, int ij, int io) :i(ii), j(ij), order(io) {}
	Bond(const Bond &ia)
	{
		i = ia.i;
		j = ia.j;
		order = ia.order;
	}
	void SetValue(int ii, int jj, int oor)
	{
		i = ii;
		j = jj;
		order = oor;
	}
	void IncreaseLabel(int num)
	{
		i += num;
		j += num;
	}
	friend ostream& operator<<(ostream &os, Bond &ia)
	{
		os << ia.i << " " << ia.j << " " << ia.order;
		return os;
	}
};

class Molecule
{
protected:
	int number;
	vector<string> name;
	double corr[MaxAtom][3];
public:
	//Formulation function
	Molecule();
	Molecule(const Molecule &ia);
	//Reload operators
	Molecule& operator=(const Molecule &ia);
	friend Molecule operator+(const Molecule &ia, const Molecule &ib);
	friend bool operator==(const Molecule &ia, const Molecule &ib);
	friend bool operator!=(const Molecule &ia, const Molecule &ib);
	friend bool operator>=(const Molecule &ia, const Molecule &ib);
	double AtomMass(const string &iname);

	string MoleculeName();
	friend ostream& operator<<(ostream &os, Molecule &ia);

	//Read&write files
	void ReadFromXYZfile(const string filename);
	void ReadFromXYZOnlyGeo(ifstream &infile, int atom_num);
	void ToXYZfile(const string filename);
	void ToXYZfileOnlyGeo(ofstream &tofile, bool judge);
	friend void ToXYZfile(const Molecule &a, const Molecule &b, string &filename, string other_info = "  Have a good day!");
	void ToMol2File(const string filename);
	void ToMol2File_AmBn(const string filename, int a_atoms, int a_numbers, int b_atoms);
	void ToMol2fileOnlyGeo(ofstream &tofile, int initial_label);
	friend void ToMol2fileAmBnType(const string filename,vector<Molecule> &As, vector<Molecule>& Bs);
	void ToNWchemFileHF(const string filename, const string basis);
	friend void ToNWchemFileHF(const Molecule &a, const Molecule &b, string &filename, const string basis = "6-31G");
	void ToNWchemFileDFT(const string filename, const string basis, const string functional);
	friend void ToNWchemFileDFT(const Molecule &a, const Molecule &b, string &filename, const string basis = "6-31G", const string functional = "b3lyp");
	void ToG09FileDFT(string &filename, string basis = "6-31g", string functional = "b3lyp");
	friend void ToG09FileDFT(Molecule &a, Molecule &b, string &filename, string basis = "6-31g", string functional = "b3lyp");
	void ReadFromGJF(string &filename, int atomNum);
	void ReadFromTinkerXYZfile(string filename);
	void ReadFromTinkerXYZGeoPart(ifstream &infile,int atom_number);
	void ToPDBfile(const string filename, int connection_matrix[][8]);
	void ToPDBfileOnlyGeo(ofstream &tofile, int initial_label);
	friend void ToPDBfileAmBnType(const string filename, vector<Molecule> &A, vector<Molecule> &B, int connect_m[][8], int connect_n[][8]);

	void AddAtom(const string atomname, double &x, double &y, double &z);
	//small functions
	int Number();
	double MoleculeMass();
	double AtomRadius(int i);
	int BondOrder(int i, int j);
	int ConnectAtoms(int i);
	vector<int> WhoConnectMe(int i, int except_label);
	string StandardAtomName(string x);
	Eigen::Vector3d MassCenter();
	friend double DistanceOfMassCenter(Molecule &ia, Molecule &ib);
	friend Molecule AbstractSomeFormNewMolecule(Molecule &ia, int begin_label,int end_label);

	void clear();
	//Simple Operate
	void PerformRot(Eigen::Matrix3d rot);
//Unfinished
	void PerformBondRot(int i, int j, double angle_radian)
	{
		cout << "Enter Bond Rot " << i << " " << j << " with angle " << int(angle_radian*180/PI) << endl;
		vector<int> RotAtoms(WhoConnectMe(j, i));
		//cout << RotAtoms[0] << endl;
		Eigen::Vector3d axis;
		axis << corr[j][0] - corr[i][0], corr[j][1] - corr[i][1], corr[j][2] - corr[i][2];
		axis.normalize();
		Eigen::AngleAxis<double> rot(angle_radian, axis);
		for (int i = 0; i < RotAtoms.size(); i++)
		{
			Eigen::Vector3d singleAtom;
			singleAtom << corr[RotAtoms[i]][0], corr[RotAtoms[i]][1], corr[RotAtoms[i]][2];
			singleAtom = rot*singleAtom;
			corr[RotAtoms[i]][0] = singleAtom(0);
			corr[RotAtoms[i]][1] = singleAtom(1);
			corr[RotAtoms[i]][2] = singleAtom(2);
		}
	}
	void PerformTrans(const Eigen::Vector3d trans);
	void PerformXTrans(const double &deltaX);
	void PerformYTrans(const double &deltaY);
	void PerformZTrans(const double &deltaZ);
	void PerformAxisRot(Eigen::Vector3d axis, double angle_radian);
	void PerformOnePointRotToXMinus(Eigen::Vector3d point);
	void PerformOnePointRotToXPlus(Eigen::Vector3d point);
	Eigen::Matrix3d EulerRot(double a1_radian, double a2_radian, double a3_radian);
	void PerformRandomRotEuler(Eigen::Vector3d Point, const double &Precision_degree);
	void PerformRandomRotEuler(Eigen::Vector3d Point, const double Precision_degree, double &ax_degree, double &ay_degree,double &az_degree);
	Eigen::Matrix3d EulerRotXY(double ax_radian, double ay_radian);
	void PerformRandomRotEulerXY(Eigen::Vector3d Point, const double &Precision_degree);
	void PerformRandomRotEulerXY(Eigen::Vector3d Point, const double &Precision_degree, double &ax_degree, double &ay_degree);


	void MCtoOrigin();
	void MCtoVector(const Eigen::Vector3d x);
	friend double ClosestDistance(Molecule &ia, Molecule &ib);
    friend double ClosestDistanceWithLabel(Molecule &ia, Molecule &ib, int &ia_label, int &ib_label);
    friend double AtomRadiusSum(Molecule &ia, Molecule &ib, int ia_label, int ib_label);
	friend void MakeAtomsSuitableDistanceMoveB(Molecule &ia, Molecule &ib, const double SmallestDistance);
    friend void MakeAtomsUniformDistanceMoveB(Molecule &ia, Molecule &ib, const double times);
	friend double RMSD(Molecule &ia, Molecule &ib);
	friend void SpaceTransform(Molecule ref, Molecule &change);
	void AligenToStandardConfig();
	friend void RandomRotPossibleBond(Molecule &ia, double Resolusion_radian)
	{
		//Firstly, judge which bond has potential to rot
		vector<Bond> rotable;
		for (int i = 0; i < ia.number; i++)
			for (int j = i + 1; j < ia.number; j++)
			{
				int order = ia.BondOrder(i, j);
				if (order == 1 && ia.ConnectAtoms(i)>1 && ia.ConnectAtoms(j)>1)
					rotable.push_back(Bond(i, j, order));
			}
		//Secondly, rot with this bond.
		int num = rotable.size();
		clock_t now = clock();
		srand(now);
		int Max = int(2 * PI / Resolusion_radian);
		double x = (rand() % (int)((Max)) + 1)*Resolusion_radian;
		for (int i = 0; i < num; i++)
		{
			ia.PerformBondRot(rotable[i].i, rotable[i].j, x);
		}
	}
};

class DoubleMolecule
{
private:
	Molecule a, b;
	double energy;
public:
	DoubleMolecule(Molecule &ia, Molecule &ib, double &ienergy);
	DoubleMolecule(const DoubleMolecule &i);
	DoubleMolecule();
	DoubleMolecule& operator=(DoubleMolecule &id);
	friend bool operator>(const DoubleMolecule &ia, const DoubleMolecule &ib);
	friend bool operator<(const DoubleMolecule &ia, const DoubleMolecule &ib);
	void Set(Molecule ia, Molecule ib, double ienergy);
	void GetInfo(Molecule &ia, Molecule &ib, double &ienergy);
	void ToXYZ(string filename);
	void output();
	double Energy();
	friend double RMSD(const DoubleMolecule ia, const DoubleMolecule ib);
};

class Fragments
{
protected:
	int frag_number;
	vector<Molecule> frags;
	int index;
public:
	Fragments();
	Fragments(Fragments &ia);
	friend ostream& operator<<(ostream &os, Fragments &ia);
	Fragments& operator++();
	void IndexToZero();
	Molecule ThisFragment();
	Molecule OneFragment(unsigned int &my_index);
	Molecule OtherFragments();
	Molecule TotalFragments();
	Eigen::Vector3d VectorFromThisToOther();
	int FragNumbers();
	void ReadFromXYZfile(const string filename, const int inumber, int matrix[][2]);
	//Geometry operation
	void PerformRot(Eigen::Matrix3d rot);
	void PerformTrans(const Eigen::Vector3d trans);
	void PerformXTrans(const double &deltaX);
	void PerformZTrans(const double &deltaZ);
	void PerformAxisRot(Eigen::Vector3d axis, double angle_radian);
	void PerformOnePointRotToXMinus(Eigen::Vector3d point);
	void PerformOnePointRotToXPlus(Eigen::Vector3d point);
};
class SolventCube
{
protected:
	Molecule solute;
	vector<Molecule> solvent;
	int solvent_num;
	double a,b,c,alpha,beta,gama;
public:
	SolventCube()
	{
		a =b=c= 0.0;
		alpha = beta = gama = 90.0;
		solvent_num = 0;
	}
	Molecule Solute()
	{
		return solute;
	}
	void ReadFromTinkerXYZfile(const string filename,int solute_atoms, int single_solvent_atoms)
	{
		ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << " to get the geometry info" << endl;
			exit(1);
		}
		int total_num;
		infile >> total_num;
		string temp;
		getline(infile, temp);
		infile >> a>>b>>c>>alpha>>beta>>gama;
		getline(infile, temp);
		
		solvent_num = (total_num - solute_atoms) / single_solvent_atoms;
		int inumber;
		string iname, iconnection;
		//input solute
		solute.clear();
		solute.ReadFromTinkerXYZGeoPart(infile,solute_atoms);
		//input solvents
		solvent.clear();
		for (int i = 0; i < solvent_num; i++)
		{
			Molecule temp_mole;
			temp_mole.ReadFromTinkerXYZGeoPart(infile, single_solvent_atoms);
			solvent.push_back(temp_mole);
		}
		infile.close();
	}
	void ToXYZfile(const string filename, bool IfOutputStandardAtomName)
	{
		ofstream tofile1(filename.c_str(), ios::out);
		if (!tofile1)
		{
			cerr << "Error to write " << filename << endl;
			exit(1);
		}
		tofile1.close();
		ofstream tofile(filename.c_str(), ios::app);
		tofile << solvent_num*solvent[0].Number() + solute.Number() << endl;
		tofile<<a<<"\t"<<b<<"\t"<<c<<"\t"<<alpha<<"\t"<<beta<<"\t"<<gama<< endl;
		solute.ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		for (int i = 0; i < solvent.size(); i++)
			solvent[i].ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		tofile.close();
	}
	//*_singal represents the direction you want to expand this cube to a larger 8 cubes. *_singal=+1 or -1, *_singal=0 means do nothing
	void ExpandCubeTo8(int x_singal,int y_singal,int z_singal)
	{
		if (abs(x_singal) == 1 && abs(y_singal) == 1 && abs(z_singal) == 1)
		{
			for (int i = 0; i != solvent_num; i++)
			{
				Molecule one_solvent[8];
				one_solvent[0] = solvent[i];
				for (int j = 1; j != 8; j++)
					one_solvent[j] = one_solvent[0];
				one_solvent[1].PerformXTrans(x_singal*a);
				one_solvent[2].PerformYTrans(y_singal*b);
				one_solvent[3].PerformXTrans(x_singal*a);
				one_solvent[3].PerformYTrans(y_singal*b);
				one_solvent[4].PerformZTrans(z_singal*c);
				one_solvent[5].PerformXTrans(x_singal*a);
				one_solvent[5].PerformZTrans(z_singal*c);
				one_solvent[6].PerformYTrans(y_singal*b);
				one_solvent[6].PerformZTrans(z_singal*c);
				one_solvent[7].PerformXTrans(x_singal*a);
				one_solvent[7].PerformYTrans(y_singal*b);
				one_solvent[7].PerformZTrans(z_singal*c);
				for (int j = 1; j != 8; j++)
					solvent.push_back(one_solvent[j]);
			}
			solvent_num = solvent_num * 8;
			a = a * 2;
			b = b * 2;
			c = c * 2;
		}
	}
	int CountSolventNumberNearSolute(double radius)
	{
		Eigen::Vector3d MC,mc;
		MC = solute.MassCenter();
		int count = 0;
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();

			double length = sqrt((MC(0)-mc(0))*(MC(0) - mc(0))+ (MC(1) - mc(1))*(MC(1) - mc(1))+ (MC(2) - mc(2))*(MC(2) - mc(2)));
			if (length <= radius)
				count += 1;
		}
		return count;
	}
	void ToXYZSolventNearSolute(const string filename,double radius,bool IfOutputStandardAtomName)
	{
		//count how many solvent molecules within
		int count = 0;
		Eigen::Vector3d MC, mc;
		MC = solute.MassCenter();
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			double length = sqrt((MC(0) - mc(0))*(MC(0) - mc(0)) + (MC(1) - mc(1))*(MC(1) - mc(1)) + (MC(2) - mc(2))*(MC(2) - mc(2)));
			if (length <= radius)
				count += 1;
		}

		ofstream tofile1(filename.c_str(), ios::out);
		if (!tofile1)
		{
			cerr << "Error to write " << filename << endl;
			exit(1);
		}
		tofile1.close();
		ofstream tofile(filename.c_str(), ios::app);
		tofile << count*solvent[0].Number() + solute.Number() << endl;
		tofile << a << "\t" << b << "\t" << c << "\t" << alpha << "\t" << beta << "\t" << gama << endl;
		solute.ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		for (int i = 0; i < solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			double length = sqrt((MC(0) - mc(0))*(MC(0) - mc(0)) + (MC(1) - mc(1))*(MC(1) - mc(1)) + (MC(2) - mc(2))*(MC(2) - mc(2)));
			if (length <= radius)
				solvent[i].ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		}
		tofile.close();
	}
	void AppendToXYZSolventNearSolute(ofstream &tofile, double radius,bool IfOutputStandardAtomName)
	{
		//count how many solvent molecules within
		int count = 0;
		Eigen::Vector3d MC, mc;
		MC = solute.MassCenter();
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			double length = sqrt((MC(0) - mc(0))*(MC(0) - mc(0)) + (MC(1) - mc(1))*(MC(1) - mc(1)) + (MC(2) - mc(2))*(MC(2) - mc(2)));
			if (length <= radius)
				count += 1;
		}

		tofile << count*solvent[0].Number() + solute.Number() << endl;
		tofile << a << "\t" << b << "\t" << c << "\t" << alpha << "\t" << beta << "\t" << gama << endl;
		solute.ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		for (int i = 0; i < solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			double length = sqrt((MC(0) - mc(0))*(MC(0) - mc(0)) + (MC(1) - mc(1))*(MC(1) - mc(1)) + (MC(2) - mc(2))*(MC(2) - mc(2)));
			if (length <= radius)
				solvent[i].ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		}
	}
	int CountSolventNumberInSolute(double MCtoMarginLength)
	{
		int count = 0;
		Eigen::Vector3d MC, mc;
		MC = solute.MassCenter();
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			if (abs(mc(0) - MC(0)) < MCtoMarginLength && abs(mc(1) - MC(1)) < MCtoMarginLength && abs(mc(2) - MC(2)) < MCtoMarginLength)
			{
				count += 1;
			}
		}
		return count;
	}
	void ToXYZSolventInSolute(const string filename, double MCtoMarginLength,bool IfOutputStandardAtomName)
	{
		//count how many solvent molecules within
		int count = 0;
		Eigen::Vector3d MC, mc;
		MC = solute.MassCenter();
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			if (abs(mc(0) - MC(0)) < MCtoMarginLength && abs(mc(1) - MC(1)) < MCtoMarginLength && abs(mc(2) - MC(2)) < MCtoMarginLength)
			{
				count += 1;
			}
		}

		ofstream tofile1(filename.c_str(), ios::out);
		if (!tofile1)
		{
			cerr << "Error to write " << filename << endl;
			exit(1);
		}
		tofile1.close();
		ofstream tofile(filename.c_str(), ios::app);
		tofile << count*solvent[0].Number() + solute.Number() << endl << endl;
		solute.ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		for (int i = 0; i < solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			if (abs(mc(0) - MC(0)) < MCtoMarginLength && abs(mc(1) - MC(1)) < MCtoMarginLength && abs(mc(2) - MC(2)) < MCtoMarginLength)
			{
				solvent[i].ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
			}
		}
		tofile.close();
	}
	void AppendToXYZSolventInSolute(ofstream &tofile, double MCtoMarginLength,bool IfOutputStandardAtomName)
	{
		//count how many solvent molecules within
		int count = 0;
		Eigen::Vector3d MC, mc;
		MC = solute.MassCenter();
		for (int i = 0; i != solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			if (abs(mc(0) - MC(0)) < MCtoMarginLength && abs(mc(1) - MC(1)) < MCtoMarginLength && abs(mc(2) - MC(2)) < MCtoMarginLength)
			{
				count += 1;
			}
		}
		tofile << count*solvent[0].Number() + solute.Number() << endl;
		tofile << a << "\t" << b << "\t" << c << "\t" << alpha << "\t" << beta << "\t" << gama << endl;
		solute.ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
		for (int i = 0; i < solvent_num; i++)
		{
			mc = solvent[i].MassCenter();
			if (abs(mc(0) - MC(0)) < MCtoMarginLength && abs(mc(1) - MC(1)) < MCtoMarginLength && abs(mc(2) - MC(2)) < MCtoMarginLength)
			{
				solvent[i].ToXYZfileOnlyGeo(tofile, IfOutputStandardAtomName);
			}
		}
	}
	void ReadFromTinkerArcGetXYZonce(ifstream &infile, int &solute_atoms, int &single_solvent_atoms)
	{
			int total_num;
			infile >> total_num;
			string temp;
			getline(infile, temp);
			infile >> a>>b>>c>>alpha>>beta>>gama;//input a
			getline(infile, temp);

			solvent_num = (total_num - solute_atoms) / single_solvent_atoms;//input solvent numbers
			int inumber;
			string iname, iconnection;
			//input solute
			solute.clear();
			solute.ReadFromTinkerXYZGeoPart(infile, solute_atoms);
			//input solvents
			solvent.clear();
			for (int i = 0; i < solvent_num; i++)
			{
				Molecule temp_mole;
				temp_mole.ReadFromTinkerXYZGeoPart(infile, single_solvent_atoms);
				solvent.push_back(temp_mole);
			}
	}
};

//Arc.analysis
void ReadFromWholeTinkerArc(const string arc_filename, const string save_filename, int solute_atoms, int each_solvent_atoms, int x_singal, int y_singal, int z_singal, double radius, double SoluteCenterToMarginLength,bool IfOutputStandardAtomName);
void Do_ReadFromWholeTinkerArc_FromTxt();
//GenerateFunction() result analysis multi-xyz file
void AligenMultiXYZfileToOneKeepMoleculeAStill(const string dir_name, int total_file_num, int molecule_A_atoms);
void Do_AligenMultiXYZ_Program();
void AlignEachXYZToStandardForm(const string dir_name, int total_file_num, int molecule_A_atoms);
void Do_AligenXYZStandard_Program();

//xyz-->pdb
void XYZToPDB_MoleculeAmBnType(const string xyz_filename, const string save_filename, int a_atoms, int a_num, int b_atoms, int b_num, int connect_m[][8], int connect_n[][8]);
void ReadConnectionInfo(const string filename, int connect[][8], int atoms);
void Do_XYZToPDB_MoleculeAmBnType();

//xyz-->mol2
void XYZToMol2_MoleculeAmBnType(const string xyz_filename, const string save_filename, int a_atoms, int a_num, int b_atoms);
void Do_XYZToMol2_MoleculeAmBnType();
#endif

