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
const double PI = 3.1415926;

template <class T> string X_ToStr(T tmp)
{
	stringstream ss;
	ss << tmp;
	return ss.str();
}
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
	Eigen::Vector3d MassCenter();
	friend double DistanceOfMassCenter(Molecule &ia, Molecule &ib);
	friend Molecule AbstractSomeFormNewMolecule(Molecule &ia, int begin_label,int end_label);

	void clear();
	//Simple Operate
	void PerformRot(Eigen::Matrix3d rot);
	void PerformTrans(const Eigen::Vector3d trans);
	void PerformXTrans(const double &deltaX);
	void PerformYTrans(const double &deltaY);
	void PerformZTrans(const double &deltaZ);
	void PerformAxisRot(Eigen::Vector3d axis, double angle_radian);
	void PerformOnePointRotToXMinus(Eigen::Vector3d point);
	Eigen::Matrix3d EulerRot(double &a1, double &a2, double &a3);
	void PerformRandomRotEuler(Eigen::Vector3d Point, const double &Precision_degree);
	void MCtoOrigin();
	void MCtoVector(const Eigen::Vector3d x);
	friend double ClosestDistance(Molecule &ia, Molecule &ib);
	friend void MakeAtomsSuitableDistanceMoveB(Molecule &ia, Molecule &ib, const double SmallestDistance);
	friend double RMSD(Molecule &ia, Molecule &ib);
	friend void SpaceTransform(Molecule ref, Molecule &change);
	void AligenToStandardConfig();
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
	friend bool operator>(DoubleMolecule &ia, DoubleMolecule &ib);
	friend bool operator<(DoubleMolecule &ia, DoubleMolecule &ib);
	void Set(Molecule ia, Molecule ib, double ienergy);
	void GetInfo(Molecule &ia, Molecule &ib, double &ienergy);
	void ToXYZ(string filename);
	void output();
	double Energy();
	friend double RMSD(DoubleMolecule &ia, DoubleMolecule &ib);
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
#endif

