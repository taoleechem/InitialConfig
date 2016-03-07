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
	void ToXYZfile(const string filename);
	friend void ToXYZfile(const Molecule &a, const Molecule &b, string &filename, string other_info = "  Have a good day!");
	void ToNWchemFileHF(const string filename, const string basis = "6-31G");
	friend void ToNWchemFileHF(const Molecule &a, const Molecule &b, string &filename, const string basis = "6-31G");
	void ToG09FileDFT(string &filename, string basis = "6-31g", string functional = "b3lyp");
	friend void ToG09FileDFT(Molecule &a, Molecule &b, string &filename, string basis = "6-31g", string functional = "b3lyp");
	void ReadFromGJF(string &filename, int atomNum);

	void AddAtom(const string atomname, double &x, double &y, double &z);
	//small functions
	int Number();
	double MoleculeMass();
	Eigen::Vector3d MassCenter();
	friend double DistanceOfMassCenter(Molecule &ia, Molecule &ib);

	void clear();
	//Simple Operate
	void PerformRot(Eigen::Matrix3d rot);
	void PerformTrans(const Eigen::Vector3d trans);
	void PerformXTrans(const double &deltaX);
	void PerformZTrans(const double &deltaZ);
	void PerformAxisRot(Eigen::Vector3d axis, double angle_radian);
	void PerformOnePointRotToXMinus(Eigen::Vector3d point);
	Eigen::Matrix3d EulerRot(double &a1, double &a2, double &a3);
	void PerformRandomRotEuler(Eigen::Vector3d Point, const double &Precision_degree);
	void MCtoOrigin();
	void MCtoVector(const Eigen::Vector3d x);
	friend void MakeAtomsShortestDistanceMoveB(Molecule &ia, Molecule &ib, const double SmallestDistance);
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
};

#endif

