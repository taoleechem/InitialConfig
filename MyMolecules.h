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
	Molecule()
	{
		number = 0;
		for (int i = 0; i != MaxAtom; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = -10242048;
	}
	Molecule(const Molecule &ia)
	{
		number = ia.number;
		name = ia.name;
		for (int i = 0; i != number; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = ia.corr[i][j];
	}
    //Reload operators
	Molecule& operator=(const Molecule &ia)
	{
		number = ia.number;
		name = ia.name;
		for (int i = 0; i != number; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = ia.corr[i][j];
		return *this;
	}
	friend Molecule operator+(const Molecule &ia,const Molecule &ib)
	{
		Molecule temp;
		temp.number = ia.number + ib.number;
		for (int i = 0; i != ia.number; i++)
		{
			temp.name.push_back(ia.name[i]);
			for (int j = 0; j != 3; j++)
				temp.corr[i][j] = ia.corr[i][j];
		}
		for (int i = 0; i != ib.number; i++)
		{
			temp.name.push_back(ib.name[i]);
			for (int j = 0; j != 3; j++)
				temp.corr[i+ia.number][j] = ib.corr[i][j];
		}
		return temp;
	}
	friend bool operator==(const Molecule &ia, const Molecule &ib)
	{
		if (ia.number == ib.number&&ia.name == ib.name)
			return true;
		else return false;
	}
	friend bool operator!=(const Molecule &ia, const Molecule &ib)
	{
		if (ia.number != ib.number||ia.name != ib.name)
			return true;
		else return false;
	}
	friend bool operator>=(const Molecule &ia, const Molecule &ib)
	{
		if (ia.number >= ib.number)
			return true;
		else 
			return false;
	}
	const double AtomMass(const string &iname)
	{
		if (iname == "H")
			return 1.007825;
		else if (iname == "C")
			return 12.00;
		else if (iname == "N")
			return 14.003070;
		else if (iname == "O")
			return 15.994910;
		else if (iname == "S")
			return 31.972070;
		else if (iname == "Cl")
			return 34.968850;
		else
			return 0.00;
	}
	bool IsBiggerAtom(const string &ia, const string &ib)
	{
		if (AtomMass(ia) > AtomMass(ib))
			return true;
		else
			return false;
	}
	string MoleculeName()
	{
		vector<string> new_name(name);
		sort(new_name.begin(),new_name.end());
		vector<string>::iterator end_unique = unique(new_name.begin(), new_name.end());
		new_name.erase(end_unique, new_name.end());
		vector<string>::iterator iter1,iter2;
		int count=0;
		string atom_name;
		for (iter1 = new_name.begin(); iter1 != new_name.end(); iter1++)
		{
			atom_name += *iter1;
			for (iter2 = name.begin(); iter2 != name.end(); iter2++)
			{
				if (*iter2 == *iter1)
					count += 1;
			}
			if(count>=2)
				atom_name += X_ToStr<int>(count);
			count = 0;
		}
		return atom_name;
	}
	friend ostream& operator<<(ostream &os, Molecule &ia)
	{
		os << ia.number << endl;
		os << "standard .xyz file  "<<ia.MoleculeName()<<endl;
		for (int i = 0; i != ia.number; i++)
		{
			os << ia.name[i];
			for (int j = 0; j != 3; j++)
			{
				if (abs(ia.corr[i][j]) < 1e-7)
					ia.corr[i][j] = 0;
				os <<"\t"<< fixed << setprecision(7)<< ia.corr[i][j];
			}
			os << endl;
		}
		return os;
	}

	//Read&write files
	void ReadFromXYZfile(const string filename)
	{
		ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << " to get the geometry info" << endl;
			exit(1);
		}
		infile >> number;
		//Clear and initialize the info
		string temp;
		getline(infile, temp);
		getline(infile, temp);
		string iname;
		double icorr;
		for (int i = 0; i != number; i++)
		{
			infile >> iname;
			name.push_back(iname);
			for (int j = 0; j != 3; j++)
			{
				infile >> icorr;
				corr[i][j] = icorr;
			}
		}
		infile.close();
	}
	void ToXYZfile(const string filename)
	{
		ofstream tofile(filename.c_str(), ios::out);
		if (!tofile)
		{
			cerr << "Error to write " << filename << endl;
			exit(1);
		}
		if (number== 0)
			cout << "Empty molecule and no info is written to " << filename << endl;
		else
		{
			tofile << number << endl << endl;
			for (int i = 0; i != number; i++)
			{
				tofile << name[i];
				for (int j = 0; j != 3; j++)
				{
					if (abs(corr[i][j]) < 1e-7)
						corr[i][j] = 0;
					tofile << "\t" << fixed << setprecision(7) << corr[i][j];
				}
				tofile << endl;
			}
		}
		tofile.close();
	}
	friend void ToXYZfile(const Molecule &a, const Molecule &b, string &filename, string other_info="  Have a good day!")
	{
		ofstream tofile(filename.c_str(), ios::out);
		if (!tofile)
		{
			cerr << "Error to write " << filename << endl;
			exit(1);
		}
		int atomNum = a.number + b.number;
		if (atomNum == 0)
			cout << "Empty molecule and no info is written to " << filename << endl;
		else
		{
			tofile << atomNum << endl;
			tofile << other_info << endl;
			for (unsigned int i = 0; i != a.number; i++)
			{
				tofile << a.name[i];
				for (int j = 0; j != 3; j++)
				{
					tofile << "\t" << fixed << setprecision(7) << a.corr[i][j];
				}
				tofile << endl;
			}
			for (unsigned int i = 0; i != b.number; i++)
			{
				tofile << b.name[i];
				for (int j = 0; j != 3; j++)
				{
					tofile << "\t" << fixed << setprecision(7) << b.corr[i][j];
				}
				tofile << endl;
			}
		}
		tofile.close();
	}
	void ToNWchemFileHF(const string filename, const string basis = "6-31G")
	{
		ofstream out(filename.c_str(), ios::out);
		out << "# ================================================================" << endl;
		out << "# # NWChem input file made by LiTao for Molecule Recognization" << endl;
		out << "# ================================================================" << endl << endl;
		out << "charge 0" << endl << endl;
		out << "geometry" << endl;
		vector<string>::iterator iter;
		int i = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
		{
			out << " " << *iter;
			for (int j = 0; j != 3; j++)
			{
				if (abs(corr[i][j]) < 1e-7)
					corr[i][j] = 0;
				out << " " << fixed << setprecision(7) << corr[i][j];
			}
			out << endl;
		}
		out << "end" << endl << endl;
		out << "basis  \"ao basis\" spherical" << endl;
		out << " * library " << basis << endl;
		out << "end" << endl << "scf" << endl << " Singlet" << endl << "end" << endl << endl;
		out << "task SCF" << endl;
		out.close();
	}
	friend void ToNWchemFileHF(const Molecule &a, const Molecule &b, string &filename, const string basis = "6-31G")
	{
		ofstream out(filename.c_str(), ios::out);
		out << "# ================================================================" << endl;
		out << "# # NWChem input file made by LiTao for Molecule Recognization" << endl;
		out << "# ================================================================" << endl << endl;
		out << "charge 0" << endl << endl;
		out << "geometry" << endl;
		for (unsigned int i = 0; i != a.number; i++)
		{
			out << " " << a.name[i];
			for (int j = 0; j != 3; j++)
			{
				out << " " << fixed << setprecision(7) << a.corr[i][j];
			}
			out << endl;
		}
		for (unsigned int i = 0; i != b.number; i++)
		{
			out << " " << b.name[i];
			for (int j = 0; j != 3; j++)
			{
				out<< " " << fixed << setprecision(7) << b.corr[i][j];
			}
			out << endl;
		}
		out << "end" << endl << endl;
		out << "basis  \"ao basis\" spherical" << endl;
		out << " * library " << basis << endl;
		out << "end" << endl << "scf" << endl << " Singlet" << endl << "end" << endl << endl;
		out << "task SCF" << endl;
		out.close();
	}
	void ToG09FileDFT(string &filename, string basis = "6-31g", string functional="b3lyp")
	{
		ofstream out((filename + ".gjf").c_str(), ios::out);
		out << "\%nprocshared=7" << endl;
		out << "\%mem=8GB" << endl;
		out << "\%chk=system.chk" << endl;
		out << "#p "<<functional <<"/" << basis << endl;
		out << endl;
		out << "Title Card Required" << endl << endl;
		out << 0 << " " << 1 << endl;
		vector<string>::iterator iter;
		int i = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
		{
			out << " " << *iter;
			for (int j = 0; j != 3; j++)
			{
				if (abs(corr[i][j]) < 1e-7)
					corr[i][j] = 0;
				out << " " << fixed << setprecision(7) << corr[i][j];
			}
			out << endl;
		}
		out << endl;
		out << endl << endl << endl << endl;
		out.close();
	}
	friend void ToG09FileDFT(Molecule &a, Molecule &b, string &filename, string basis = "6-31g", string functional = "b3lyp")
	{
		ofstream out((filename + ".gjf").c_str(), ios::out);
		out << "\%nprocshared=7" << endl;
		out << "\%mem=8GB" << endl;
		out << "\%chk=system.chk" << endl;
		out << "#p "<<functional<<"/" << basis << endl;
		out << endl;
		out << "Title Card Required" << endl << endl;
		out << 0 << " " << 1 << endl;
		vector<string>::iterator iter;
		int i = 0;
		for (iter = a.name.begin(), i = 0; iter != a.name.end(); i++, iter++)
		{
			out << " " << a.name[i];
			for (int j = 0; j != 3; j++)
			{
				if (abs(a.corr[i][j]) < 1e-7)
					a.corr[i][j] = 0;
				out << " " << fixed << setprecision(7) << a.corr[i][j];
			}
			out << endl;
		}
		for (iter = b.name.begin(), i = 0; iter != b.name.end(); i++, iter++)
		{
			out << " " << b.name[i];
			for (int j = 0; j != 3; j++)
			{
				if (abs(b.corr[i][j]) < 1e-7)
					b.corr[i][j] = 0;
				out << " " << fixed << setprecision(7) << b.corr[i][j];
			}
			out << endl;
		}
		out << endl;
		out << endl << endl << endl << endl;
		out.close();
	}
	void ReadFromGJF(string &filename, int atomNum)
	{
		ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << " to get the geometry info" << endl;
			exit(1);
		}
		//Clear and initialize the info
		if (number != 0)
			clear();
		string temp;
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		getline(infile, temp);
		for (int i = 0; i != atomNum; i++)
		{
			string iname;
			infile >> iname;
			name.push_back(iname);
			for (int j = 0; j != 3; j++)
			{
				double icorr;
				infile >> icorr;
				corr[i][j] = icorr;
			}
		}
		number = name.size();
		infile.close();
	}

	void AddAtom(const string atomname, double &x, double &y, double &z)
	{
		corr[number][0] = x;
		corr[number][1] = y;
		corr[number][2] = z;
		number += 1;
		name.push_back(atomname);
	}
	//small functions
	int Number() { return number; }
	double MoleculeMass()
	{
		double totalMass = 0;
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			totalMass += AtomMass(*iter);
		return totalMass;
	}
	Eigen::Vector3d MassCenter()
	{
		double totalMass = 0;
		Eigen::Vector3d massCenter;
		double x[3];
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			totalMass += AtomMass(*iter);
		int i = 0;
		double multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += AtomMass(*iter)*corr[i][0];
		x[0] = multi / totalMass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += AtomMass(*iter)*corr[i][1];
		x[1] = multi / totalMass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += AtomMass(*iter)*corr[i][2];
		x[2] = multi / totalMass;
		massCenter << x[0], x[1], x[2];
		return massCenter;
	}
	friend double DistanceOfMassCenter(Molecule &ia,  Molecule &ib)
	{
		Eigen::Vector3d x1 = ia.MassCenter();
		Eigen::Vector3d x2 = ib.MassCenter();
		x1 = x1 - x2;
		return sqrt(x1(0)*x1(0) + x1(1)*x1(1) + x1(2)*x1(2));
	}
	
	void clear()
	{
		number = 0;
		name.clear();
		for (int i = 0; i != MaxAtom; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = -10242048;
	}
	//Simple Operate
	void PerformRot(Eigen::Matrix3d rot)
	{
		Eigen::Vector3d singleAtom;
		for (int i = 0; i != number; i++)
		{
			singleAtom << corr[i][0], corr[i][1], corr[i][2];
			singleAtom = rot*singleAtom;
			corr[i][0] = singleAtom(0);
			corr[i][1] = singleAtom(1);
			corr[i][2] = singleAtom(2);
		}
	}
	void PerformTrans(const Eigen::Vector3d trans)//This is to add a vector to molecular coordiantes
	{
		Eigen::Vector3d singleAtom;
		for (int i = 0; i != number; i++)
		{
			singleAtom << corr[i][0], corr[i][1], corr[i][2];
			singleAtom = trans + singleAtom;
			corr[i][0] = singleAtom(0);
			corr[i][1] = singleAtom(1);
			corr[i][2] = singleAtom(2);
		}
	}
	void PerformXTrans(const double &deltaX)//This is to make all atoms translate a deltaX in x coordinate.
	{
		for (int i = 0; i != number; i++)
			corr[i][0] += deltaX;
	}
	void PerformZTrans(const double &deltaZ)
	{
		for (int i = 0; i != number; i++)
			corr[i][2] += deltaZ;
	}
	void PerformAxisRot(Eigen::Vector3d axis, double angle_radian)//Rot with one axis
	{
		Eigen::AngleAxis<double> rot(angle_radian,axis);
		Eigen::Vector3d singleAtom;
		for (int i = 0; i != number; i++)
		{
			singleAtom << corr[i][0], corr[i][1], corr[i][2];
			singleAtom = rot*singleAtom;
			corr[i][0] = singleAtom(0);
			corr[i][1] = singleAtom(1);
			corr[i][2] = singleAtom(2);
		}
	}
	void PerformOnePointRotToXMinus(Eigen::Vector3d point)
	{
		double a = point(0), b = point(1), c = point(2);
		double alpha =acos( b / (sqrt(b*b + c*c)) );
		double beta = 1*acos( sqrt(b*b+c*c)/sqrt(a*a+b*b+c*c))+PI/2;
		Eigen::Vector3d e_x(1, 0, 0);
		Eigen::Vector3d e_z(0, 0, 1);
		Eigen::AngleAxis<double> rot1(alpha,e_x);
		Eigen::AngleAxis<double> rot2(beta, e_z);
		Eigen::Vector3d singleAtom;
		for (int i = 0; i != number; i++)
		{
			singleAtom << corr[i][0], corr[i][1], corr[i][2];
			singleAtom = rot2*(rot1*singleAtom);
			corr[i][0] = singleAtom(0);
			corr[i][1] = singleAtom(1);
			corr[i][2] = singleAtom(2);
		}
	}
	const Eigen::Matrix3d EulerRot(double &a1, double &a2, double &a3)
	{
		Eigen::Matrix3d m;
		m = Eigen::AngleAxisd(a1, Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(a2, Eigen::Vector3d::UnitY())* Eigen::AngleAxisd(a3, Eigen::Vector3d::UnitZ());
		return m;
	}
	void PerformRandomRotEuler(Eigen::Vector3d Point,const double &Precision_degree)
	{
		double x[3];
		const int nums = (int)(360 / Precision_degree);
		clock_t now = clock();
		srand(now);
		for (int i = 0; i != 3; i++)
			x[i] = (rand() % nums + 1)*Precision_degree;
		//std::default_random_engine generator(now);
		/*	std::default_random_engine generator(now);
		std::uniform_int_distribution<int> dis(0, nums);
		for(int i=0;i!=3;i++)
		{
		x[i] = dis(generator)*precision;
		}
		*/
		Eigen::Matrix3d randomRotation;
		randomRotation = EulerRot(x[0], x[1], x[2]);
		PerformTrans(-1 * Point);
		PerformRot(randomRotation);
		PerformTrans(Point);
	}
	void MCtoOrigin()
	{
		Eigen::Vector3d mc;
		mc = MassCenter();
		for (int i = 0; i != number; i++)
		{
			corr[i][0] -= mc(0);
			corr[i][1] -= mc(1);
			corr[i][2] -= mc(2);
		}
	}
	void MCtoVector(const Eigen::Vector3d x)
	{
		Eigen::Vector3d mc;
		mc = x - MassCenter();
		for (int i = 0; i != number; i++)
		{
			corr[i][0] += mc(0);
			corr[i][1] += mc(1);
			corr[i][2] += mc(2);
		}
	}
	//待完善

};
class JudgeMolecule:public Molecule
{
protected:
	double score;
public:
	JudgeMolecule() :Molecule() {  score = 0; }
	JudgeMolecule(JudgeMolecule &ia)
	{
		number = ia.number;
		name = ia.name;
		score = ia.score;
		for (int i = 0; i != number; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = ia.corr[i][j];
	}
	JudgeMolecule(Molecule &ia, double &ib) :Molecule(ia) { score = ib;}
	double Score() { return score; }
	JudgeMolecule& operator=(JudgeMolecule &ia)
	{
		number = ia.number;
		name = ia.name;
		score = ia.score;
		for (int i = 0; i != number; i++)
			for (int j = 0; j != 3; j++)
				corr[i][j] = ia.corr[i][j];
		return *this;
	}
};
class DoubleMolecule
{
private:
	Molecule a, b;
	double energy;
public:
	DoubleMolecule(Molecule &ia, Molecule &ib, double &ienergy)
	{
		a = ia;
		b = ib;
		energy = ienergy;
	}
	DoubleMolecule(const DoubleMolecule &i)
	{
		a = i.a;
		b = i.b;
		energy = i.energy;
	}
	DoubleMolecule()
	{
		Molecule ia, ib;
		a = ia;
		b = ib;
		energy = 0;
	}
	DoubleMolecule& operator=(DoubleMolecule &id)
	{
		a = id.a;
		b = id.b;
		energy = id.energy;
		return *this;
	}
	bool friend operator>(DoubleMolecule &ia, DoubleMolecule &ib)
	{
		if (ia.Energy() > ib.Energy())
			return true;
		else
			return false;
	}
	bool friend operator<(DoubleMolecule &ia, DoubleMolecule &ib)
	{
		if (ia.Energy() < ib.Energy())
			return true;
		else
			return false;
	}
	void Set(Molecule ia, Molecule ib, double ienergy)
	{
		a = ia;
		b = ib;
		energy = ienergy;
	}
	void GetInfo(Molecule &ia, Molecule &ib, double &ienergy)
	{
		ia = a;
		ib = b;
		ienergy = energy;
	}
	void ToXYZ(string filename, int count=0)
	{
		if(count==0)
			ToXYZfile(a, b, filename, X_ToStr<double>(energy)+" (energy)");
		else
			ToXYZfile(a, b, filename, X_ToStr<double>(energy) + " (energy),  " + X_ToStr<int>(count));

	}
	void output()
	{
		cout << "Energy: " << energy << endl;
		cout << a << endl;
		cout << b << endl;
	}
	double Energy()
	{
		return energy;
	}
};

class Fragments
{
protected:
	int frag_number;
	vector<Molecule> frags;
	int index;
public:
	Fragments() { frag_number = 0; index = 0; };
	Fragments(Fragments &ia)
	{
		frag_number = ia.frag_number;
		index = ia.index;
		for (int i = 0; i != ia.frag_number; i++)
			frags.push_back(ia.frags[i]);
	}
	friend ostream& operator<<(ostream &os, Fragments &ia)
	{
		os << ia.frag_number << " fragments" << endl;
		for (int i = 0; i != ia.frag_number; i++)
		{
			os <<ia.frags[i];
		}
		return os;
	}
	Fragments& operator++()
	{
		if (index < frag_number)
			index += 1;
		return *this;
	}
	void IndexToZero() { index = 0; }
	Molecule ThisFragment()
	{
		if (index < frag_number)
			return frags[index];
		else
			return frags[frag_number - 1];
	}
	Molecule OneFragment(unsigned int &my_index)
	{
		if (index < frag_number)
			return frags[my_index];
		else
			return frags[frag_number - 1];
	}
	Molecule OtherFragments()
	{
		Molecule temp;
		for (int i = 0; i != frag_number; i++)
		{
			if (i != index)
				temp = temp + frags[i];
		}
		return temp;
	}
	Molecule TotalFragments()
	{
		Molecule temp;
		for (int i = 0; i != frag_number; i++)
		{
				temp = temp + frags[i];
		}
		return temp;
	}
	Eigen::Vector3d VectorFromThisToOther()
	{
		return OtherFragments().MassCenter() - ThisFragment().MassCenter();
	}
	int FragNumbers()
	{
		return frag_number;
	}
	void ReadFromXYZfile(const string filename, const int inumber, int matrix[][2])
	{
		frag_number = inumber;
		index = 0;
		ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << " to get the geometry info" << endl;
			exit(1);
		}
		int number;
		infile >> number;
		string temp;
		getline(infile, temp);
		getline(infile, temp);
		string iname;
		double ix,iy,iz;
		Molecule temp_m;
		for (int i = 0; i != inumber; i++)
		{
			temp_m.clear();
			for (int j = matrix[i][0]; j <= matrix[i][1]; j++)
			{
				infile >> iname >> ix >> iy >> iz;
				temp_m.AddAtom(iname, ix, iy, iz);
			}
			frags.push_back(temp_m);
		}
		infile.close();
	}
};
void GenerateFunction();

/*

const int MAX_ATOM_NUMBER = 37;
const string ATOM[MAX_ATOM_NUMBER]= {"Null","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};
*/

/*
//已调试Atom类
class Atom
{
protected:
string name;
int number;
int charge;
public:
//构造函数
Atom()
{
name = ATOM[0];
number = 0;
charge = 0;
}
Atom(const string &a)
{
charge = 0;
int i = 0;
for (i = 0; i != MAX_ATOM_NUMBER; i++)
{
if (a == ATOM[i])
{
name = ATOM[i];
number = i;
break;
}
}
if (i == MAX_ATOM_NUMBER)
{
name = ATOM[0];
number = 0;
}
}
Atom(const string &a, int MyCharge)
{
charge = MyCharge;
int i = 0;
for (i = 0; i != MAX_ATOM_NUMBER; i++)
{
if (a == ATOM[i])
{
name = ATOM[i];
number = i;
break;
}
}
if (i == MAX_ATOM_NUMBER)
{
name = ATOM[0];
number = 0;
}
}
Atom(const Atom &a)
{
name = a.name;
number = a.number;
charge = a.charge;
}
//操作符重载
friend ostream &operator <<(ostream &os, const Atom &ia)
{
os << ia.name << "\t" << ia.number << "\t" << ia.charge;
return os;
}
friend istream &operator >>(istream &in,  Atom &ia)
{
int icharge = 0, inumber = 0;
in >> ia.name;
int i = 0;
for (i = 0; i != MAX_ATOM_NUMBER; i++)
{
if (ia.name == ATOM[i])
{
inumber = i;
break;
}
}
if (i == MAX_ATOM_NUMBER)
{
ia.name = ATOM[0];
inumber = 0;
}
ia.number = inumber;
ia.charge = icharge;
return in;
}
friend bool operator==(const Atom &ia, const Atom &ib)
{
if (ia.name == ib.name&&ia.charge == ib.charge)
return true;
else
return false;
}
friend bool operator!=(const Atom &ia, const Atom &ib)
{
if (ia.name != ib.name || ia.charge != ib.charge)
return true;
else
return false;
}
friend bool operator>(const Atom &ia,const Atom &ib)
{
if (ia.number > ib.number)
return true;
else
return false;
}
friend bool operator<(const Atom &ia, const Atom &ib)
{
if (ia.number < ib.number)
return true;
else
return false;
}
Atom& operator=(const Atom &ia)
{
name = ia.name;
number = ia.number;
charge = ia.charge;
return *this;
}
Atom& operator=(const string &a)
{
charge = 0;
int i = 0;
for (i = 0; i != MAX_ATOM_NUMBER; i++)
{
if (a == ATOM[i])
{
name = ATOM[i];
number = i;
}
}
if (i == MAX_ATOM_NUMBER)
{
name = ATOM[0];
number = 0;
}
return *this;
}
//简单函数
int Number() { return number; }
int Charge() { return charge; }
void ChangeCharge(int ic) { charge = ic; }
string Name() { return name; }
void clear()
{
name.clear();
number = 0;
charge = 0;
}

};

//已调试Point3D类
class Point3D
{
protected:
double corr[3];
public:
//构造函数
Point3D()
{
for (int i = 0; i != 3; i++)
corr[i] = -10242048;
}
Point3D(const Point3D &a)
{
for (int i = 0; i != 3; i++)
corr[i] = a.corr[i];
}
Point3D(const double &x, const double &y, const double &z)
{
corr[0] = x;
corr[1] = y;
corr[2] = z;
}
Point3D(Eigen::Vector3d &X)
{
corr[0] = X(0);
corr[1] = X(1);
corr[2] = X(2);
}
//操作符重载
bool operator==(const Point3D &a)
{
if (abs(corr[0]-a.corr[0])<10e-6 && abs(corr[1]- a.corr[1])<10e-6 && abs(corr[2]-a.corr[2])<10e-6)
return true;
else
return false;
}
void operator=(const Point3D &a)
{
for (int i = 0; i != 3;i++)
corr[i] = a.corr[i];
}
Point3D operator+(Eigen::Vector3d &X)
{
Point3D thispoint;
thispoint.corr[0] += X(0);
thispoint.corr[1] += X(1);
thispoint.corr[2] += X(2);
return thispoint;
}
Point3D operator-(Eigen::Vector3d &X)
{
Point3D thispoint;
thispoint.corr[0] -= X(0);
thispoint.corr[1] -= X(1);
thispoint.corr[2] -= X(2);
return thispoint;
}
friend Point3D MiddlePoint(Point3D &a, Point3D &b)
{
Point3D thispoint;
for (int i = 0; i != 3; i++)
thispoint.corr[i] = (a.corr[i]+b.corr[i])/2;
return thispoint;
}
void clear()
{
for (int i = 0; i != 3; i++)
corr[i] = -10242048;
}
friend ostream &operator<<(ostream &os, const Point3D &ip)
{
os << ip.corr[0] << "\t" << ip.corr[1] << "\t" << ip.corr[2];
return os;
}
friend istream &operator>>(istream &in, Point3D &ip)
{
in >> ip.corr[0]>>ip.corr[1]>>ip.corr[2];
return in;
}
//简单函数
};

class AtomMolecule:public Atom,public Point3D
{
public:
//构造函数
AtomMolecule(){}
AtomMolecule(Atom &ia, Point3D &ip):Atom(ia),Point3D(ip){}
AtomMolecule(string &ia,double &ix,double &iy,double &iz):Atom(ia),Point3D(ix,iy,iz){}
AtomMolecule(string &ia, int &ic, double &ix, double &iy, double &iz) :Atom(ia,ic), Point3D(ix, iy, iz) {}
AtomMolecule(string &ia, Eigen::Vector3d &ix):Atom(ia),Point3D(ix){}
AtomMolecule(AtomMolecule &ia)
{
name = ia.name;
number = ia.number;
charge = ia.charge;
for (int i = 0; i != 3; i++)
corr[i] = ia.corr[i];
}
//操作符重载
friend ostream& operator<<(ostream &os, const AtomMolecule &ia)
{
os << ia.name << "\t" << ia.corr[0] << "\t" << ia.corr[1] << "\t" << ia.corr[2];
return os;
}
friend istream &operator>>(istream &in, AtomMolecule &ia)
{
int icharge = 0, inumber = 0;
in >> ia.name;
int i = 0;
for (i = 0; i != MAX_ATOM_NUMBER; i++)
{
if (ia.name == ATOM[i])
{
inumber = i;
break;
}
}
if (i == MAX_ATOM_NUMBER)
{
ia.name = ATOM[0];
inumber = 0;
}
ia.number = inumber;
ia.charge = icharge;
in >> ia.corr[0] >> ia.corr[1] >> ia.corr[2];
return in;
}
bool operator==(AtomMolecule &ia)
{
if (name == ia.name && charge == ia.charge &&corr[0] == ia.corr[0] && corr[1] == ia.corr[1] && corr[2] == ia.corr[2])
return true;
else
return false;
}
void operator=(AtomMolecule &ia)
{
name = ia.name;
number = ia.number;
charge = ia.charge;
for (int i = 0; i != 3; i++)
corr[i] = ia.corr[i];
}
AtomMolecule &operator+(Eigen::Vector3d ix)
{
;
}
void clear()
{
Atom::clear();
Point3D::clear();
}

//简单函数

};
*/