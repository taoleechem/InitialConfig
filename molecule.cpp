#include "molecule.h"
#include "rmsd.h" 

Molecule::Molecule()
{
	number = 0;
	for (int i = 0; i != MaxAtom; i++)
		for (int j = 0; j != 3; j++)
			corr[i][j] = -10242048;
}
Molecule::Molecule(const Molecule &ia)
{
	number = ia.number;
	name = ia.name;
	for (int i = 0; i != number; i++)
		for (int j = 0; j != 3; j++)
			corr[i][j] = ia.corr[i][j];
}
Molecule& Molecule::operator=(const Molecule &ia)
{
	number = ia.number;
	name = ia.name;
	for (int i = 0; i != number; i++)
		for (int j = 0; j != 3; j++)
			corr[i][j] = ia.corr[i][j];
	return *this;
}
Molecule operator+(const Molecule &ia, const Molecule &ib)
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
			temp.corr[i + ia.number][j] = ib.corr[i][j];
	}
	return temp;
}
bool operator==(const Molecule &ia, const Molecule &ib)
{
	if (ia.number == ib.number&&ia.name == ib.name)
		return true;
	else return false;
}
bool operator!=(const Molecule &ia, const Molecule &ib)
{
	if (ia.number != ib.number || ia.name != ib.name)
		return true;
	else return false;
}
bool operator>=(const Molecule &ia, const Molecule &ib)
{
	if (ia.number >= ib.number)
		return true;
	else
		return false;
}
double Molecule::AtomMass(const string &iname)
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
string Molecule::MoleculeName()
{
	vector<string> new_name(name);
	sort(new_name.begin(), new_name.end());
	vector<string>::iterator end_unique = unique(new_name.begin(), new_name.end());
	new_name.erase(end_unique, new_name.end());
	vector<string>::iterator iter1, iter2;
	int count = 0;
	string atom_name;
	for (iter1 = new_name.begin(); iter1 != new_name.end(); iter1++)
	{
		atom_name += *iter1;
		for (iter2 = name.begin(); iter2 != name.end(); iter2++)
		{
			if (*iter2 == *iter1)
				count += 1;
		}
		if (count >= 2)
			atom_name += X_ToStr<int>(count);
		count = 0;
	}
	return atom_name;
}
ostream& operator<<(ostream &os, Molecule &ia)
{
	os << ia.number << endl;
	os << "standard .xyz file  " << ia.MoleculeName() << endl;
	for (int i = 0; i != ia.number; i++)
	{
		os << ia.name[i];
		for (int j = 0; j != 3; j++)
		{
			if (abs(ia.corr[i][j]) < 1e-7)
				ia.corr[i][j] = 0;
			os << "\t" << fixed << setprecision(7) << ia.corr[i][j];
		}
		os << endl;
	}
	return os;
}
void Molecule::ReadFromXYZfile(const string filename)
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
void Molecule::ToXYZfile(const string filename)
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	if (number == 0)
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
void ToXYZfile(const Molecule &a, const Molecule &b, string &filename, string other_info)
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
void Molecule::ToNWchemFileHF(const string filename, const string basis)
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
void ToNWchemFileHF(const Molecule &a, const Molecule &b, string &filename, const string basis)
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
			out << " " << fixed << setprecision(7) << b.corr[i][j];
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
void Molecule::ToG09FileDFT(string &filename, string basis , string functional)
{
	ofstream out((filename + ".gjf").c_str(), ios::out);
	out << "\%nprocshared=7" << endl;
	out << "\%mem=8GB" << endl;
	out << "\%chk=system.chk" << endl;
	out << "#p " << functional << "/" << basis << endl;
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
void ToG09FileDFT(Molecule &a, Molecule &b, string &filename, string basis, string functional )
{
	ofstream out((filename + ".gjf").c_str(), ios::out);
	out << "\%nprocshared=7" << endl;
	out << "\%mem=8GB" << endl;
	out << "\%chk=system.chk" << endl;
	out << "#p " << functional << "/" << basis << endl;
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
void Molecule::ReadFromGJF(string &filename, int atomNum)
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
void Molecule::AddAtom(const string atomname, double &x, double &y, double &z)
{
	corr[number][0] = x;
	corr[number][1] = y;
	corr[number][2] = z;
	number += 1;
	name.push_back(atomname);
}
//small functions
int Molecule::Number() { return number; }
double Molecule::MoleculeMass()
{
	double totalMass = 0;
	vector<string>::iterator iter;
	for (iter = name.begin(); iter != name.end(); iter++)
		totalMass += AtomMass(*iter);
	return totalMass;
}
Eigen::Vector3d Molecule::MassCenter()
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
double DistanceOfMassCenter(Molecule &ia, Molecule &ib)
{
	Eigen::Vector3d x1 = ia.MassCenter();
	Eigen::Vector3d x2 = ib.MassCenter();
	x1 = x1 - x2;
	return sqrt(x1(0)*x1(0) + x1(1)*x1(1) + x1(2)*x1(2));
}

void Molecule::clear()
{
	number = 0;
	name.clear();
	for (int i = 0; i != MaxAtom; i++)
		for (int j = 0; j != 3; j++)
			corr[i][j] = -10242048;
}
//Simple Operate
void Molecule::PerformRot(Eigen::Matrix3d rot)
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
void Molecule::PerformTrans(const Eigen::Vector3d trans)//This is to add a vector to molecular coordiantes
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
void Molecule::PerformXTrans(const double &deltaX)//This is to make all atoms translate a deltaX in x coordinate.
{
	for (int i = 0; i != number; i++)
		corr[i][0] += deltaX;
}
void Molecule::PerformZTrans(const double &deltaZ)
{
	for (int i = 0; i != number; i++)
		corr[i][2] += deltaZ;
}
void Molecule::PerformAxisRot(Eigen::Vector3d axis, double angle_radian)//Rot with one axis
{
	Eigen::AngleAxis<double> rot(angle_radian, axis);
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
void Molecule::PerformOnePointRotToXMinus(Eigen::Vector3d point)
{
	double a = point(0), b = point(1), c = point(2);
	double alpha = acos(b / (sqrt(b*b + c*c)));
	double beta = 1 * acos(sqrt(b*b + c*c) / sqrt(a*a + b*b + c*c)) + PI / 2;
	Eigen::Vector3d e_x(1, 0, 0);
	Eigen::Vector3d e_z(0, 0, 1);
	Eigen::AngleAxis<double> rot1(alpha, e_x);
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
Eigen::Matrix3d Molecule::EulerRot(double &a1, double &a2, double &a3)
{
	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(a1, Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(a2, Eigen::Vector3d::UnitY())* Eigen::AngleAxisd(a3, Eigen::Vector3d::UnitZ());
	return m;
}
void Molecule::PerformRandomRotEuler(Eigen::Vector3d Point, const double &Precision_degree)
{
	double x[3];
	const int nums = (int)(360 / Precision_degree);
	
	clock_t now = clock();
	/*
	srand(now);
	for (int i = 0; i != 3; i++)
	x[i] = (rand() % nums + 1)*Precision_degree;
	*/

	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(0, nums);
	for(int i=0;i!=3;i++)
	{
	x[i] = dis(generator)*Precision_degree;
	}
	
	Eigen::Matrix3d randomRotation;
	randomRotation = EulerRot(x[0], x[1], x[2]);
	PerformTrans(-1 * Point);
	PerformRot(randomRotation);
	PerformTrans(Point);
}
void Molecule::MCtoOrigin()
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
void Molecule::MCtoVector(const Eigen::Vector3d x)
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
void MakeAtomsSuitableDistanceMoveB(Molecule &ia, Molecule &ib, const double SmallestDistance)
{
	int label1 = 0, label2 = 0;
	double distance = sqrt((ia.corr[0][0] - ib.corr[0][0])*(ia.corr[0][0] - ib.corr[0][0]) + (ia.corr[0][1] - ib.corr[0][1])*(ia.corr[0][1] - ib.corr[0][1]) + (ia.corr[0][2] - ib.corr[0][2])*(ia.corr[0][2] - ib.corr[0][2]));
	for (int i = 0; i != ia.number; i++)
		for (int j = 0; j != ib.number; j++)
		{
			double temp = sqrt((ia.corr[i][0] - ib.corr[j][0])*(ia.corr[i][0] - ib.corr[j][0]) + (ia.corr[i][1] - ib.corr[j][1])*(ia.corr[i][1] - ib.corr[j][1]) + (ia.corr[i][2] - ib.corr[j][2])*(ia.corr[i][2] - ib.corr[j][2]));
			if (temp < distance)
			{
				distance = temp;
				label1 = i;
				label2 = j;
			}
		}
	Eigen::Vector3d MoveVector;
	MoveVector << ia.corr[label1][0] - ib.corr[label2][0], ia.corr[label1][1] - ib.corr[label2][1], ia.corr[label1][2] - ib.corr[label2][2];
	if (abs(distance - SmallestDistance)>1e-6)
		ib.PerformTrans((1 - SmallestDistance / distance)*MoveVector);
}
double RMSD(Molecule &ia, Molecule &ib)
{
	if (ia.number != ib.number)
		return 1;
	else
	{
		double ref_xlist[MaxAtom][3];
		double mov_xlist[MaxAtom][3];
		int n_list=ia.number;
		double rmsd=0;
	//Initilization of 2 list

	for (int i = 0; i != ia.number; i++)
	for (int j = 0; j != 3; j++)
	{
	ref_xlist[i][j] = ia.corr[i][j];
	mov_xlist[i][j] = ib.corr[i][j];
	}

	fast_rmsd(ref_xlist, mov_xlist, n_list, rmsd);
	return rmsd;
	}
}




DoubleMolecule::DoubleMolecule(Molecule &ia, Molecule &ib, double &ienergy)
{
	a = ia;
	b = ib;
	energy = ienergy;
}
DoubleMolecule::DoubleMolecule(const DoubleMolecule &i)
{
	a = i.a;
	b = i.b;
	energy = i.energy;
}
DoubleMolecule::DoubleMolecule()
{
	Molecule ia, ib;
	a = ia;
	b = ib;
	energy = 0;
}
DoubleMolecule& DoubleMolecule::operator=(DoubleMolecule &id)
{
	a = id.a;
	b = id.b;
	energy = id.energy;
	return *this;
}
bool operator>(DoubleMolecule &ia, DoubleMolecule &ib)
{
	if (ia.Energy() > ib.Energy())
		return true;
	else
		return false;
}
bool operator<(DoubleMolecule &ia, DoubleMolecule &ib)
{
	if (ia.Energy() < ib.Energy())
		return true;
	else
		return false;
}
void DoubleMolecule::Set(Molecule ia, Molecule ib, double ienergy)
{
	a = ia;
	b = ib;
	energy = ienergy;
}
void DoubleMolecule::GetInfo(Molecule &ia, Molecule &ib, double &ienergy)
{
	ia = a;
	ib = b;
	ienergy = energy;
}
void DoubleMolecule::ToXYZ(string filename)
{
		ToXYZfile(a, b, filename, X_ToStr<double>(energy) + " (energy)");
}
void DoubleMolecule::output()
{
	cout << "Energy: " << energy << endl;
	cout << a << endl;
	cout << b << endl;
}
double DoubleMolecule::Energy()
{
	return energy;
}
double RMSD(DoubleMolecule &ia, DoubleMolecule &ib)
{
	Molecule t1, t2;
	t1 = ia.a + ia.b;
	t2 = ib.a + ib.b;
	return RMSD(t1, t2);
}


Fragments::Fragments() { frag_number = 0; index = 0; };
Fragments::Fragments(Fragments &ia)
{
	frag_number = ia.frag_number;
	index = ia.index;
	for (int i = 0; i != ia.frag_number; i++)
		frags.push_back(ia.frags[i]);
}
ostream& operator<<(ostream &os, Fragments &ia)
{
	os << ia.frag_number << " fragments" << endl;
	for (int i = 0; i != ia.frag_number; i++)
	{
		os << ia.frags[i];
	}
	return os;
}
Fragments& Fragments::operator++()
{
	if (index < frag_number)
		index += 1;
	return *this;
}
void Fragments::IndexToZero() { index = 0; }
Molecule Fragments::ThisFragment()
{
	if (index < frag_number)
		return frags[index];
	else
		return frags[frag_number - 1];
}
Molecule Fragments::OneFragment(unsigned int &my_index)
{
	if (index < frag_number)
		return frags[my_index];
	else
		return frags[frag_number - 1];
}
Molecule Fragments::OtherFragments()
{
	Molecule temp;
	for (int i = 0; i != frag_number; i++)
	{
		if (i != index)
			temp = temp + frags[i];
	}
	return temp;
}
Molecule Fragments::TotalFragments()
{
	Molecule temp;
	for (int i = 0; i != frag_number; i++)
	{
		temp = temp + frags[i];
	}
	return temp;
}
Eigen::Vector3d Fragments::VectorFromThisToOther()
{
	return OtherFragments().MassCenter() - ThisFragment().MassCenter();
}
int Fragments::FragNumbers()
{
	return frag_number;
}
void Fragments::ReadFromXYZfile(const string filename, const int inumber, int matrix[][2])
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
	double ix, iy, iz;
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
//Geometry operation
void Fragments::PerformRot(Eigen::Matrix3d rot)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformRot(rot);
}
void Fragments::PerformTrans(const Eigen::Vector3d trans)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformTrans(trans);
}
void Fragments::PerformXTrans(const double &deltaX)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformXTrans(deltaX);
}
void Fragments::PerformZTrans(const double &deltaZ)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformZTrans(deltaZ);
}
void Fragments::PerformAxisRot(Eigen::Vector3d axis, double angle_radian)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformAxisRot(axis, angle_radian);
}
void Fragments::PerformOnePointRotToXMinus(Eigen::Vector3d point)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformOnePointRotToXMinus(point);
}
