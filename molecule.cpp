/*
 * =====================================================================================
 *
 *       Filename:  molecule.cpp
 *
 *    Description:  Implementation of molecule.h
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

#include "molecule.h"
#include "rmsd.h" 
#include <iostream>

//Generate <double> random number from 0-MaxValue
double RandomNumber(double MaxValue)
{
	clock_t now = clock();
	/*
	srand(now);
	for (int i = 0; i != 3; i++)
	x[i] = (rand() % nums + 1)*Precision_degree;
	*/
	std::default_random_engine generator(now);
	std::uniform_real_distribution<double> dis(0, MaxValue);
	return dis(generator);

}
int RandomNumberInt(int MaxValue)
{
	clock_t now = clock();
	/*
	srand(now);
	for (int i = 0; i != 3; i++)
	x[i] = (rand() % nums + 1)*Precision_degree;
	*/
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(0, MaxValue);
	return dis(generator);

}

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

void Molecule::Get2DArray(double position[][3], int top_n)
{
	if (top_n > number)
		top_n = number;
	for (int i = 0; i < top_n; i++)
		for (int j = 0; j < 3; j++)
			position[i][j] = corr[i][j];
}

double Molecule::AtomMass(const string &name)
{
	if (name == "H")
		return 1.007825;
	else if (name=="B")
		return 10.811;
	else if (name == "C")
		return 12.00;
	else if (name == "N")
		return 14.003070;
	else if (name == "O")
		return 15.994910;
	else if (name=="F")
		return 18.998403;
	else if (name=="P")
		return 30.973761;
	else if (name == "S")
		return 31.972070;
	else if (name == "Cl")
		return 34.968850;
	else if (name=="Br")
		return 79.904;
	else if (name=="I")
		return 126.90447;
	else if (name =="Fe")
		return 55.845;
	else if (name=="Al")
		return 26.981538;
	else if (name=="Mg")
		return 24;
	//This part is to deal with the non-standard atom name
	else if ((name.size() == 2) && (name[1] <= 'Z'&&name[1] >= 'A'))
	{
		string iname = X_ToStr<char>(name[0]);
		return AtomMass(iname);
	}
	else
		return 0;
	
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
	clear();
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
void Molecule::ReadFromXYZOnlyGeo(ifstream &infile, int atom_num)
{
	clear();
	number = atom_num;
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
}
//change tinker xyz to xyz atomatically
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
			if (name[i].size() > 1)
			{
				if (name[i][1] <= 'Z'&&name[i][1] >= 'A')
					tofile << name[i][0];
				else
					tofile << name[i];
			}
			else
				tofile << name[i];
			for (int j = 0; j != 3; j++)
			{
				if (abs(corr[i][j]) < 1e-7)
					corr[i][j] = 0;
				tofile << " " << fixed << setprecision(7) << corr[i][j];
			}
			tofile << endl;
		}
	}
	tofile.close();
}
//bool judge=true means output the standard atom name(C), false means output the original name(CT)
void Molecule::ToXYZfileOnlyGeo(ofstream &tofile, bool IfStandardAtomName)
{
		for (int i = 0; i != number; i++)
		{

			if (IfStandardAtomName == true && name[i].size() > 1)
			{
				if (name[i][1] <= 'Z'&&name[i][1] >= 'A')
					tofile << name[i][0];
				else
					tofile << name[i];
			}
			else
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
void Molecule::ToMol2File(const string filename)
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	if (number == 0)
		cout << "Empty molecule and no info is written to " << filename << endl;
	//scan and save bond connection info
	vector<Bond> bonds;
	for (int i = 0; i < number; i++)
		for (int j = i + 1; j < number; j++)
		{
			int order = BondOrder(i, j);
			Bond bond_a;
			if (order != 0)
			{
				bond_a.SetValue(i+1, j+1, order);
				bonds.push_back(bond_a);
			}
		}
	int bond_number = bonds.size();
	//output *.mol2 file;
	tofile << "@<TRIPOS>MOLECULE" << endl;
	tofile << "Molecule Name" << endl;
	tofile << number<<" " <<bond_number<< endl;
	tofile << "SMALL" << endl; 
	tofile << "NO_CHARGES" << endl<<endl;
	tofile << "@<TRIPOS>ATOM" << endl;
	for (int i = 0; i < number; i++)
	{
		tofile << i + 1 << " " <<StandardAtomName(name[i]) << i + 1 << "\t" <<fixed<<setprecision(4)<< corr[i][0] << "\t" << corr[i][1] << "\t" << corr[i][2] << " " << StandardAtomName(name[i]) << endl;
	}
	tofile << "@<TRIPOS>BOND" << endl;
	for (int i = 0; i < bond_number; i++)
	{
		tofile << i + 1 << " " << bonds[i] << endl;
	}
	tofile.close();
}
void Molecule::ToMol2File_AmBn(const string filename, int a_atoms, int a_numbers,  int b_atoms)
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	if (number == 0)
		cout << "Empty molecule and no info is written to " << filename << endl;
	int b_numbers = (number - a_numbers*a_atoms) / b_atoms;
	//find bond info of A
	vector<Bond> Abonds;
	for (int i = 0; i < a_atoms; i++)
		for (int j = i + 1; j < a_atoms; j++)
		{
			int order = BondOrder(i, j);
			Bond bond_a;
			if (order != 0)
			{
				bond_a.SetValue(i + 1, j + 1, order);
				//cout << bond_a << endl;
				Abonds.push_back(bond_a);
			}
		}
	vector<Bond> Bbonds;
	for (int i = a_numbers*a_atoms; i < a_numbers*a_atoms+b_atoms; i++)
		for (int j = i + 1; j < a_numbers*a_atoms + b_atoms; j++)
		{
			int order = BondOrder(i, j);
			Bond bond_b;
			if (order != 0)
			{
				bond_b.SetValue(i + 1, j + 1, order);
				//cout << bond_b << endl;
				Bbonds.push_back(bond_b);
			}
		}
	int bond_number = Abonds.size()*a_numbers+Bbonds.size()*b_numbers;
	//output *.mol2 file;
	tofile << "@<TRIPOS>MOLECULE" << endl;
	tofile << "Molecule Name" << endl;
	tofile << number << " " << bond_number << endl;
	tofile << "SMALL" << endl;
	tofile << "NO_CHARGES" << endl << endl;
	tofile << "@<TRIPOS>ATOM" << endl;
	for (int i = 0; i < number; i++)
	{
		tofile << i + 1 << " " << StandardAtomName(name[i]) << i + 1 << "\t" << fixed << setprecision(4) << corr[i][0] << "\t" << corr[i][1] << "\t" << corr[i][2] << " " << StandardAtomName(name[i]) << endl;
	}
	tofile << "@<TRIPOS>BOND" << endl;
	int i_bond = 0;
	for (int i = 0; i < a_numbers; i++)
	{
		for (int j = 0; j < Abonds.size(); j++)
		{
			tofile << i_bond + 1<<" ";
			i_bond += 1;
			tofile << Abonds[j] << endl;
			Abonds[j].IncreaseLabel(a_atoms);
		}
	}
	for (int i = 0; i < b_numbers; i++)
	{
		for (int j = 0; j < Bbonds.size(); j++)
		{
			tofile << i_bond + 1<<" ";
			i_bond += 1;
			tofile << Bbonds[j] << endl;
			Bbonds[j].IncreaseLabel(b_atoms);
		}
	}
	tofile.close();
}
void Molecule::ToMol2fileOnlyGeo(ofstream &tofile, int initial_label)
{
	for (int i = 0; i < number; i++)
	{
		tofile << i + initial_label << " " << StandardAtomName(name[i]) << i + initial_label << "\t" << fixed << setprecision(4) << corr[i][0] << "\t" << corr[i][1] << "\t" << corr[i][2] << " " << StandardAtomName(name[i]) << endl;
	}
}
void ToMol2fileAmBnType(const string filename, vector<Molecule> &As, vector<Molecule>& Bs)
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	vector<Bond> all_bonds;
	int a_number = As.size();
	int b_number = Bs.size();
	int a_atoms = As[0].number;
	int b_atoms = Bs[0].number;
	for (int i_molecule = 0; i_molecule < a_number; i_molecule++)
	{
		for (int i = 0; i < a_atoms; i++)
			for (int j = i + 1; j < a_atoms; j++)
			{
				int order = As[i_molecule].BondOrder(i, j);
				Bond temp;
				if (order != 0)
				{
					temp.SetValue(i +a_atoms*i_molecule + 1, j+a_atoms*i_molecule + 1, order);
					all_bonds.push_back(temp);
				}
			}
	}
	for (int i_molecule = 0; i_molecule < b_number; i_molecule++)
	{
		for (int i = 0; i < b_atoms; i++)
			for (int j = i + 1; j < b_atoms; j++)
			{
				int order = Bs[i_molecule].BondOrder(i, j);
				Bond temp;
				if (order != 0)
				{
					temp.SetValue(i+1+ b_atoms*i_molecule+a_atoms*a_number, j+1+ b_atoms*i_molecule + a_atoms*a_number, order);
					all_bonds.push_back(temp);
				}
			}
	}
	int bond_number = all_bonds.size();
	//output *.mol2 file;
	tofile << "@<TRIPOS>MOLECULE" << endl;
	tofile << "Molecule Name" << endl;
	tofile << a_number*a_atoms+b_atoms*b_number << " " << bond_number << endl;
	tofile << "SMALL" << endl;
	tofile << "NO_CHARGES" << endl << endl;
	tofile << "@<TRIPOS>ATOM" << endl;
	int i_count = 1;
	for (int i_molecule = 0; i_molecule < a_number; i_molecule++)
	{
		As[i_molecule].ToMol2fileOnlyGeo(tofile, i_count);
		i_count += a_atoms;
	}
	for (int i_molecule = 0; i_molecule < b_number; i_molecule++)
	{
		Bs[i_molecule].ToMol2fileOnlyGeo(tofile, i_count);
		i_count += b_atoms;
	}
	tofile << "@<TRIPOS>BOND" << endl;
	for (int i = 0; i != bond_number; i++)
		tofile << i + 1 << " " << all_bonds[i] << endl;
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
void Molecule::ToNWchemFileDFT(const string filename, const string basis, const string functional)
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
	out << "basis" << endl;
	out << " * library " << basis << endl;
	out << "end" << endl << endl;
	out << "dft" << endl;
	out << "  xc " << functional << endl;
	out << "  mult 1" << endl;
	out << "end" << endl << endl;
	out << "task dft energy" << endl;
	out.close();
}
void ToNWchemFileDFT(const Molecule &a, const Molecule &b, string &filename, const string basis, const string functional)
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
	out << "basis" << endl;
	out << " * library " << basis << endl;
	out << "end" << endl << endl;
	out << "dft" << endl;
	out << "  xc " << functional << endl;
	out << "  mult 1" << endl;
	out << "end" << endl << endl;
	out << "task dft energy" << endl;
	out.close();
}
void Molecule::ToG09FileDFT(string &filename, string basis, string functional, string othercommand)
{
	ofstream out((filename + ".gjf").c_str(), ios::out);
	out << "\%nprocshared=7" << endl;
	out << "\%mem=8GB" << endl;
	out << "\%chk=system.chk" << endl;
	out << "# " << functional << " " << basis <<" "<<othercommand<< endl;
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
void ToG09FileDFT(Molecule &a, Molecule &b, string &filename, string basis, string functional, string othercommand)
{
	ofstream out((filename + ".gjf").c_str(), ios::out);
	out << "\%nprocshared=7" << endl;
	out << "\%mem=8GB" << endl;
	out << "\%chk=system.chk" << endl;
	out << "# " << functional << " " << basis <<" "<<othercommand<< endl;
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
void Molecule::ReadFromTinkerXYZfile(const string filename)
{
	ifstream infile(filename.c_str());
	if (!std::cout)
	{
		cerr << "Error to open " << filename << " to get the geometry info" << endl;
		exit(1);
	}
	infile >> number;
	//Clear and initialize the info
	if (number != 0)
		name.clear();
	string temp;
	getline(infile, temp);
	int inumber;
	string iname, iconnection;
	for (int i = 0; i != number; i++)
	{
		infile >> inumber >> iname;
		name.push_back(iname);
		for (int j = 0; j != 3; j++)
		{
			double icorr;
			infile >> icorr;
			corr[i][j] = icorr;
		}
		getline(infile, iconnection);
	}
	infile.close();
}
void Molecule::ReadFromTinkerXYZGeoPart(ifstream &infile, int atom_number)
{
		number = atom_number;
		//Clear and initialize the info
		if (number != 0)
			name.clear();
		int inumber;
		string iname, iconnection;
		for (int i = 0; i != number; i++)
		{
			infile >> inumber >> iname;
			name.push_back(iname);
			for (int j = 0; j != 3; j++)
			{
				double icorr;
				infile >> icorr;
				corr[i][j] = icorr;
			}
			getline(infile, iconnection);
	      }
}
void Molecule::ToPDBfile(const string filename, int connection_matrix[][8])
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	tofile << "COMPND    " << MoleculeName() << endl;
	tofile << "AUTHOR    GENERATED BY " << "Tao Li" << endl;
	for (int i = 0; i != number; i++)
	{
		tofile << "HETATM    " << i+1 << "	 " << name[i] << "   LIG     1      " <<fixed<<setprecision(3)<< corr[i][0] << "  " << corr[i][1] << "  " << corr[i][2] << "  1.00  0.00           " << name[i] << endl;
	}
	for (int i = 0; i != number; i++)
	{
		tofile << "CONECT    "<<i+1;
		for (int j = 0; j != 8; j++)
		{
			if (connection_matrix[i][j] != 0)
				tofile << "    " << connection_matrix[i][j];
		}
		tofile << endl;
	}
	tofile << "END" << endl;
	tofile.close();
}
//initial_label comes from 0
void Molecule::ToPDBfileOnlyGeo(ofstream &tofile, int initial_label)
{
	for (int i = 0; i != number; i++)
	{
		tofile << "HETATM    " << initial_label+i+1 << "  " << name[i] << "   LIG     1      "  <<fixed<<setprecision(3)<< corr[i][0] << "  " << corr[i][1] << "  "  << corr[i][2] << "  1.00  0.00          " << name[i] << endl;
		cout << "HETATM    " << initial_label + i + 1<<"  " << name[i] << "   LIG     1      "  << fixed<< setprecision(3) <<corr[i][0] << "  "  << corr[i][1] << "  "  << corr[i][2] << "  1.00  0.00          " << name[i] << endl;
	}
}
void ToPDBfileAmBnType(const string filename, vector<Molecule> &A, vector<Molecule> &B, int connect_m[][8], int connect_n[][8])
{
	ofstream tofile(filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << filename << endl;
		exit(1);
	}
	int m = A.size();
	int n = B.size();
	int a_atoms = A[0].Number();
	int b_atoms = B[0].Number();
	tofile << "COMPND    AnalysisResult" << endl;
	cout << "COMPND    AnalysisResult"  << endl;
	tofile << "AUTHOR    GENERATED BY " << "Tao Li" << endl;
	cout << "AUTHOR    GENERATED BY " << "Tao Li" << endl;

	for (int i = 0; i != m; i++)
		A[i].ToPDBfileOnlyGeo(tofile, i*a_atoms);
	for (int i = 0; i != n; i++)
		B[i].ToPDBfileOnlyGeo(tofile, m*a_atoms +(i+m-1)*b_atoms);

	for (int i = 0; i != m; i++)
	{
		for (int j = 0; j != a_atoms; j++)
		{
			tofile << "CONECT    " << i*a_atoms + j + 1;
			cout << "CONECT    " << i*a_atoms + j + 1;
			for (int k = 0; k != 8; k++)
				if (connect_m[j][k] != 0)
				{
					tofile << "    " << connect_m[j][k] + i*a_atoms;
					cout << "    " << connect_m[j][k] + i*a_atoms;
				}
			tofile << endl;
			cout << endl;
		}
			
	}
	for (int i = m; i != m+n; i++)
	{
		for (int j = 0; j != B[0].number; j++)
		{
			tofile << "CONECT    " << m*a_atoms +(i-m)*b_atoms + j + 1;
			cout << "CONECT    " << m*a_atoms + (i - m)*b_atoms + j + 1;
			for (int k = 0; k != 8; k++)
				if (connect_n[j][k] != 0)
				{
					tofile << "    " << connect_n[j][k] + m*a_atoms + (i - m)*b_atoms;
					cout << "    " << connect_n[j][k] + m*a_atoms + (i - m)*b_atoms;
				}
			tofile << endl;
			cout << endl;
		}

	}
	tofile << "MASTER        0    0    0    0    0    0    0    0   " << a_atoms*m + b_atoms*n << "    0   " << a_atoms*m + b_atoms*n << "    0" << endl;
	cout << "MASTER        0    0    0    0    0    0    0    0   " << a_atoms*m + b_atoms*n << "    0   " << a_atoms*m + b_atoms*n << "    0" << endl;
	tofile << "END" << endl;
	cout << "END" << endl;
	tofile.close();
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
//This is to return radius of No.i atom
double Molecule::AtomRadius(int i)
{
	int j;
	for (j = 0; j != MaxElementOrder; j++)
	{
		if (name[i] == ELEMENT[j])
			break;
	}
	if (j < MaxElementOrder)
		return ELEMENT_RADIUS[j];
	else
		return 0.70;
}
//This is to return bond order of No.i atom & No.j atom
int Molecule::BondOrder(int i, int j)
{
	double distance = sqrt((corr[i][0] - corr[j][0])*(corr[i][0] - corr[j][0]) + (corr[i][1] - corr[j][1])*(corr[i][1] - corr[j][1])+ (corr[i][2] - corr[j][2])*(corr[i][2] - corr[j][2]));
	double sum = AtomRadius(i) + AtomRadius(j);
	double order = distance / sum;
	if (order > 0.62&&order < 0.72)
		return 4;
	else if (order > 0.72&&order < 0.82)
		return 3;
	else if (order > 0.82&&order < 0.92)
		return 2;
	else if (order > 0.92&&order < 1.25)
		return 1;
	else
		return 0;
}
int Molecule::ConnectAtoms(int i)
{
	int temp = 0;
	for (int j = 0; j != number; j++)
	{
		if (j == i)
			continue;
		else
		{
			if (BondOrder(i, j) > 0)
				temp += 1;
		}
	}
	return temp;
}

vector<int> Molecule::WhoConnectMe(int i, int except_label)
{
	vector<int> temp;
	for (int j = 0; j < number; j++)
	{
		if (j != except_label&&BondOrder(i, j) > 0)
		{
			temp.push_back(j);
			//cout << j << " connect " << i << " besides " << except_label << endl;
		}
	}
	return temp;
}

string Molecule::StandardAtomName(string x)
{
	string temp;
	if (x.size() == 2 && x[1] >= 'A'&&x[1] <= 'Z')
		temp += x[0];
	else
		temp = x;
	return temp;
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
//This is to form a new Molecule by abstract information of molecule ia atom[begin_label]-->atom[end_label], begin_label begins from 0
Molecule AbstractSomeFormNewMolecule(Molecule &ia, int begin_label, int end_label)
{
	Molecule temp;
	if (begin_label < 0 || end_label >= ia.Number())
		return temp;
	else
	{
		temp.number = end_label - begin_label + 1;
		temp.name.clear();
		for (int i = 0; i != temp.number; i++)
		{
			temp.name.push_back(ia.name[i + begin_label]);
			for (int j = 0; j != 3; j++)
				temp.corr[i][j] = ia.corr[i + begin_label][j];
		}
		return temp;
	}
}
//label is 1,2,...,last one
Molecule Molecule::NewMoleculeFromThis(int nums, int label[])
{
	Molecule temp;
	temp.number = nums;
	for (int i = 0; i < nums; i++)
	{
		temp.name.push_back(name[label[i] - 1]);
		for (int j = 0; j < 3;j++)
			temp.corr[i][j] = corr[label[i] - 1][j];
	}
	return temp;
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
void Molecule::PerformYTrans(const double &deltaY)//This is to make all atoms translate a deltaY in y coordinate.
{
	for (int i = 0; i != number; i++)
		corr[i][1] += deltaY;
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
	double alpha = atan2(c, b);
	//double beta = 1 * acos(sqrt(b*b + c*c) / sqrt(a*a + b*b + c*c)) + PI / 2;
	double beta = PI - atan2(sqrt(b*b + c*c), a);
	Eigen::Vector3d e_x(-1, 0, 0);
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
void Molecule::PerformOnePointRotToXPlus(Eigen::Vector3d point)
{
	double a = point(0), b = point(1), c = point(2);
	double alpha = atan2(c, b);
	double beta = atan2(sqrt(b*b + c*c), a);
	Eigen::Vector3d e_x(-1, 0, 0);
	Eigen::Vector3d e_z(0, 0, -1);
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
//angle in unit of radian
Eigen::Matrix3d Molecule::EulerRot(double a1_radian, double a2_radian, double a3_radian)
{
	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(a1_radian, Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(a2_radian, Eigen::Vector3d::UnitY())* Eigen::AngleAxisd(a3_radian, Eigen::Vector3d::UnitZ());
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
	//cout << x[0] << "\t" << x[1] << "\t" << x[2] << endl;
	Eigen::Matrix3d randomRotation;
	randomRotation = EulerRot(x[0]*PI/180, x[1]*PI/180, x[2]*PI/180);
	PerformTrans(-1 * Point);
	PerformRot(randomRotation);
	PerformTrans(Point);
}
//This could return the inertial angle value
void Molecule::PerformRandomRotEuler(Eigen::Vector3d Point, const double Precision_degree, double &ax_degree, double &ay_degree, double &az_degree)
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
	for (int i = 0; i != 3; i++)
	{
		x[i] = dis(generator)*Precision_degree;
	}
	//cout << x[0] << "\t" << x[1] << "\t" << x[2] << endl;
	Eigen::Matrix3d randomRotation;
	ax_degree = x[0]; ay_degree = x[1]; az_degree = x[2];
	randomRotation = EulerRot(x[0] * PI / 180, x[1] * PI / 180, x[2] * PI / 180);
	PerformTrans(-1 * Point);
	PerformRot(randomRotation);
	PerformTrans(Point);
}
Eigen::Matrix3d Molecule::EulerRotXY(double ax_radian, double ay_radian)
{
	Eigen::Matrix3d m;
	m = Eigen::AngleAxisd(ax_radian, Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(ay_radian, Eigen::Vector3d::UnitY());
	return m;
}
void Molecule::PerformRandomRotEulerXY(Eigen::Vector3d Point, const double &Precision_degree)
{
	double x[2];
	const int nums = (int)(360 / Precision_degree);

	clock_t now = clock();
	/*
	srand(now);
	for (int i = 0; i != 3; i++)
	x[i] = (rand() % nums + 1)*Precision_degree;
	*/

	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(0, nums);
	for (int i = 0; i != 2; i++)
	{
		x[i] = dis(generator)*Precision_degree;
	}
	//cout << x[0] << "\t" << x[1]  << endl;
	Eigen::Matrix3d randomRotation;
	randomRotation = EulerRotXY(x[0]*PI/180, x[1]*PI/180);
	PerformTrans(-1 * Point);
	PerformRot(randomRotation);
	PerformTrans(Point);
}
void Molecule::PerformRandomRotEulerXY(Eigen::Vector3d Point, const double &Precision_degree, double &ax_degree, double &ay_degree)
{
	double x[2];
	const int nums = (int)(360 / Precision_degree);

	clock_t now = clock();
	/*
	srand(now);
	for (int i = 0; i != 3; i++)
	x[i] = (rand() % nums + 1)*Precision_degree;
	*/

	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(0, nums);
	for (int i = 0; i != 2; i++)
	{
		x[i] = dis(generator)*Precision_degree;
	}
	//cout << x[0] << "\t" << x[1]  << endl;
	Eigen::Matrix3d randomRotation;
	ax_degree = x[0]; ay_degree = x[1];
	randomRotation = EulerRotXY(x[0] * PI / 180, x[1] * PI / 180);
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
double ClosestDistance(Molecule &ia, Molecule &ib)
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
	return distance;
}
double ClosestDistanceWithLabel(Molecule &ia, Molecule &ib, int &ia_label, int &ib_label)
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
      ia_label=label1;
      ib_label=label2;
	  return distance;
}
double AtomRadiusSum(Molecule &ia, Molecule &ib, int ia_label, int ib_label)
{
    int i,j;
    double radius_a,radius_b;
    for(i=0;i!=MaxElementOrder;i++)
    {
        if(ia.name[ia_label]==ELEMENT[i])
        break;
    }
    radius_a=ELEMENT_RADIUS[i];
     for(j=0;j!=MaxElementOrder;j++)
    {
        if(ib.name[ib_label]==ELEMENT[j])
        break;
    }
    radius_b=ELEMENT_RADIUS[j];
    return radius_a+radius_b;
}
void MakeAtomsSuitableDistanceMoveB(Molecule &ia, Molecule &ib, const double SmallestDistance)
{
	Eigen::Vector3d MC_A, MC_B,MC, MoveVector;
	MC_A = ia.MassCenter();
	while (ClosestDistance(ia, ib) < SmallestDistance)
	{
		//move B out of A
		MC_B = ib.MassCenter();
		MC = MC_B - MC_A;
		double length = sqrt(MC(0)*MC(0) + MC(1)*MC(1) + MC(2)*MC(2));
		MoveVector = 0.005 / length*MC;
		ib.PerformTrans(MoveVector);
	}
        while (ClosestDistance(ia,ib)>SmallestDistance)
        {
		//move B towards A
		MC_B=ib.MassCenter();
		MC=MC_A-MC_B;
 		double length = sqrt(MC(0)*MC(0) + MC(1)*MC(1) + MC(2)*MC(2));
                MoveVector = 0.005 / length*MC;
                ib.PerformTrans(MoveVector);
        }
}
void MakeAtomsUniformDistanceMoveB(Molecule &ia, Molecule &ib, const double times)
{
      Eigen::Vector3d MC_A, MC_B,MC, MoveVector;
	  MC_A = ia.MassCenter();
      int label1=0, label2=0;
      double closestdis=ClosestDistanceWithLabel(ia, ib,label1, label2);
	while (ClosestDistanceWithLabel(ia, ib,label1, label2) < times*AtomRadiusSum(ia,ib,label1,label2))
	{
		//move B out of A
		MC_B = ib.MassCenter();
		MC = MC_B - MC_A;
		double length = sqrt(MC(0)*MC(0) + MC(1)*MC(1) + MC(2)*MC(2));
		MoveVector = 0.005 / length*MC;
		ib.PerformTrans(MoveVector);
	}
        while (ClosestDistanceWithLabel(ia, ib,label1, label2)>times*AtomRadiusSum(ia,ib,label1,label2))
        {
		//move B towards A
		MC_B=ib.MassCenter();
		MC=MC_A-MC_B;
 		double length = sqrt(MC(0)*MC(0) + MC(1)*MC(1) + MC(2)*MC(2));
                MoveVector = 0.005 / length*MC;
                ib.PerformTrans(MoveVector);
        }
}
double RMSD( Molecule &ia,  Molecule &ib)
{
	if (ia.number != ib.number)
		return 100;
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
//This is to transform Molecule 2 coords to align well with Molecule 1, align first 3 atoms well
void SpaceTransform(Molecule ref, Molecule &change)
{
	;
}
void Molecule::AligenToStandardConfig()
{
	Eigen::Vector3d atom1, atom2, atom3;
	atom1 << corr[0][0],corr[0][1],corr[0][2];
	PerformTrans(-1 * atom1);//Make atom1 to origin
	atom2 << corr[1][0], corr[1][1], corr[1][2];
	PerformOnePointRotToXMinus(atom2); //make atom2 to x- axis
	atom3 << corr[2][0], corr[2][1], corr[2][2];
	double alpha = atan2(corr[2][2],corr[2][1]);
	Eigen::Vector3d axis;
	axis << -1, 0, 0;
	PerformAxisRot(axis, alpha);
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
bool operator>(const DoubleMolecule &ia, const DoubleMolecule &ib)
{
	if (ia.energy > ib.energy)
		return true;
	else
		return false;
}
bool operator<(const DoubleMolecule &ia, const DoubleMolecule &ib)
{
	if (ia.energy < ib.energy)
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
	cout << a<< b << endl;
}
double DoubleMolecule::Energy()
{
	return energy;
}
double RMSD(const DoubleMolecule ia, const DoubleMolecule ib)
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
void Fragments::PerformOnePointRotToXPlus(Eigen::Vector3d point)
{
	for (int i = 0; i != frag_number; i++)
		frags[i].PerformOnePointRotToXPlus(point);
}



//tinker .arc file analysis
void ReadFromWholeTinkerArc(const string arc_filename, const string save_filename, int solute_atoms, int each_solvent_atoms, int x_singal,int y_singal,int z_singal, double radius, double SoluteCenterToMarginLength, bool IfOutputStandardAtomName)
{
	//result save file name
	ofstream tofile(save_filename.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << save_filename << endl;
		exit(1);
	}
	//Save xyz info
	ofstream toxyz1((save_filename+"_near.arc").c_str(), ios::out);
	if (!toxyz1)
	{
		cerr << "Error to write " << save_filename << endl;
		exit(1);
	}
	ofstream toxyz2((save_filename+"_in.arc").c_str(), ios::out);
	if (!toxyz2)
	{
		cerr << "Error to write " << save_filename << endl;
		exit(1);
	}
	//read from Arc File
	ifstream infile(arc_filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << arc_filename << " to get the geometry info" << endl;
		exit(1);
	}
	tofile << "//This file is to read Tinker XYZ information from a huge .arc file and save some results" << endl;
	tofile << "//In each line, the first number represents the number of solvent molecules within the radius(" << radius << ") of solute mass center, the second number reperests the number of solvent molecules within the framework(" << SoluteCenterToMarginLength << ") of solute molecules" << endl;
	//Build a New class and read information from arc once and get many results
	SolventCube A;
	int x, y;
	while (!infile.eof())
	{
		A.ReadFromTinkerArcGetXYZonce(infile, solute_atoms, each_solvent_atoms);
		//Expand A to a larger one
		A.ExpandCubeTo8(x_singal,y_singal,z_singal);
		x = A.CountSolventNumberNearSolute(radius);
		y = A.CountSolventNumberInSolute(SoluteCenterToMarginLength);
		A.AppendToXYZSolventNearSolute(toxyz1, radius, IfOutputStandardAtomName);
		A.AppendToXYZSolventInSolute(toxyz2, SoluteCenterToMarginLength, IfOutputStandardAtomName);
		cout << x<< "\t" << y<< endl;
		tofile << x << "\t" <<y << endl;
	}
	infile.close();
	toxyz1.close();
	toxyz2.close();
	tofile.close();
}
void Do_ReadFromWholeTinkerArc_FromTxt()
{
	cout << "############################################" << endl;
	cout << "# Arc.Analysis program developed by Tao Li #" << endl;
	cout << "############################################" << endl << endl;;
	cout << "1. You should input a filename of task file, for example: /home/usr/task.txt" << endl;
	cout << "In the task file, you should write as below:" << endl << endl;
	cout << "//This is task file (this line is necessory in task file as a function of explanation)" << endl;
	cout << "water.arc  //arc filename where the program read from" << endl;
	cout << "result.txt //save filename, this program will output three files, *, *_near.arc, *_in.arc." << endl;
	cout << "108	3	//number of whole solute atoms, number of each solvent atoms(eg. H2O --> 3)" << endl;
	cout << "1	-1	-1	//this represents the direction you want to expand your box to 8 boxes.(eg. 1, -1, -1) If doesn't need expand, input 0 0 0" << endl;
	cout << "8.0		//this is the radius of the ball from the mass center of solute molecule" << endl;
	cout << "3.5		//this is the length of mass center to margin of the solute molecue." << endl;
	cout << "0		//0 means output atom name like CT HW, 1 means output standard atom name like C H" << endl << endl << endl;
	cout << "Please enter the task filename(eg. task.txt)" << endl;
	string task_filename;
	cin >> task_filename;
	ifstream infile(task_filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << task_filename << " to get the geometry info" << endl;
		exit(1);
	}
	
	string temp;
	string arc_filename, save_filename;
	int solute_atoms, each_solvent_atoms, x_singal, y_singal, z_singal;
	double radius, SoluteCenterToMarginLength;
	getline(infile, temp);
	infile >> arc_filename;
	getline(infile, temp);
	infile >> save_filename;
	getline(infile, temp);
	infile >> solute_atoms >> each_solvent_atoms;
	getline(infile, temp);
	infile >> x_singal >> y_singal >> z_singal;
	getline(infile, temp);
	infile >> radius;
	getline(infile, temp);
	infile >> SoluteCenterToMarginLength;
	getline(infile, temp);
	bool IfOutputStandardAtomName;
	infile >> IfOutputStandardAtomName;
	cout << endl << endl;
	cout << "This is your input file:"<< endl;
	cout << "arc_filename: " << arc_filename << endl;
	cout << "save_filename: " << save_filename << endl;
	cout << "solute_atoms: " << solute_atoms << ", each_solvent_atoms" << endl;
	cout << "x, y, z singal: " << x_singal << "\t" << y_singal << "\t" << z_singal << endl;
	cout << "radius: " << radius << endl;
	cout << "length from solute center to margin: " << SoluteCenterToMarginLength << endl;
	cout << "If output standard atom names: " << IfOutputStandardAtomName << endl<<endl;
	cout << "Program is busy now, please wait for good news" << endl;
	ReadFromWholeTinkerArc(arc_filename, save_filename, solute_atoms, each_solvent_atoms, x_singal, y_singal, z_singal, radius, SoluteCenterToMarginLength, IfOutputStandardAtomName);
	infile.close();
}
//This is to analysis the result(many *.xyz) of GenerateFunction(). Read from many *.xyz files, keep molecule A at one position, output different position of molecule B into one *.xyz file
void AligenMultiXYZfileToOneKeepMoleculeAStill(const string dir_name, int total_file_num, int molecule_A_atoms)
{
	Molecule ref,change;
	ref.ReadFromXYZfile((dir_name + "/0.xyz"));
	string savefile_name = dir_name + "/final.xyz";

	ofstream tofile(savefile_name.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << savefile_name << endl;
		exit(1);
	}
	tofile << molecule_A_atoms + total_file_num*(ref.Number() - molecule_A_atoms) << endl;
	tofile << "Analysis from " << total_file_num << " .xyz files" << endl;
	ref.ToXYZfileOnlyGeo(tofile, true);
	for (int i = 1; i != total_file_num; i++)
	{
		change.ReadFromXYZfile((dir_name + "/" + X_ToStr<int>(i) + ".xyz"));
		SpaceTransform(ref, change);
		Molecule temp;
		temp = AbstractSomeFormNewMolecule(change, molecule_A_atoms, ref.Number() - 1);
		temp.ToXYZfileOnlyGeo(tofile, true);
	}
	tofile.close();
}
void Do_AligenMultiXYZ_Program()
{
	cout << "####################################################################" << endl;
	cout << "# Result.analysis for GenerateFunction program developed by Tao Li #" << endl;
	cout << "####################################################################" << endl << endl;;
	cout << "Input a dir name where you save many .xyz file(0.xyz, 1.xyz, ...), eg: SaveConfigs" << endl;
	string dir_name;
	cin >> dir_name;
	cout << "Input total .xyz file numbers(from 0 to MaxValue):" << endl;
	int num;
	cin >> num;
	cout << "Input the number of atoms of first molecule: " << endl;
	int atoms;
	cin >> atoms;
	cout << "Program dealing..." << endl;
	AligenMultiXYZfileToOneKeepMoleculeAStill(dir_name, num, atoms);
	cout << "Done!" << endl;
}
void AlignEachXYZToStandardForm(const string dir_name, int total_file_num, int molecule_A_atoms)
{
	Molecule ref, change;
	ref.ReadFromXYZfile((dir_name + "/0.xyz"));
	string savefile_name = dir_name + "/final.xyz";

	ofstream tofile(savefile_name.c_str(), ios::out);
	if (!tofile)
	{
		cerr << "Error to write " << savefile_name << endl;
		exit(1);
	}
	tofile << molecule_A_atoms + total_file_num*(ref.Number() - molecule_A_atoms) << endl;
	tofile << "Analysis from " << total_file_num << " .xyz files" << endl;
	ref.AligenToStandardConfig();
	ref.ToXYZfileOnlyGeo(tofile, true);
	for (int i = 1; i != total_file_num; i++)
	{
		change.ReadFromXYZfile((dir_name + "/" + X_ToStr<int>(i) + ".xyz"));
		change.AligenToStandardConfig();
		Molecule temp;
		temp = AbstractSomeFormNewMolecule(change, molecule_A_atoms, ref.Number() - 1);
		temp.ToXYZfileOnlyGeo(tofile, true);
	}
	tofile.close();
}
void Do_AligenXYZStandard_Program()
{
	cout << "####################################################################" << endl;
	cout << "# Result.analysis for GenerateFunction program developed by Tao Li #" << endl;
	cout << "####################################################################" << endl << endl;;
	cout << "Input a dir name where you save many .xyz file(0.xyz, 1.xyz, ...), eg: SaveConfigs" << endl;
	string dir_name;
	cin >> dir_name;
	cout << "Input total .xyz file numbers(from 0 to MaxValue):" << endl;
	int num;
	cin >> num;
	cout << "Input the number of atoms of first molecule: " << endl;
	int atoms;
	cin >> atoms;
	cout << "Program dealing..." << endl;
	AlignEachXYZToStandardForm(dir_name, num, atoms);
	cout << "Done!" << endl;
}


void XYZToPDB_MoleculeAmBnType(const string xyz_filename, const string save_filename, int a_atoms, int a_num, int b_atoms, int b_num, int connect_m[][8], int connect_n[][8])
{
	ifstream infile(xyz_filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << xyz_filename << " to get the geometry info" << endl;
		exit(1);
	}
	int total_num;
	infile >> total_num;
	string temps;
	getline(infile, temps);
	getline(infile, temps);
	//judge if input ok
	if (total_num == a_atoms*a_num + b_atoms*b_num)
	{
		cout << "total num difference, xyz: " << total_num << ", you input: " << a_atoms*a_num + b_atoms*b_num;
	}
	vector<Molecule> As;
	vector<Molecule> Bs;
	Molecule temp;
	for (int i = 0; i != a_num; i++)
	{
		temp.ReadFromXYZOnlyGeo(infile, a_atoms);
		cout << temp;
		As.push_back(temp);
	}
	for (int i = 0; i != b_num; i++)
	{
		temp.ReadFromXYZOnlyGeo(infile, b_atoms);
		cout << temp;
		Bs.push_back(temp);
	}
	//write
	ToPDBfileAmBnType(save_filename, As, Bs, connect_m, connect_n);
	infile.close();

}
void ReadConnectionInfo(const string filename, int connect[][8],int atoms)
{
	for (int i = 0; i != MaxAtom; i++)
		for (int j = 0; j != 8; j++)
			connect[i][j] = 0;
	ifstream infile(filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << filename << " to get the geometry info" << endl;
		exit(1);
	}
	string temp;
	int temp_num;
	for (int i = 0; i != atoms; i++)
	{
		for (int j = 0; j != 8; j++)
		{
			infile >> temp_num;
			if (temp_num != 0)
				connect[i][j] = temp_num;
			else
			{
				getline(infile, temp);
				break;
			}
		}
	}
	infile.close();
	for (int i = 0; i != atoms; i++)
	{
		for (int j = 0; j != 8; j++)
		{
			if (connect[i][j] != 0)
				cout << connect[i][j] << "\t";
			else
			{
				cout << endl;
				break;
			}
		}
		
	}

}
void Do_XYZToPDB_MoleculeAmBnType()
{
	cout << "#################################" << endl;
	cout << "# xyz-->pdb developed by Tao Li #" << endl;
	cout << "#################################" << endl << endl;
	cout << "Please input xyz filename:" << endl;
	string xyz_filename,save_filename;
	cin >> xyz_filename;
	save_filename = xyz_filename + ".pdb";
	cout << "Please input atom numbers of molecule A and B (eg: 7 4):" << endl;
	int a_atoms, b_atoms, a_num, b_num;
	cin >> a_atoms >> b_atoms;
	cout << "Please input numbers of molecule A and B (eg: 1 36):" << endl;
	cin >> a_num >> b_num;
	cout << "Please input a filename to save connection info of A:" << endl;
	string filename1, filename2;
	cin >> filename1;
	int connect_m[MaxAtom][8], connect_n[MaxAtom][8];
	ReadConnectionInfo(filename1, connect_m, a_atoms);
	cout << "Please input a filename to save connection info of B:" << endl;
	cin >> filename2;
	
	ReadConnectionInfo(filename2, connect_n, b_atoms);
	cout << "Program dealing..." << endl;
	XYZToPDB_MoleculeAmBnType(xyz_filename, save_filename, a_atoms, a_num, b_atoms, b_num, connect_m, connect_n);
	cout << "Done!" << endl;
}

void XYZToMol2_MoleculeAmBnType(const string xyz_filename, const string save_filename, int a_atoms, int a_num, int b_atoms)
{
	ifstream infile(xyz_filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << xyz_filename << " to get the geometry info" << endl;
		exit(1);
	}
	int total_num;
	infile >> total_num;
	string temps;
	getline(infile, temps);
	getline(infile, temps);
	//judge if input ok
	int b_num = (total_num - a_atoms*a_num) / b_atoms;
	vector<Molecule> As;
	vector<Molecule> Bs;
	Molecule temp;
	for (int i = 0; i != a_num; i++)
	{
		temp.ReadFromXYZOnlyGeo(infile, a_atoms);
		As.push_back(temp);
	}
	for (int i = 0; i != b_num; i++)
	{
		temp.ReadFromXYZOnlyGeo(infile, b_atoms);
		Bs.push_back(temp);
	}
	//write
	ToMol2fileAmBnType(save_filename, As, Bs);
	infile.close();

}
void Do_XYZToMol2_MoleculeAmBnType()
{
	cout << "#################################" << endl;
	cout << "# xyz-->mol2 developed by Tao Li #" << endl;
	cout << "#################################" << endl << endl;
	cout << "Please input xyz filename:" << endl;
	string xyz_filename, save_filename;
	cin >> xyz_filename;
	save_filename = xyz_filename + ".mol2";
	cout << "Please input atom numbers of molecule A and B (eg: 7 4):" << endl;
	int a_atoms, b_atoms, a_num, b_num;
	cin >> a_atoms >> b_atoms;
	cout << "Please input numbers of molecule A(eg: 1):" << endl;
	cin >> a_num;
	XYZToMol2_MoleculeAmBnType(xyz_filename, save_filename, a_atoms, a_num, b_atoms);
	cout << "Done!" << endl;
}
