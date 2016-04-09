#include "molecule.h"
#include "configs.h"
#include "rmsd.h"

void change(Molecule reference, Molecule &tip5p_model)
{
	Eigen::Matrix3d U;
	double ref[3][3], mov[3][3];
	int dim = 3;
	reference.Get2DArray(ref, dim);
	tip5p_model.Get2DArray(mov, dim);
	double  U2[3][3], mov_com[3], m2r[3], rmsd;
	calculate_rotation_rmsd(ref, mov, 3, mov_com, m2r, U2, rmsd);
	U << U2[0][0], U2[0][1], U2[0][2], U2[1][0], U2[1][1], U2[1][2], U2[2][0], U2[2][1], U2[2][2];
	Eigen::Vector3d MCA, MCB;
	MCB = tip5p_model.MassCenter();
	tip5p_model.PerformTrans(-1 * MCB);
	tip5p_model.PerformRot(U);
	tip5p_model.PerformTrans(MCB);
	cout << tip5p_model << endl;
}

void WriteWaters(string filename, string tip5p_file, int filename_single_num)
{
	Molecule result;
	Molecule All;
	All.ReadFromTinkerXYZfile(filename);
	cout << All;
	Molecule tip5p;
	tip5p.ReadFromXYZfile(tip5p_file);
	cout << tip5p;
	int water_num = All.Number() / filename_single_num;
	cout << "There are " << water_num << " H2O in " << filename << endl;
	ifstream infile(filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << filename << " to get the geometry info" << endl;
		exit(1);
	}
	string temp;
	getline(infile, temp);
	Molecule x;
	for (int i = 0; i < water_num; i++)
	{
		cout << "Treating No." << i + 1 << " H2O" << endl;
		x.ReadFromTinkerXYZGeoPart(infile, filename_single_num);
		Molecule tip5p1 = tip5p;
		change(x, tip5p1);
		if (result.Number() == 0)
			result = tip5p1;
		else
			result = result + tip5p1;
	}
	result.ToXYZfile("result.xyz");
	infile.close();
}


