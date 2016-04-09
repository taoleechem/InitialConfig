#include "molecule.h"
#include "configs.h"
#include "rmsd.h"

void change(Molecule reference, Molecule tip5p_model, Molecule &generate_one)
{
	generate_one = reference;
	Eigen::Matrix3d U;
	double ref[3][3], mov[3][3];
	int dim = 3;
	reference.Get2DArray(ref, dim);
	tip5p_model.Get2DArray(mov, dim);
	cout<<"Get Array from ref:"<<endl;
	for(int i=0; i<3; i++)
	{
		for(int j=0;j<3;j++)
			cout<<ref[i][j]<<" ";
		cout<<endl;
	}
	cout<<"Get Array from TIP5P:"<<endl;
	for(int i=0; i<3; i++)
	{
		for(int j=0;j<3;j++)						                        
			cout<<mov[i][j]<<" ";
		cout<<endl;
	}
	double  U2[3][3], mov_com[3], m2r[3], rmsd;
	calculate_rotation_rmsd(ref, mov, 3, mov_com, m2r, U2, rmsd);
	U << U2[0][0], U2[0][1], U2[0][2], U2[1][0], U2[1][1], U2[1][2], U2[2][0], U2[2][1], U2[2][2];
	Eigen::Vector3d MCA, MCB;
	MCA = reference.MassCenter();
	MCB = tip5p_model.MassCenter();
	cout<<"MC of tip5p is "<<MCB<<endl;
	cout<<"before move, tip5p is\n"<<tip5p_model<<endl;
	tip5p_model.PerformTrans(-1 * MCB);
	reference.PerformTrans(-1*MCA);
	tip5p_model.PerformRot(U);
	cout<<"After rot, the 2(r, tip5p) are:"<<endl;
	cout<<reference<<endl;
	cout<<tip5p_model<<endl;
	tip5p_model.PerformTrans(MCA);
	cout<<"Move Back, tip5p is:"<<endl;
	cout<<tip5p_model<<endl;
	int label[2] = { 4,5 };
	generate_one = generate_one + tip5p_model.NewMoleculeFromThis(2, label);
}

void WriteWaters(string filename, string tip5p_file, int filename_single_num)
{
	Molecule result;
	Molecule All;
	All.ReadFromTinkerXYZfile(filename);
	cout << All;
	Molecule tip5p;
	tip5p.ReadFromTinkerXYZfile(tip5p_file);
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
		Molecule tip5p1 = tip5p, generate1;
		change(x, tip5p1, generate1);
		if (result.Number() == 0)
			result = generate1;
		else
			result = result + generate1;
	}
	result.ToXYZfile("result.xyz");
	infile.close();
}


