#include "configs.h"
#include "molecule.h"
#include <time.h>
//#define _NWCHEM_
#define _GAUSSIAN_
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
static double G09energy(Molecule a, string basis = "6-31g", string functional = "b3lyp")
{
#ifdef _WIN32
	clock_t now = clock();
	srand(now);
	double x=(rand() % 1000 + 1);
	return -1*abs(x);
#else
#ifdef _NWCHEM_
	string filename("DATA/NW.nw");
	a.ToNWchemFileDFT(filename,basis,functional);
	//Use shell script to solve the dft energy
	double total_energy = 0;
	system("./NW_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#ifdef _GAUSSIAN__
	string filename("DATA/system");
	a.ToG09FileDFT(filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#endif

}
static double G09energy(Molecule a, Molecule b, string basis = "6-31g", string functional = "b3lyp")
{
#ifdef _WIN32
	clock_t now = clock();
	srand(now);
	double x= (rand() % 1000 + 1);
	return -1 * abs(x);
#else
#ifdef _NWCHEM_
	string filename("DATA/NW.nw");
	ToNWchemFileDFT(a, b, filename, basis,functional);
	//Use shell script to solve the dft energy
	double total_energy = 0;
	system("./NW_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#ifdef _GAUSSIAN_
	string filename("DATA/system");
	ToG09FileDFT(a, b, filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#endif
}

static void GenerateFunction(int matrix[][2],int index, int matrix2[][2],int index2, const int OutPutNumber,const string xyz_filename1,const string xyz_filename2)
{
	const double RotPrecision = 50;
	const double B1_default_value = 3.00;
	const int EachPairSaveNumber = 6;
        const int MaxRotTimes=EachPairSaveNumber*10;
	const double RMSD_Precision = 0.35;
	Fragments FA, FB;
	FA.ReadFromXYZfile(xyz_filename1, index, matrix);
	FB.ReadFromXYZfile(xyz_filename2, index2, matrix2);
	cout << "Configuration of Molecule A:" << endl;
	cout << FA<<endl;
	cout << "Configuration of Molecule B:" << endl;
	cout << FB<<endl;
	const double RestEnergies = G09energy(FA.TotalFragments()) + G09energy(FB.TotalFragments());
	vector<DoubleMolecule> SaveSuitableCofigs;

	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			cout << "Do calculate of " << FA.ThisFragment().MoleculeName() << " No." << i << ", and " << FB.ThisFragment().MoleculeName() << " No." << j << endl;
			//Molecule A1 = FA.ThisFragment(); Molecule A2 = FA.OtherFragments(); Molecule B1 = FB.ThisFragment(); Molecule B2 = FB.OtherFragments();
			//translate A1-MC to origin
			Eigen::Vector3d MC_A1 = FA.ThisFragment().MassCenter();
			FA.PerformTrans(-1 * MC_A1);
			//keep A1-MC to origin, put A2-MC to x- axis
			Eigen::Vector3d MC_A2 = FA.OtherFragments().MassCenter();
			FA.PerformOnePointRotToXMinus(MC_A2);
			//translate B1-MC to x+ axis(default distance is 3.00), B2 at random position
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;//3.0, 0, 0
			Eigen::Vector3d MC_B1 = FB.ThisFragment().MassCenter();
			FB.PerformTrans(-1 * MC_B1 + Default_B1);
			//rot B at point B1_MC randomly
			Molecule A = FA.TotalFragments();
			Molecule B = FB.TotalFragments();
			MC_A1 = FA.ThisFragment().MassCenter();
			MC_B1 = FB.ThisFragment().MassCenter();
			vector<DoubleMolecule> TempConfigs;
			TempConfigs.clear();
			for (int k = 0; k != MaxRotTimes; k++)
			{
				//We should  rot A at the same time to make sure all suitable configurations happen!
				A.PerformRandomRotEuler(MC_A1, RotPrecision*1.3);
				B.PerformRandomRotEuler(MC_B1, RotPrecision);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
				double  potential = G09energy(A, B) - RestEnergies;
				//cout << potential << "\t";
				DoubleMolecule temp_save;
				temp_save.Set(A, B, potential);
				TempConfigs.push_back(temp_save);
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
			//Save Least Energy 5 configs and avoid rmsd similar one.
			int output_count = 0;
			for (int ii = 0; output_count < EachPairSaveNumber&&ii<MaxRotTimes; ii++)
			{
					int total_size = SaveSuitableCofigs.size();
					//initilize SaveSuitableConfigs
					if (total_size == 0)
					{
						SaveSuitableCofigs.push_back(TempConfigs[ii]);
						output_count += 1;
					}
					int jj = 0;
					for (jj = 0; jj < total_size; jj++)
					{
						double temp_x = RMSD(SaveSuitableCofigs[jj], TempConfigs[ii]);
						if (abs(temp_x)< RMSD_Precision)
							break;
					}
					if (jj == total_size)
					{
						SaveSuitableCofigs.push_back(TempConfigs[ii]);
						output_count += 1;
					}
			}
		}
	}

	//sort SaveSuitableConfigs
	int total = SaveSuitableCofigs.size();
	for (int i = 0; i < total-1; i++)
		for (int j = i + 1; j < total; j++)
		{
			if (SaveSuitableCofigs[i] > SaveSuitableCofigs[j])
			{
				DoubleMolecule temp;
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

static void GetGroupDevideInfoFromFile(const string filename,int con[][2], int label)
{
	ifstream infile(filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << filename << " to get the geometry info" << endl;
		exit(1);
	}
	string temp;
	for (int i = 0; i != label; i++)
	{
		infile >> con[i][0] >> con[i][1];
		getline(infile, temp);
	}
	infile.close();
}

void Do_GenerateFunction_Program_FromFile(string filename)
{
	
	ifstream infile(filename.c_str());
	if (!cout)
	{
		cerr << "Error to open " << filename << " to get the geometry info" << endl;
		exit(1);
	}
	string temp;
	getline(infile, temp);
	string xyzfile1, xyzfile2;
	infile >> xyzfile1;
	getline(infile, temp);
	infile >> xyzfile2;
	getline(infile, temp);
	cout << "Read molecule A from " << xyzfile1 << ", read molecule B from " << xyzfile2 << endl;
	int index1, index2;
	infile >> index1;
	getline(infile, temp);
	infile >> index2;
	getline(infile, temp);
	string groupdevide1, groupdevide2;
	infile >> groupdevide1;
	getline(infile, temp);
	infile >> groupdevide2;
	getline(infile, temp);
	int con1[MaxAtom][2], con2[MaxAtom][2];
	GetGroupDevideInfoFromFile(groupdevide1, con1, index1);
	GetGroupDevideInfoFromFile(groupdevide2, con2, index2);
	cout << "Molecule A Devide Group Method:" << endl;
	for (int i = 0; i != index1; i++)
	{
		cout << con1[i][0] << "\t" << con1[i][1] << endl;
	}
	cout << "Molecule B Devide Group Method:" << endl;
	for (int i = 0; i != index2; i++)
	{
		cout << con2[i][0] << "\t" << con2[i][1] << endl;
	}
	int OutputNum;
	infile >> OutputNum;
	cout << "Generate Max = " << OutputNum << " configurations to ./SaveConfigs" << endl;
	GenerateFunction(con1, index1, con2, index2, OutputNum, xyzfile1, xyzfile2);
	infile.close();
}


