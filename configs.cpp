/*
 * =====================================================================================
 *
 *       Filename:  configs.cpp
 *
 *    Description:  Implementation of configs.h
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

#include "configs.h"
#include "molecule.h"
#include <time.h>
#include <queue>

//#define _NWCHEM_
#define _GAUSSIAN_
static double RandomNumber(double MaxValue)
{
	clock_t now = clock();
	srand(now);
	double x = (rand() % (int)((MaxValue*1000)) + 1);
	return abs(x)/1000;
}
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
	double x=(rand() % 1000 + 1)/50000;
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
#ifdef _GAUSSIAN_
	string filename("DATA/system");
	a.ToG09FileDFT(filename, basis, functional);
	//Use shell script to solve the scf energy
	double total_energy = -330;
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
	double x= (rand() % 1000 + 1)/50000;
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
	double total_energy = -660;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#endif
}

static void GenerateFunction(int matrix[][2], int index, int matrix2[][2], int index2, const int OutPutNumber, const string xyz_filename1, const string xyz_filename2)
{
	const double RotPrecision = 50;
	const double B1_default_value = 3.00;
	const int EachPairSaveNumber = 6;
	const int MaxRotTimes = EachPairSaveNumber * 10;
	const double RMSD_Precision = 0.35;
	Fragments FA, FB;
	FA.ReadFromXYZfile(xyz_filename1, index, matrix);
	FB.ReadFromXYZfile(xyz_filename2, index2, matrix2);
	cout << "Configuration of Molecule A:" << endl;
	cout << FA << endl;
	cout << "Configuration of Molecule B:" << endl;
	cout << FB << endl;
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
	for (int i = 0; i < total - 1; i++)
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

//Version 1.0, default 3.0 distance, potential partition function
static void GenerateFunction2(int matrix[][2],int index, int matrix2[][2],int index2, const int OutPutNumber,const string xyz_filename1,const string xyz_filename2)
{
	cout << "Enter Calculating..." << endl;
	const double RotPrecision = 20;
	const double B1_default_value = 3.00;
        const double Radius_Times=1.50;
	const double RMSD_Precision = 0.40;
	//for each pair config(ij[k]), rot * times
	const int EachSaveConfigRotTimes = 8;
	Fragments FA, FB;
	FA.ReadFromXYZfile(xyz_filename1, index, matrix);
	FB.ReadFromXYZfile(xyz_filename2, index2, matrix2);
	cout << "Configuration of Molecule A:" << endl;
	cout << FA << endl;
	cout << "Configuration of Molecule B:" << endl<<endl;
	cout << FB << endl;
	const double RestEnergies = G09energy(FA.TotalFragments()) + G09energy(FB.TotalFragments());
        cout<<"At rest, energy of 2 molecules is: "<<RestEnergies<<endl;

	//We need to find the sepcific EachPairSaveNumber(i,j) and MaxRotTimes(i,j) for each specific group combination according to partition function
	int MaxRotTimes[MAXFRAGMENT*2][MAXFRAGMENT];
	int EachPairSaveNumber[MAXFRAGMENT*2][MAXFRAGMENT];
	double PartitionFunction[MAXFRAGMENT*2][MAXFRAGMENT];
	//try and to find partition function
	double potential[MAXFRAGMENT*2][MAXFRAGMENT];//have no unit
	double total_partition = 0;
	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			Eigen::Vector3d MC_A1 = FA.ThisFragment().MassCenter();
			FA.PerformTrans(-1 * MC_A1);
			Eigen::Vector3d MC_A2 = FA.OtherFragments().MassCenter();
			FA.PerformOnePointRotToXMinus(MC_A2);//A this part at O, other part at x-
			Eigen::Vector3d MC_B1 = FB.ThisFragment().MassCenter();
			FB.PerformTrans(-1 * MC_B1);
			Eigen::Vector3d MC_B2 = FB.OtherFragments().MassCenter();
			FB.PerformOnePointRotToXPlus(MC_B2);//B this part at O, other part at x+
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;
			FB.PerformTrans(Default_B1);//B move at direction(1,0,0) to 3.00
			Molecule A = FA.TotalFragments();
			Molecule B = FB.TotalFragments();
			MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
			MC_A1 = A.MassCenter();
			MC_B1 = B.MassCenter();
			Eigen::Vector3d e_x;
			e_x << 1, 0, 0;
			double temp_potential[EachSaveConfigRotTimes];
			potential[i][j] = 0;
			for (int k = 0; k < EachSaveConfigRotTimes; k++)
			{
				Molecule tA = A;
				Molecule tB = B;
				tA.PerformAxisRot(e_x, RandomNumber(2 * PI));
				tB.PerformAxisRot(e_x, RandomNumber(2 * PI));
				temp_potential[k] = (G09energy(tA, tB) - RestEnergies)*HARTREE / K_B_BOLTZMAN / ROOM_TEMPERATURE;
				potential[i][j] += temp_potential[k];
			}
			potential[i][j] = potential[i][j] / EachSaveConfigRotTimes;
			cout << X_ToStr<int>(i) << "," << X_ToStr<int>(j) << " group-combination has potential: " << potential[i][j] <<"(unitless, V*Hartree/K_b/T)"<< endl;
			total_partition += exp(-1 * potential[i][j]);
		}
	}
    cout<<endl;
	//transfer potential matrix --> partition function matrix
    cout<<"Here is the partition function of each group-combination pair:"<<endl;
	for (int i = 0; i != FA.FragNumbers(); i++)
	{
        for (int j = 0; j != FB.FragNumbers(); j++)
	    {
		    PartitionFunction[i][j] = exp(-1 * potential[i][j]) / total_partition;
            cout<<PartitionFunction[i][j]<<"\t";
            if (PartitionFunction[i][j]>0.5)
                 EachPairSaveNumber[i][j] = 0.5*OutPutNumber;
		    else if (PartitionFunction[i][j]*OutPutNumber > 2)
			    EachPairSaveNumber[i][j] = PartitionFunction[i][j]*OutPutNumber;
		    else
			    EachPairSaveNumber[i][j] = 2;
		    MaxRotTimes[i][j] = EachPairSaveNumber[i][j] * EachSaveConfigRotTimes;
        }
        cout<<endl;
    }
	cout << endl << endl;

	//Begin operation 
	vector<DoubleMolecule> SaveSuitableCofigs;
	FA.IndexToZero();
	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			cout << "#Do calculate of " << FA.ThisFragment().MoleculeName() << " No." << i << ", and " << FB.ThisFragment().MoleculeName() << " No." << j << endl;
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
			for (int k = 0; k != MaxRotTimes[i][j]; k++)
			{
				//We should  rot A at the same time to make sure all suitable configurations happen!
                Molecule tA=A;
                Molecule tB=B;
				tA.PerformRandomRotEuler(MC_A1, RotPrecision*1.15);
				tB.PerformRandomRotEuler(MC_B1, RotPrecision);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				double  potential = G09energy(tA, tB) - RestEnergies;
				//cout << potential << "\t";
				DoubleMolecule temp_save;
				temp_save.Set(tA, tB, potential);
				cout << "Generate No."<<k<<" configuration with energy "<< temp_save.Energy() << endl;
				TempConfigs.push_back(temp_save);
			}
         
			//Sort configurations 
			for (int i1 = 0; i1 < MaxRotTimes[i][j]; i1++)
				for (int jj = i1+1; jj < MaxRotTimes[i][j]; jj++)
				{
					if (TempConfigs[jj] < TempConfigs[i1])
					{
                        			DoubleMolecule temp_config;
						temp_config = TempConfigs[jj];
						TempConfigs[jj] = TempConfigs[i1];
						TempConfigs[i1] = temp_config;
					}
				}

            cout<<"#After sorting,"<<endl;
            for(int ii = 0; ii < MaxRotTimes[i][j]; ii++)
            {
                cout<<"No."<<ii<<" config has energy "<<TempConfigs[ii].Energy()<<endl;
            }
			//Save Least Energy 5 configs and avoid rmsd similar one.
            cout << "Save " << EachPairSaveNumber[i][j] << " least energy configuration:" << endl;
			int output_count = 0;
			for (int ii = 0; output_count < EachPairSaveNumber[i][j]&&ii<MaxRotTimes[i][j]; ii++)
			{
					int total_size = SaveSuitableCofigs.size();
					int jj = 0;
					for (jj = 0; jj < total_size; jj++)
					{
						double temp_x = RMSD(SaveSuitableCofigs[jj], TempConfigs[ii]);
						if (abs(temp_x)< RMSD_Precision)
							break;
					}
					//if this config is different from other, saves
					if (jj == total_size)
					{
						SaveSuitableCofigs.push_back(TempConfigs[ii]);
						output_count += 1;
						//output this one to ./SaveConfigs/temp/ dir
						TempConfigs[ii].ToXYZ("SaveConfigs/temp/" + X_ToStr<int>(i) + "_" + X_ToStr<int>(j) + "_" + X_ToStr<int>(output_count) + ".xyz");
						TempConfigs[ii].output();
					}
			}
			cout <<"-------------------------------------------------" <<endl << endl;
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
	cout << "#Here output " << ((total > OutPutNumber) ? OutPutNumber : total) << " .xyz files to ./SaveConfigs/ as the final result" << endl;
	for (int i = 0; i < SaveSuitableCofigs.size() && i < OutPutNumber; i++)
		{
            SaveSuitableCofigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
            cout<<"No."<<i<<" least energy energy is: "<<SaveSuitableCofigs[i].Energy()<<endl;
        }
        
        //combine these .xyz files to one .xyz file 
        int total_file_num;
        if(total>OutPutNumber)
            total_file_num=OutPutNumber;
        else 
            total_file_num=total;
        AlignEachXYZToStandardForm("SaveConfigs", total_file_num, FA.TotalFragments().Number());
        cout<<"#Have generate an analysis .xyz file"<<endl;
		//translate .xyz --> .mol2
		XYZToMol2_MoleculeAmBnType("SaveConfigs/final.xyz","SaveConfigs/final.mol2", FA.TotalFragments().Number(),1, FB.TotalFragments().Number());
}



//Version 2.0, default 2.8 distance, has chance to change rotation of intramolecular structure. NBO charge exclusion
static void GenerateFunction3(int matrix[][2],int index, int matrix2[][2],int index2, const int OutPutNumber,const string xyz_filename1,const string xyz_filename2,bool Rotable1, bool Rotable2)
{
	cout << "Enter Calculating..." << endl;
	const double RotPrecision = 60;
	const double B1_default_value = 2.80;
        const double Radius_Times=1.50;
	const double RMSD_Precision = 0.40;
	//for each pair config(ij[k]), rot * times
	const int EachSaveConfigRotTimes = 10;
	Fragments FA, FB;
	FA.ReadFromXYZfile(xyz_filename1, index, matrix);
	FB.ReadFromXYZfile(xyz_filename2, index2, matrix2);
	cout << "Configuration of Molecule A:" << endl;
	cout << FA << endl;
	cout << "Configuration of Molecule B:" << endl<<endl;
	cout << FB << endl;
	const double RestEnergies = G09energy(FA.TotalFragments()) + G09energy(FB.TotalFragments());
        cout<<"At rest, energy of 2 molecules is: "<<RestEnergies<<endl;

	//We need to find the sepcific EachPairSaveNumber(i,j) and MaxRotTimes(i,j) for each specific group combination according to partition function
	int MaxRotTimes[MAXFRAGMENT*2][MAXFRAGMENT];
	int EachPairSaveNumber[MAXFRAGMENT*2][MAXFRAGMENT];
	double PartitionFunction[MAXFRAGMENT*2][MAXFRAGMENT];
	//try and to find partition function
	double potential[MAXFRAGMENT*2][MAXFRAGMENT];//have no unit
	double total_partition = 0;
	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			Eigen::Vector3d MC_A1 = FA.ThisFragment().MassCenter();
			FA.PerformTrans(-1 * MC_A1);
			Eigen::Vector3d MC_A2 = FA.OtherFragments().MassCenter();
			FA.PerformOnePointRotToXMinus(MC_A2);//A this part at O, other part at x-
			Eigen::Vector3d MC_B1 = FB.ThisFragment().MassCenter();
			FB.PerformTrans(-1 * MC_B1);
			Eigen::Vector3d MC_B2 = FB.OtherFragments().MassCenter();
			FB.PerformOnePointRotToXPlus(MC_B2);//B this part at O, other part at x+
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;
			FB.PerformTrans(Default_B1);//B move at direction(1,0,0) to 3.00
			Molecule A = FA.TotalFragments();
			Molecule B = FB.TotalFragments();
			MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
			MC_A1 = A.MassCenter();
			MC_B1 = B.MassCenter();
			Eigen::Vector3d e_x;
			e_x << 1, 0, 0;
			double temp_potential[5];
			potential[i][j] = 0;
			for (int k = 0; k < 5; k++)
			{
				Molecule tA = A;
				Molecule tB = B;
				tA.PerformAxisRot(e_x, RandomNumber(2 * PI));
				tB.PerformAxisRot(e_x, RandomNumber(2 * PI));
				temp_potential[k] = (G09energy(tA, tB) - RestEnergies)*HARTREE / K_B_BOLTZMAN / ROOM_TEMPERATURE;
				potential[i][j] += temp_potential[k];
			}
			potential[i][j] = potential[i][j] / 5;
			cout << X_ToStr<int>(i) << "," << X_ToStr<int>(j) << " group-combination has potential: " << potential[i][j] <<"(unitless, V*Hartree/K_b/T)"<< endl;
			total_partition += exp(-1 * potential[i][j]);
		}
	}
    cout<<endl;
	//transfer potential matrix --> partition function matrix
    cout<<"Here is the partition function of each group-combination pair:"<<endl;
	for (int i = 0; i != FA.FragNumbers(); i++)
	{
        for (int j = 0; j != FB.FragNumbers(); j++)
	    {
		    PartitionFunction[i][j] = exp(-1 * potential[i][j]) / total_partition;
            cout<<PartitionFunction[i][j]<<"\t";
            if (PartitionFunction[i][j]>0.5)
                 EachPairSaveNumber[i][j] = 0.5*OutPutNumber;
		    else if (PartitionFunction[i][j]*OutPutNumber > 2)
			    EachPairSaveNumber[i][j] = PartitionFunction[i][j]*OutPutNumber;
		    else
			    EachPairSaveNumber[i][j] = 2;
		    MaxRotTimes[i][j] = EachPairSaveNumber[i][j] * EachSaveConfigRotTimes;
        }
        cout<<endl;
    }
	cout << endl << endl;

	//Begin operation 
	vector<DoubleMolecule> SaveSuitableCofigs;
	FA.IndexToZero();
	for (int i = 0; i != FA.FragNumbers(); i++, ++FA)
	{
		FB.IndexToZero();
		for (int j = 0; j != FB.FragNumbers(); j++, ++FB)
		{
			cout << "#Do calculate of " << FA.ThisFragment().MoleculeName() << " No." << i << ", and " << FB.ThisFragment().MoleculeName() << " No." << j << endl;
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
			for (int k = 0; k != MaxRotTimes[i][j]; k++)
			{
				//We should  rot A at the same time to make sure all suitable configurations happen!
               			Molecule tA=A;
              			Molecule tB=B;
				if (Rotable1)
					RandomRotPossibleBond(tA, PI);
				if (Rotable2)
					RandomRotPossibleBond(tB, PI);
				tA.PerformRandomRotEuler(MC_A1, RotPrecision*1.15);
				tB.PerformRandomRotEuler(MC_B1, RotPrecision);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				double  potential = G09energy(tA, tB) - RestEnergies;
				//cout << potential << "\t";
				DoubleMolecule temp_save;
				temp_save.Set(tA, tB, potential);
				cout << "Generate No."<<k<<" configuration with energy "<< temp_save.Energy() << endl;
				TempConfigs.push_back(temp_save);
			}
         
			//Sort configurations 
			for (int i1 = 0; i1 < MaxRotTimes[i][j]; i1++)
				for (int jj = i1+1; jj < MaxRotTimes[i][j]; jj++)
				{
					if (TempConfigs[jj] < TempConfigs[i1])
					{
                        			DoubleMolecule temp_config;
						temp_config = TempConfigs[jj];
						TempConfigs[jj] = TempConfigs[i1];
						TempConfigs[i1] = temp_config;
					}
				}

            cout<<"#After sorting,"<<endl;
            for(int ii = 0; ii < MaxRotTimes[i][j]; ii++)
            {
                cout<<"No."<<ii<<" config has energy "<<TempConfigs[ii].Energy()<<endl;
            }
			//Save Least Energy 5 configs and avoid rmsd similar one.
            cout << "Save " << EachPairSaveNumber[i][j] << " least energy configuration:" << endl;
			int output_count = 0;
			for (int ii = 0; output_count < EachPairSaveNumber[i][j]&&ii<MaxRotTimes[i][j]; ii++)
			{
					int total_size = SaveSuitableCofigs.size();
					int jj = 0;
					for (jj = 0; jj < total_size; jj++)
					{
						double temp_x = RMSD(SaveSuitableCofigs[jj], TempConfigs[ii]);
						if (abs(temp_x)< RMSD_Precision)
							break;
					}
					//if this config is different from other, saves
					if (jj == total_size)
					{
						SaveSuitableCofigs.push_back(TempConfigs[ii]);
						output_count += 1;
						//output this one to ./SaveConfigs/temp/ dir
						TempConfigs[ii].ToXYZ("SaveConfigs/temp/" + X_ToStr<int>(i) + "_" + X_ToStr<int>(j) + "_" + X_ToStr<int>(output_count) + ".xyz");
						TempConfigs[ii].output();
					}
			}
			cout <<"-------------------------------------------------" <<endl << endl;
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
	cout << "#Here output " << ((total > OutPutNumber) ? OutPutNumber : total) << " .xyz files to ./SaveConfigs/ as the final result" << endl;
	for (int i = 0; i < SaveSuitableCofigs.size() && i < OutPutNumber; i++)
		{
            SaveSuitableCofigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
            cout<<"No."<<i<<" least energy energy is: "<<SaveSuitableCofigs[i].Energy()<<endl;
        }
        
        //combine these .xyz files to one .xyz file 
        int total_file_num;
        if(total>OutPutNumber)
            total_file_num=OutPutNumber;
        else 
            total_file_num=total;
        AlignEachXYZToStandardForm("SaveConfigs", total_file_num, FA.TotalFragments().Number());
        cout<<"#Have generate an analysis .xyz file"<<endl;
		//translate .xyz --> .mol2
		XYZToMol2_MoleculeAmBnType("SaveConfigs/final.xyz","SaveConfigs/final.mol2", FA.TotalFragments().Number(),1, FB.TotalFragments().Number());
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
	int con1[MAXFRAGMENT][2], con2[MAXFRAGMENT][2];
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
	getline(infile, temp);
	bool Rotable1,Rotable2;
	infile >> Rotable1>>Rotable2;
	cout << "Generate Max = " << OutputNum << " configurations to ./SaveConfigs" << endl;
	GenerateFunction3(con1, index1, con2, index2, OutputNum, xyzfile1, xyzfile2,Rotable1,Rotable2);
	infile.close();
}


void RandomGenerate(const int OutPutNumber, const string xyz_filename1, const string xyz_filename2)
{
    cout << "Enter Formation..." << endl;
	const double RotPrecision = 20;
	const double B1_default_value = 3.00;
	const double Radius_Times=1.50;
	const double RMSD_Precision = 0.35;
    Molecule A,B;
    A.ReadFromXYZfile(xyz_filename1);
    B.ReadFromXYZfile(xyz_filename2);
    cout<<A<<B<<endl;
    vector<DoubleMolecule> SaveConfigs;
    Eigen::Vector3d MC_A, MC_B;
    DoubleMolecule temp_config;
    for(int i=0;i<OutPutNumber;i++)
    {
        MC_A = A.MassCenter();
		MC_B = B.MassCenter();
        A.PerformRandomRotEuler(MC_A, RotPrecision);
        B.PerformRandomRotEuler(MC_B, RotPrecision);
        MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
        temp_config.Set(A,B,0.0);
        if(i==0)
            SaveConfigs.push_back(temp_config);
        else
            {
                int j=0;
                for(j=0;j<i;j++)
                 {
                    double temp_x = RMSD(temp_config, SaveConfigs[j]);
					if (abs(temp_x)< RMSD_Precision)
						break;
                 }
                 if(j==i)
                    SaveConfigs.push_back(temp_config);
                 else 
                    i--;
            }
        
    }
    for(int i=0;i<OutPutNumber;i++)
    {
        SaveConfigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
    } 
    //combine these .xyz files to one .xyz file 

      AlignEachXYZToStandardForm("SaveConfigs", OutPutNumber, A.Number());
        cout<<"#Have generate an analysis .xyz file"<<endl;
}
void Do_RandomGenerate_FromFile(string filename)
{
    cout<<"################################################"<<endl;
    cout<<"# You should prepare a file: random_task.txt   #"<<endl;
    cout<<"# ##This is the first line of random_task.txt  #"<<endl;
    cout<<"# InitiConfig/A.xyz                            #"<<endl;
    cout<<"# InitiConfig/B.xyz                            #"<<endl;
    cout<<"# 50   //This is output number in ./SaveConfigs#"<<endl;
    cout<<"################################################"<<endl;
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
	int OutputNum;
	infile >> OutputNum;
	cout << "Generate Max = " << OutputNum << " configurations to ./SaveConfigs" << endl;
	RandomGenerate(OutputNum, xyzfile1, xyzfile2);
	infile.close();
}
