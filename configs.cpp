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

#define _NWCHEM_
//#define _GAUSSIAN_

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
static double G09energy(Molecule a, string basis = "6-31g", string functional = "b3lyp", string othercommand = "")
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
	a.ToG09FileDFT(filename, basis, functional, othercommand);
	//Use shell script to solve the scf energy
	double total_energy = -330;
	system("./g09_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	return total_energy;
#endif
#endif

}
static double G09energy(Molecule a, Molecule b, string basis = "6-31g", string functional = "b3lyp", string othercommand = "")
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
	ToG09FileDFT(a, b, filename, basis, functional, othercommand);
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
class EulerAngle5
{
public:
	int ax, ay, bx, by, bz, modei, modej;
	EulerAngle5(int i1, int i2, int i3, int i4, int i5,int m1, int m2) :ax(i1), ay(i2), bx(i3), by(i4), bz(i5),modei(m1),modej(m2) {}
	friend ostream& operator<<(ostream &os, EulerAngle5 &A)
	{
		os <<A.modei<<"\t"<<A.modej<<"\t"<< A.ax << "\t" << A.ay << "\t" << A.bx << "\t" << A.by << "\t" << A.bz;
		return os;
	}

};
//Maps is used to save each calculated configuration
class Maps
{
private:
	vector<EulerAngle5> EulerAngle;
	vector<double> energy;
public:
	Maps() {}
	Maps(Maps &iA)
	{
		EulerAngle = iA.EulerAngle;
		energy = iA.energy;
	}
	void AddPoint(int modei, int modej, int ax, int ay, int bx, int by, int bz, double iE)
	{
		EulerAngle5 temp(ax, ay, bx, by, bz, modei, modej);
		EulerAngle.push_back(temp);
		energy.push_back(iE);
	}
	friend ostream& operator<<(ostream &os, Maps &A)
	{
		vector<EulerAngle5>::iterator iter1;
		vector<double>::iterator iter2;
		for (iter1 = A.EulerAngle.begin(), iter2 = A.energy.begin(); iter1 != A.EulerAngle.end() && iter2 != A.energy.end(); iter1++, iter2++)
		{
			os << *iter1 << "\t" << *iter2 << endl;
		}
		return os;
	}
	double Energy(int label)
	{
		if (label >= 0 && label < energy.size())
			return energy[label];
		else
			return 10242048;
	}
	int WhichLabel(int modei,int modej, int ax, int ay, int bx, int by, int bz)
	{
		vector<EulerAngle5>::iterator iter1;
		int i = 0;
		for (i = 0, iter1 = EulerAngle.begin(); iter1 != EulerAngle.end(); iter1++, i++)
		{
			if ((*iter1).modej == modej&& (*iter1).modei == modei&&abs((*iter1).ax - ax)<=1 && abs((*iter1).ay - ay)<=1 && abs((*iter1).bx - bx)<=1 && abs((*iter1).by - by)<=1 && abs((*iter1).bz - bz)<=1)
				return i;
		}
		return -1;
	}
};
void GenerateFunction3(int matrix[][2],int index, int matrix2[][2],int index2, const int OutPutNumber,const string xyz_filename1,const string xyz_filename2,int Temperature, bool Rotable1, bool Rotable2, double B1_default_value)
{
	cout << "Temperature is " << Temperature << " K" << endl;
	cout << "Distance between closest atoms of 2 molecules is "<< B1_default_value<<endl;
	cout << "Enter Calculating..." << endl;
	Maps SaveCalculations;//save each calculation value
	if (Temperature == 0)
		Temperature = ROOM_TEMPERATURE;
	const double RotPrecision = 20;
	//const double B1_default_value = 2.80;
	if(B1_default_value ==0)
		B1_default_value = 2.80;
 	   const double Radius_Times=1.50;
	const double RMSD_Precision = 0.5;
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
	//calculate average potential between each combination mode
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
			double temp_potential[6];
			potential[i][j] = 0;
			for (int k = 0; k < 6; k++)
			{
				Molecule tA = A;
				Molecule tB = B;
				tB.PerformAxisRot(e_x, k*PI /3);
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				temp_potential[k] = (G09energy(tA, tB) - RestEnergies)*HARTREE / K_B_BOLTZMAN / Temperature;
				SaveCalculations.AddPoint(i, j, 0, 0, k * 60, 0, 0, temp_potential[k]);//Save each calculation value
				potential[i][j] += temp_potential[k];
			}
			potential[i][j] = potential[i][j] / 6;
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
        }
        cout<<endl;
    }
	cout << "Here is the partition number N_{ij} of each group-combination pair:" << endl;
	for (int i = 0; i != FA.FragNumbers(); i++)
	{
		for (int j = 0; j != FB.FragNumbers(); j++)
		{
			if (PartitionFunction[i][j]>0.5)
				EachPairSaveNumber[i][j] = 0.5*OutPutNumber;
			else if (PartitionFunction[i][j] * OutPutNumber > 2)
				EachPairSaveNumber[i][j] = PartitionFunction[i][j] * OutPutNumber;
			else
				EachPairSaveNumber[i][j] = 2;
			cout << EachPairSaveNumber[i][j] << "\t";
			MaxRotTimes[i][j] = EachPairSaveNumber[i][j] * EachSaveConfigRotTimes;
		}
		cout << endl;
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
			//translate B1-MC to x+ axis(default distance is 3.00), B2 at x- axis
			Eigen::Vector3d MC_B1 = FB.ThisFragment().MassCenter();
			FB.PerformTrans(-1 * MC_B1);
			Eigen::Vector3d MC_B2 = FB.OtherFragments().MassCenter();
			FB.PerformOnePointRotToXPlus(MC_B2);//B this part at O, other part at x+
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;
			FB.PerformTrans(Default_B1);//B move at direction(1,0,0) to 3.00
			//rot B at point B1_MC randomly
			Molecule A = FA.TotalFragments();
			Molecule B = FB.TotalFragments();
			MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
			MC_A1 = FA.ThisFragment().MassCenter();
			MC_B1 = FB.ThisFragment().MassCenter();
			vector<DoubleMolecule> TempConfigs;
			TempConfigs.clear();
			for (int k = 0; k != MaxRotTimes[i][j]; k++)
			{
				//We should  rot A at the same time to make sure all suitable configurations happen!
               	Molecule tA=A;
              	Molecule tB=B;
				double ax, ay, bx, by, bz;
				tA.PerformRandomRotEulerXY(MC_A1, RotPrecision,ax,ay);
				tB.PerformRandomRotEuler(MC_B1, RotPrecision,bx,by,bz);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				int search_label = SaveCalculations.WhichLabel(i, j, ax, ay, bx, by, bz);
				double  potential;
				if (search_label == -1)
				{
					potential = G09energy(tA, tB) - RestEnergies;
					SaveCalculations.AddPoint(i, j, ax, ay, bx, by, bz, potential);
				}
				else
				{
					potential = SaveCalculations.Energy(search_label);
					cout<<"Find existing value"<<endl;
					k--;
					continue;
				}
					
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
	//and do bond rotation analysis at same time
	if (Rotable1 ==true || Rotable2 ==true )
	{
		cout << "Do bond rotation to top " << OutPutNumber << " configurations" << endl;
	}
	cout << "#Here output " << ((total > OutPutNumber) ? OutPutNumber : total) << " .xyz files to ./SaveConfigs/ as the final result" << endl;
	for (int i = 0; i < SaveSuitableCofigs.size() && i < OutPutNumber; i++)
		{
			if (Rotable1 == true || Rotable2 == true)
			{
				Molecule t1, t2;
				double ie;
				SaveSuitableCofigs[i].GetInfo(t1, t2, ie);
				//rot analysis
				if (Rotable1 == true)
					RandomRotPossibleBond(t1, PI);
				if (Rotable2 == true)
					RandomRotPossibleBond(t2, PI);
				MakeAtomsSuitableDistanceMoveB(t1, t2, B1_default_value);
				double tE = G09energy(t1, t2) - RestEnergies;
				//random number
				clock_t now = clock();
				/*
				srand(now);
				for (int i = 0; i != 3; i++)
				x[i] = (rand() % nums + 1)*Precision_degree;
				*/
				std::default_random_engine generator(now);
				std::uniform_real_distribution<double> dis(0, 1);
				double randomD = dis(generator);
				double probability=exp((ie - tE)*HARTREE / K_B_BOLTZMAN / 3000);
				cout << "\t after rot, potential is: " << tE<<" while before potential is "<<ie<<endl;
				cout<<" \t probability from exp(E) is "<<probability <<" random number is "<<randomD<< endl;
				if (randomD< probability)
				{
					SaveSuitableCofigs[i].Set(t1, t2, tE);
					cout << "\t Sucessfully rot bond to No." << i << " configuration" << endl;
				}
			}
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
        cout<<"#Have generate an analysis .xyz and .mol2 file"<<endl;
		//translate .xyz --> .mol2
		XYZToMol2_MoleculeAmBnType("SaveConfigs/final.xyz","SaveConfigs/final.mol2", FA.TotalFragments().Number(),1, FB.TotalFragments().Number());
		cout << "At last, show saved calculations:" << endl;
		cout << SaveCalculations << endl;
		cout << "Done" << endl;
}

//usr could change basis/functional and other command
void GenerateFunction4(int matrix[][2], int index, int matrix2[][2], int index2, const int OutPutNumber, const string xyz_filename1, const string xyz_filename2, int Temperature, bool Rotable1, bool Rotable2, double B1_default_value, string functional = "b3lyp", string basis = "6-31G", string othercommand="")
{
	cout << "Temperature is " << Temperature << " K" << endl;
	cout << "Distance between closest atoms of 2 molecules is " << B1_default_value << endl;
	cout << "Enter Calculating..." << endl;
	Maps SaveCalculations;//save each calculation value
	if (Temperature == 0)
		Temperature = ROOM_TEMPERATURE;
	const double RotPrecision = 20;
	//const double B1_default_value = 2.80;
	if (B1_default_value == 0)
		B1_default_value = 2.80;
	const double Radius_Times = 1.50;
	const double RMSD_Precision = 0.5;
	//for each pair config(ij[k]), rot * times
	const int EachSaveConfigRotTimes = 8;
	Fragments FA, FB;
	FA.ReadFromXYZfile(xyz_filename1, index, matrix);
	FB.ReadFromXYZfile(xyz_filename2, index2, matrix2);
	cout << "Configuration of Molecule A:" << endl;
	cout << FA << endl;
	cout << "Configuration of Molecule B:" << endl << endl;
	cout << FB << endl;
	const double RestEnergies = G09energy(FA.TotalFragments(), basis, functional, othercommand) + G09energy(FB.TotalFragments(),basis, functional, othercommand);
	cout << "At rest, energy of 2 molecules is: " << RestEnergies << endl;

	//We need to find the sepcific EachPairSaveNumber(i,j) and MaxRotTimes(i,j) for each specific group combination according to partition function
	int MaxRotTimes[MAXFRAGMENT * 2][MAXFRAGMENT];
	int EachPairSaveNumber[MAXFRAGMENT * 2][MAXFRAGMENT];
	double PartitionFunction[MAXFRAGMENT * 2][MAXFRAGMENT];
	//try and to find partition function
	double potential[MAXFRAGMENT * 2][MAXFRAGMENT];//have no unit
	double total_partition = 0;
	//calculate average potential between each combination mode
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
			double temp_potential[6];
			potential[i][j] = 0;
			for (int k = 0; k < 6; k++)
			{
				Molecule tA = A;
				Molecule tB = B;
				tB.PerformAxisRot(e_x, k*PI / 3);
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				temp_potential[k] = (G09energy(tA, tB, basis, functional, othercommand) - RestEnergies)*HARTREE / K_B_BOLTZMAN / Temperature;
				SaveCalculations.AddPoint(i, j, 0, 0, k * 60, 0, 0, temp_potential[k]);//Save each calculation value
				potential[i][j] += temp_potential[k];
			}
			potential[i][j] = potential[i][j] / 6;
			cout << X_ToStr<int>(i) << "," << X_ToStr<int>(j) << " group-combination has potential: " << potential[i][j] << "(unitless, V*Hartree/K_b/T)" << endl;
			total_partition += exp(-1 * potential[i][j]);
		}
	}
	cout << endl;
	//transfer potential matrix --> partition function matrix
	cout << "Here is the partition function of each group-combination pair:" << endl;
	for (int i = 0; i != FA.FragNumbers(); i++)
	{
		for (int j = 0; j != FB.FragNumbers(); j++)
		{
			PartitionFunction[i][j] = exp(-1 * potential[i][j]) / total_partition;
			cout << PartitionFunction[i][j] << "\t";
		}
		cout << endl;
	}
	cout << "Here is the partition number N_{ij} of each group-combination pair:" << endl;
	for (int i = 0; i != FA.FragNumbers(); i++)
	{
		for (int j = 0; j != FB.FragNumbers(); j++)
		{
			if (PartitionFunction[i][j]>0.5)
				EachPairSaveNumber[i][j] = 0.5*OutPutNumber;
			else if (PartitionFunction[i][j] * OutPutNumber > 2)
				EachPairSaveNumber[i][j] = PartitionFunction[i][j] * OutPutNumber;
			else
				EachPairSaveNumber[i][j] = 2;
			cout << EachPairSaveNumber[i][j] << "\t";
			MaxRotTimes[i][j] = EachPairSaveNumber[i][j] * EachSaveConfigRotTimes;
		}
		cout << endl;
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
			//translate B1-MC to x+ axis(default distance is 3.00), B2 at x- axis
			Eigen::Vector3d MC_B1 = FB.ThisFragment().MassCenter();
			FB.PerformTrans(-1 * MC_B1);
			Eigen::Vector3d MC_B2 = FB.OtherFragments().MassCenter();
			FB.PerformOnePointRotToXPlus(MC_B2);//B this part at O, other part at x+
			Eigen::Vector3d Default_B1;
			Default_B1 << B1_default_value, 0, 0;
			FB.PerformTrans(Default_B1);//B move at direction(1,0,0) to 3.00
										//rot B at point B1_MC randomly
			Molecule A = FA.TotalFragments();
			Molecule B = FB.TotalFragments();
			MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
			MC_A1 = FA.ThisFragment().MassCenter();
			MC_B1 = FB.ThisFragment().MassCenter();
			vector<DoubleMolecule> TempConfigs;
			TempConfigs.clear();
			for (int k = 0; k != MaxRotTimes[i][j]; k++)
			{
				//We should  rot A at the same time to make sure all suitable configurations happen!
				Molecule tA = A;
				Molecule tB = B;
				double ax, ay, bx, by, bz;
				tA.PerformRandomRotEulerXY(MC_A1, RotPrecision, ax, ay);
				tB.PerformRandomRotEuler(MC_B1, RotPrecision, bx, by, bz);
				//Here need to adjust B to a suitable position that the closest distance between atoms of A and B is 3.0
				MakeAtomsSuitableDistanceMoveB(tA, tB, B1_default_value);
				int search_label = SaveCalculations.WhichLabel(i, j, ax, ay, bx, by, bz);
				double  potential;
				if (search_label == -1)
				{
					potential = G09energy(tA, tB, basis, functional, othercommand) - RestEnergies;
					SaveCalculations.AddPoint(i, j, ax, ay, bx, by, bz, potential);
				}
				else
				{
					potential = SaveCalculations.Energy(search_label);
					cout << "Find existing value" << endl;
					k--;
					continue;
				}

				//cout << potential << "\t";
				DoubleMolecule temp_save;
				temp_save.Set(tA, tB, potential);
				cout << "Generate No." << k << " configuration with energy " << temp_save.Energy() << endl;
				TempConfigs.push_back(temp_save);
			}

			//Sort configurations 
			for (int i1 = 0; i1 < MaxRotTimes[i][j]; i1++)
				for (int jj = i1 + 1; jj < MaxRotTimes[i][j]; jj++)
				{
					if (TempConfigs[jj] < TempConfigs[i1])
					{
						DoubleMolecule temp_config;
						temp_config = TempConfigs[jj];
						TempConfigs[jj] = TempConfigs[i1];
						TempConfigs[i1] = temp_config;
					}
				}

			//Save Least Energy 5 configs and avoid rmsd similar one.
			cout << "Save " << EachPairSaveNumber[i][j] << " least energy configuration:" << endl;
			int output_count = 0;
			for (int ii = 0; output_count < EachPairSaveNumber[i][j] && ii<MaxRotTimes[i][j]; ii++)
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
			cout << "-------------------------------------------------" << endl << endl;
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
	//and do bond rotation analysis at same time
	if (Rotable1 == true || Rotable2 == true)
	{
		cout << "Do bond rotation to top " << OutPutNumber << " configurations" << endl;
	}
	cout << "#Here output " << ((total > OutPutNumber) ? OutPutNumber : total) << " .xyz files to ./SaveConfigs/ as the final result" << endl;
	for (int i = 0; i < SaveSuitableCofigs.size() && i < OutPutNumber; i++)
	{
		if (Rotable1 == true || Rotable2 == true)
		{
			Molecule t1, t2;
			double ie;
			SaveSuitableCofigs[i].GetInfo(t1, t2, ie);
			//rot analysis
			if (Rotable1 == true)
				RandomRotPossibleBond(t1, PI);
			if (Rotable2 == true)
				RandomRotPossibleBond(t2, PI);
			MakeAtomsSuitableDistanceMoveB(t1, t2, B1_default_value);
			double tE = G09energy(t1, t2) - RestEnergies;
			//random number
			clock_t now = clock();
			/*
			srand(now);
			for (int i = 0; i != 3; i++)
			x[i] = (rand() % nums + 1)*Precision_degree;
			*/
			std::default_random_engine generator(now);
			std::uniform_real_distribution<double> dis(0, 1);
			double randomD = dis(generator);
			double probability = exp((ie - tE)*HARTREE / K_B_BOLTZMAN / 3000);
			cout << "\t after rot, potential is: " << tE << " while before potential is " << ie << endl;
			cout << " \t probability from exp(E) is " << probability << " random number is " << randomD << endl;
			if (randomD< probability)
			{
				SaveSuitableCofigs[i].Set(t1, t2, tE);
				cout << "\t Sucessfully rot bond to No." << i << " configuration" << endl;
			}
		}
		SaveSuitableCofigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
		cout << "No." << i << " least energy energy is: " << SaveSuitableCofigs[i].Energy() << endl;
	}


	//combine these .xyz files to one .xyz file 
	int total_file_num;
	if (total>OutPutNumber)
		total_file_num = OutPutNumber;
	else
		total_file_num = total;
	AlignEachXYZToStandardForm("SaveConfigs", total_file_num, FA.TotalFragments().Number());
	cout << "#Have generate an analysis .xyz and .mol2 file" << endl;
	//translate .xyz --> .mol2
	XYZToMol2_MoleculeAmBnType("SaveConfigs/final.xyz", "SaveConfigs/final.mol2", FA.TotalFragments().Number(), 1, FB.TotalFragments().Number());
	cout << "At last, show saved calculations:" << endl;
	cout << SaveCalculations << endl;
	cout << "Done" << endl;
}

void GetGroupDevideInfoFromFile(const string filename,int con[][2], int label)
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
	int temperature;
	infile >> temperature;
	getline(infile, temp);
	double B1_default_value;
	infile >> B1_default_value;
	getline(infile, temp);
	bool Rotable1,Rotable2;
	infile >> Rotable1>>Rotable2;
	string basis, functional, othercommand;
	infile >> functional;
	if(functional.size()<2)
	{
		functional = "b3lyp";
		basis = "6-31g";
		othercommand = " ";
	}
	else
	{

		infile >> basis;
		getline(infile, othercommand);
	}
	cout << "Generate Max = " << OutputNum << " configurations to ./SaveConfigs" << endl;
	GenerateFunction4(con1, index1, con2, index2, OutputNum, xyzfile1, xyzfile2,temperature, Rotable1,Rotable2, B1_default_value,functional, basis, othercommand);
	infile.close();
}


void RandomGenerate(const int OutPutNumber, const string xyz_filename1, const string xyz_filename2)
{
    cout << "Enter Formation..." << endl;
	const double RotPrecision = 5;
	const double B1_default_value = 2.80;
	const double Radius_Times=1.50;
	const double RMSD_Precision = 1.00;
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
        A.PerformRandomRotEulerXY(MC_A, RotPrecision);
        B.PerformRandomRotEuler(MC_B, RotPrecision);
        MakeAtomsSuitableDistanceMoveB(A, B, B1_default_value);
        temp_config.Set(A,B,0.0);
		SaveConfigs.push_back(temp_config);
    }
    for(int i=0;i<OutPutNumber;i++)
    {
        SaveConfigs[i].ToXYZ("SaveConfigs/" + X_ToStr<int>(i) + ".xyz");
    } 
    //combine these .xyz files to one .xyz file 

      AlignEachXYZToStandardForm("SaveConfigs", OutPutNumber, A.Number());
	  XYZToMol2_MoleculeAmBnType("SaveConfigs/final.xyz", "SaveConfigs/final.mol2", A.Number(), 1, B.Number());
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
