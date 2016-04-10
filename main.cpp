/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Enter program
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
#include "configs.h"
#include "rmsd.h"

void FunctionAdvisor()
{
	cout << "This program provides some functions: " << endl;
	cout << "1. Generate many (50 as default) 2 molecules reaction complex .xyz files" << endl;
    cout << "2. Purely Random Generate many 2 molecules reaction complex .xyz files"<<endl;
	cout << "3. Analysis the result from 1, to get one .xyz file with one molecule A and multi-B, to find the space distribution of B near A" << endl;
	cout << "4. Do a tinker .arc file analysis" << endl;
	cout << "5. xyz --> mol2 with AmBn type" << endl;
	cout << "Please enter one numer to enter this sub program: " << endl;
	int num;
	cin >> num;
	if (num == 1)
		Do_GenerateFunction_Program_FromFile("task.txt");
    else if(num==2)
        Do_RandomGenerate_FromFile("random_task.txt");
	else if (num == 3)
		Do_AligenXYZStandard_Program();
	else if (num == 4)
		Do_ReadFromWholeTinkerArc_FromTxt();
	else if (num == 5)
		Do_XYZToMol2_MoleculeAmBnType();
}


int main(int argc, char* argv[])
{
		if(argc==1)
		FunctionAdvisor();
		else if (argc == 2)
		{
		Do_GenerateFunction_Program_FromFile(X_ToStr<char*>(argv[1]));
		}
	
	system("PAUSE");
    return 0;
}



