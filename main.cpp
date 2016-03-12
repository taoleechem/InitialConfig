#include "molecule.h"
#include "configs.h"

void FunctionAdvisor()
{
	cout << "This program provides some functions: " << endl;
	cout << "1. Generate many (50 as default) 2 molecules reaction complex .xyz files" << endl;
	cout << "2. Analysis the result from 1, to get one .xyz file with one molecule A and multi-B, to find the space distribution of B near A" << endl;
	cout << "3. Do a tinker .arc file analysis" << endl;
	cout << "Please enter one numer to enter this sub program: " << endl;
	int num;
	cin >> num;
	if (num == 1)
		Do_GenerateFunction_Program_FromFile("task.txt");
	else if (num == 2)
		Do_AligenXYZStandard_Program();
	else if (num == 3)
		Do_ReadFromWholeTinkerArc_FromTxt();
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




