#include <iostream>
#include<vector>
#include<string.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include <fstream>
#include<cmath>
#include<sstream>
#include<iomanip>
using namespace std;
const int Max = 10000;
const double rfac=1.25;
struct coordi
{
    string name;
    double x;
    double y;
    double z;
};
int n;
typedef pair<int,double>tw;
string filename_luna;
coordi coo1[Max];
//vector<int> adj[Max];        //     1~N
vector<tw> QQ[Max];          //      1~N
double dist[Max][Max];
double rdist[Max][Max];
string elem[109]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K",
                  "Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
                  "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs",
                  "Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
                  "W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa",
                  "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt"};
                  //109 ge yuan su
double radius[109]={0.30, 1.16, 1.23, 0.89, 0.88,0.77, 0.70, 0.66, 0.58, 0.55,1.40, 1.36, 1.25, 1.17, 1.10,
                    1.11, 0.99, 1.58, 2.03, 1.74,1.44, 1.32, 1.20, 1.13, 1.17,1.16, 1.16, 1.15, 1.17, 1.25,
                    1.25, 1.22, 1.21, 1.17, 1.14,1.89, 2.25, 1.92, 1.62, 1.45,1.34, 1.29, 1.23, 1.24, 1.25,
                    1.28, 1.34, 1.41, 1.50, 1.40,1.41, 1.37, 1.33, 2.09, 2.35,1.98, 1.69, 1.65, 1.65, 1.64,
                    1.64, 1.66, 1.85, 1.61, 1.59,1.59, 1.58, 1.57, 1.56, 1.70,1.56, 1.44, 1.34, 1.30, 1.28,
                    1.26, 1.26, 1.29, 1.34, 1.44,1.55, 1.54, 1.52, 1.53, 1.52,1.53, 2.45, 2.02, 1.70, 1.63,
                    1.46, 1.40, 1.36, 1.25, 1.57,1.58, 1.54, 1.53, 1.84, 1.61,1.50, 1.49, 1.38, 1.36, 1.26,
                    1.20, 1.16, 1.14, 1.06 };
void remove_v(vector<tw>&v,tw val)
{
    vector<tw>::iterator ite;
    for(ite=v.begin();ite!=v.end();)
        {
        if((*ite).first==val.first)
            ite=v.erase(ite);
        else
            ++ite;
    }

}
void create_conn()
{
    double x,y,z;
    for(int i=1;i<n;i++)
     {
         for(int j=i+1;j<=n;j++)
         {
             x=coo1[i].x-coo1[j].x;
             y=coo1[i].y-coo1[j].y;
             z=coo1[i].z-coo1[j].z;
             dist[i][j]=sqrt(x*x+y*y+z*z);
             dist[j][i]=dist[i][j];
         }
     }
    string temp;
    int inuch[n+1];
    memset(inuch ,0,sizeof(inuch));
    for(int i=1;i<=n;i++)
    {
        temp=coo1[i].name;
        for(int j=0;j<109;j++)
        {
          if(temp==elem[j])
          {
             inuch[i]=j;
             break;
          }
        }
    }
    for(int i=1;i<n;i++)
    {
        int ii=inuch[i];
        for(int j=i+1;j<=n;j++)
        {
            int jj=inuch[j];
            rdist[i][j]=dist[i][j]/(radius[ii]+radius[jj]);
            rdist[j][i]=rdist[i][j];
        }
    }
 for(int i=1;i<n;i++)
    {
        for(int j=i+1;j<=n;j++)
        {
            if(rdist[i][j]>rfac)
                continue;
   //         adj[i].push_back(j);
   //         adj[j].push_back(i);
            QQ[i].push_back(make_pair(j,1.0));
            QQ[j].push_back(make_pair(i,1.0));
        }
    }
}
void read_xyz(string filename_xyz)
{
      char line[100]={0};
      ifstream fin(filename_xyz.c_str(),ios::in);
      fin>>n;
      cout<<n<<endl;
      fin.getline(line,sizeof(line));
      fin.getline(line,sizeof(line));
      string name_temp;
      double x ,y ,z;
      for (int i=1;i<=n;i++)
    {
        fin>>name_temp>>x>>y>>z;
        coo1[i].name=name_temp;
        coo1[i].x=x;
        coo1[i].y=y;
        coo1[i].z=z;
    }
}
void print_conn(int a,int na,int b,string xyzfilename, string savefilename)
{
    int ni;
    int nj;
    double aa=a;
    double bb=b;
    double T=a*na;
    ofstream fout(savefilename.c_str(),ios::out);
	fout << "#n B3LYP/6-31G(d) SP geom=connectivity" << endl << endl;
	fout << " Title" << endl<<endl;
	fout << "0 1" << endl;
	//
	char line[100] = { 0 };
	ifstream fin(xyzfilename.c_str(), ios::in);
	fin >> n;
	string temp;
	getline(fin, temp);
	getline(fin, temp);
	for (int i = 1; i <= n; i++)
	{
		getline(fin, temp);
		fout << temp<<endl;
	}
	fin.close();
	fout << endl;
    for(int i=1;i<=n;i++)
    {
         double ii=i;
         if(ii<T)
            {
                ni=na-floor((T-ii)/aa);
            }
            else
            {
                ni=na+ceil((ii-T)/bb);
            }
               fout<<i<<"  ";
        for(unsigned int j=0;j<QQ[i].size();j++)
        {
            double jj;
            jj=QQ[i][j].first;
             if(jj<=T)
            {
                nj=na-floor((T-jj)/aa);
            }
            else
            {
                nj=na+ceil((jj-T)/bb);
            }
            if(ni==nj)
            {
                  if(i<QQ[i][j].first)
                    {
                        fout<< setiosflags(ios::fixed)<<setprecision(1) <<QQ[i][j].first<<"  "<< setiosflags(ios::fixed)<<setprecision(1) <<QQ[i][j].second<<"  ";
                    }
            }
        }
        fout<<endl;
    }
	fout.close();

}
int main()
{
	cout << "Please enter xyz filename:" << endl;
	string xyzfilename;
	cin >> xyzfilename;
	cout << "Please enter atoms of molecule A" << endl;
	int a_atoms, a_num, b_atoms;
	cin >> a_atoms;
	cout << "Please enter atoms of molecule B" << endl;
	cin >> b_atoms;
	cout << "Please enter numbers of molecule A" << endl;
	cin >> a_num;

       read_xyz(xyzfilename);
	   string savefilename = xyzfilename + ".gjf";
       create_conn();
       print_conn( a_atoms, a_num, b_atoms, xyzfilename,savefilename);
	   return 0;
}
