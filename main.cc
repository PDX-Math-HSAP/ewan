
#include "CSRclass.h"
using namespace std;


int main()
{
  string filePath,destination,name,answer;
  bool dataInFile;
  bool useData = false;
  cout << "File must either consist of lines of the form 'i j', with i and j positive integers, or 'i j d', with data d a numeric type with 0,1,multiplication and addition. \n\r";
  cout << "Enter file path:\n\r";
  cin >> filePath;
  cout << "Does you file include data?\n\r";
  cin >> answer;
  if (answer=="yes")
  {
    dataInFile = true;
  }
  if (dataInFile)
  {
    cout << "Should I use this data?  If not, then each non zero entry will be assumed to be 1.\n\r";
    cin >> answer;
    if (answer == "yes")
    {
      useData = true;
    }
  }
  cout << "enter destination file path\n\r";
  cin >> destination;
  cout << "enter name of network\n\r"; 
  cin >> name; 
  cout << "importing\n\r";
  const CSRMat<double> bob(filePath,dataInFile, useData,"lap");
  size_t i=0;
  cout<< "initializing community list\n\r";
  CSRMat<double> proj = identity(bob,"right");
  bool keepGoing = true;
  while (keepGoing)
  {
	  cout<<"pass "<<i<<":\n\r";
    string a = to_string(i++);
    cout<<"->computing quotient laplacian of size "<< proj.cols<< " from laplacian of size "<< proj.rows <<"\n\r";
    CSRMat<double> quot = spMM(spMt(proj),spMM(bob,proj));
    cout<<"->removing isolated vertices\n\r";
    quot.removeZeros();
    cout<<"->building matchlist from "<< quot.rows <<" communities\n\r";
    CSRMat<size_t> vertEdge = generateRowByNnz(quot);
    CSRMat<size_t> edgeVert = generateNnzByRow(quot);
    Vec<size_t> match = parLuby(quot, edgeVert, vertEdge);
    cout<<"->generating projection into "<< quot.rows - match.size <<" communities\n\r";
    CSRMat<double> proj_temp = spMM(proj,contract(quot, edgeVert, vertEdge ,match));
    cout<<"->updating community list to "<< proj_temp.cols <<" communities\n\r";
    destroy(proj);
    proj = proj_temp;
    
    keepGoing = match.size>5;
    
    if ((i-1)%5==0 && keepGoing)
    {
      cout<<"->exporting\n\r";
      writeEdgeList(proj, destination + "/" + name + "com"+a+".txt","all", false, false);
    }
    cout<<"->checking loop condition  \n\r";
    
    if (!keepGoing)
    {
      writeEdgeList(proj, destination + "/" + name + "com"+a+".txt","all", false, false);
      writeEdgeList(quot, destination + "/" + name + "quotient"+a+".txt","all", false, false);
    }
    destroy(quot);
  }

  cout<<"Finished\n";

}