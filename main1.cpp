#include"class.h"
#define NUMBER 100000000
#define COUNTER 5000

void test1()
{
vector<CFabrics*> fabric;
	vector<CComplexVector*> w;
	fabric.push_back(new CFabrics1);
	fabric.push_back(new CFabrics2);
	ifstream in("in.txt"); 
	string str;
	while(getline(in,str)){
		int Type; 
		string outfile;
		vector<ComplexNumber> vectmp;
		stringstream ss(str);
		ss >> Type >> outfile;
		for(ComplexNumber tmp;ss >> tmp;){
			vectmp.push_back(tmp);
		}
		CComplexVector *t=fabric[Type-1]->Create();
		*t=vectmp;
		t->input(outfile);
		w.push_back(t);
	}
	for(size_t i=0;i<w.size();i++){
		delete w[i];
	}
	for(size_t i=0;(i<fabric.size());i++)  {delete fabric[i];}
}

void test2()
{
try {
   vector<CFabrics*> fabric;
   vector<CComplexVector*> v;
   fabric.push_back(new CFabrics1);
   fabric.push_back(new CFabrics2);
   
   ifstream in("dataset.txt"); string str;
	while(getline(in,str)){
		vector<ComplexNumber> vectmp;
		stringstream ss(str);
		for(ComplexNumber tmp;ss >> tmp;){
			vectmp.push_back(tmp);
		}

		CComplexVector *t=fabric[0]->Create();
		*t=vectmp;
		//t->input(outfile);
		v.push_back(t);
	}
   CComplexVector::ParallelTest(v);
   cout<<"Number of vectors = "<<v.size()<<endl;
   
    for(size_t i=0;(i<v.size());i++)  {
		delete v[i];
	}
	for(size_t i=0;(i<fabric.size());i++)  {
		delete fabric[i];
	}
   
  } catch(...) {cout << "error\n" <<endl;}
}

void test3()
{
try {
 CComplexVector2 a,b,c;
 a.resize(NUMBER);
 b.resize(NUMBER);
 c.resize(NUMBER);
 ComplexNumber var;
 //#pragma omp parallel for
 cout << "Setting Random values in a vector."<<endl;
 for(int i=0;i<NUMBER;i++)
 {
  float re = -100 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(200)));
  float im = -100 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(200)));
  ComplexNumber tmp(re,im);
  a[i]=tmp;
  re = -100 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(200)));
  im = -100 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(200)));
  ComplexNumber tmp1(re,im);
  b[i]=tmp1;
 }
 
 c=a+b;
 var=a*b;
 
} catch(...) {cout << "error\n" <<endl;}
}

int main(){
	test3();
	return 0;
}
