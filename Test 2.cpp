// Test 2.cpp : main project file.
#include "stdafx.h"
#include <cstdio>   //include printf()    
#include <cstdlib>  // rand() srand()   
#include <cmath>    //include cos(),sin()   
#include <ctime>    //include time()       
#include <limits> 
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream> //include cout cin    
using namespace std;
#define solution_size 30 
#define MOA_popsize 30                           
#define MOA_maxgen 1000 
int R=6;
double lowBound=-R/2,upperBound=R/2;  
#define rnd(low,uper) (( (double)rand() * ( uper - (low) ) ) / (double)RAND_MAX + (low))   
#define RUNTIMES 30 
#define L floor(pow(MOA_popsize,0.5))
vector<int> f;
vector<int> d; 

//////////////////////////////////////////////////////////////////////////////////////////
//                    readDataset Function
//  This function reads from dataset matrix and convertes to an array
//////////////////////////////////////////////////////////////////////////////////////////

void readDataset(){   // This function open .dat files and put them in an array. DAT files are our datasets
    ifstream numberfile; 
    int i; //needed to iterate between file and vector/array
    numberfile.open("d_nug30.dat");//opens file
    
    while(numberfile.good()){
        //will generally go untill end of file. it is better if you how long the file is
        numberfile >> i; //gets number 
        d.push_back(i);//adds number to end of the array
    }
    numberfile.close();
	numberfile.open("f_nug30.dat");//opens file
    
    while(numberfile.good()){
        //will generally go untill end of file. it is better if you how long the file is
        numberfile >> i; //gets number 
        f.push_back(i);//adds number to end of the array
    }
    numberfile.close();

	return;

}

//////////////////////////////////////////////////////////////////////////////////////////
//                    MergSort Class
//  This class sorts a vector and finds sesquences of every gene
//////////////////////////////////////////////////////////////////////////////////////////

class MergSort 
{
public:
	MergSort(vector<double> );
	void sort();
	vector<int> sequence();
private:
	int size; // vector size
	vector<double> data;
	vector<int> row;
	void sortSubVector(int,int);
	void merge(int,int,int,int);
};

MergSort::MergSort(vector<double> vectorInput)
{
	size=(int) vectorInput.size();
	data.resize(size);
	copy(vectorInput.begin(),vectorInput.end(),data.begin());  // copy vectorInput to data
	for(int i=0;i<size;i++){
		row.push_back(i); //sequence vector
	}
}

void MergSort::sort()
{
	sortSubVector(0,size-1); //recursively sort entire vector
	return;
}

void MergSort::sortSubVector(int low,int high)
{
	if ((high-low)>=1)
	{
		int middle1=(low+high)/2;
		int middle2=middle1+1;
		sortSubVector(low,middle1);
		sortSubVector(middle2,high);
		merge(low,middle1,middle2,high);
	}
	return;
}

void MergSort::merge(int left,int middle1,int middle2,int right)
{
	int leftIndex=left;
	int rightIndex=middle2;
	int combinedIndex=left;  //index into temporary working vector
	vector<double> combined(size);  //working vector
	vector<int> combinedRow(size);

	while(leftIndex<=middle1 && rightIndex<=right)
	{
		if(data[leftIndex]<=data[rightIndex]){
			combined[combinedIndex]=data[leftIndex];
			combinedRow[combinedIndex++]=row[leftIndex++];
	}
		else{
			combined[combinedIndex]=data[rightIndex];
			combinedRow[combinedIndex++]=row[rightIndex++];
		}
	}
	
	if(leftIndex==middle2)
	{
		while(rightIndex<=right){
			combined[combinedIndex]=data[rightIndex];
			combinedRow[combinedIndex++]=row[rightIndex++];
		}
	}
	else
	{
		while(leftIndex<=middle1){
			combined[combinedIndex]=data[leftIndex];
			combinedRow[combinedIndex++]=row[leftIndex++];
		}
	}
	
	for(int i=left;i<=right;i++){
		data[i]=combined[i];
		row[i]=combinedRow[i];
	}
	return;
}
vector<int> MergSort::sequence()
{
	vector<int> tempRow(size,0);
	for(int i=0;i<size;i++){
		tempRow[row[i]]=i;
	}
	return tempRow;
}
 
//////////////////////////////////////////////////////////////////////////////////////////
//                    Individual Struct
//  This struct defines a individual and calculates its fitness
//////////////////////////////////////////////////////////////////////////////////////////
struct Individual
{
	vector<double> indi;   
    vector<double> v;    
	vector<double> force; 
    double fitness; 
	void calcFitness();
	Individual():fitness(numeric_limits<double>::max()),force(solution_size),v(solution_size),indi(solution_size)
					
	{
		for(int j=0;j<solution_size;j++)   
		{          
			indi.at(j) = rnd(lowBound,upperBound);   
			v.at(j) = 0;
		}
	}


};

void Individual::calcFitness()   //In this function Problem is defined and calculated
{ 

	MergSort sortVector(indi);
	sortVector.sort();
	vector<int> sortRow(sortVector.sequence());
	double sumFitness=0;
	for(int i=0;i<solution_size;i++){
		for(int j=0;j<solution_size;j++){
			sumFitness+=(d[(sortRow[i]*solution_size)+sortRow[j]]*f[(i*solution_size)+j]);// f and row are 1D array but they are addressed in 2D style
		}
	}

	if (sumFitness<fitness)
		fitness=sumFitness;
	
	return;
}   
//////////////////////////////////////////////////////////////////////////////////////////
//                    Individual Magnetic Class
//  This class search through a vector of Individual struct by magnetic method
//////////////////////////////////////////////////////////////////////////////////////////

class Individualmagnetic{   //Using Individuals from Individual class for searching through the search space
public:   
	double gbest; 
    vector<Individual> MOA_pop;
	//Individual MOA_pop[MOA_popsize]; 
	vector<double> results;
	Individualmagnetic():gbest(numeric_limits<double>::max()),results(MOA_maxgen),MOA_pop(MOA_popsize)
	{
		
	}
//public:    
	void init();   
    void search(double,double);    
};  


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Individualmagnetic::init(){   //Getting fitness of every individual in population for more calculations

	double tempSum=0;
	vector<double> tempFitness(MOA_popsize);  
	double Max=(numeric_limits<double>::min());
	double Min=(numeric_limits<double>::max());
	for(int j=0;j<MOA_popsize;j++){ 
		MOA_pop[j].calcFitness(); 
		tempFitness.at(j)=(MOA_pop[j].fitness);
		if (tempFitness.at(j)>Max)
			Max=tempFitness.at(j);
		if (tempFitness.at(j)<Min)
			Min=tempFitness.at(j);
		//cout<<Min<<' ';
		tempSum+=MOA_pop[j].fitness;
	}
	//double Max = *max_element (tempFitness.begin(), tempFitness.end());// WARNINGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
	//double Min = *min_element (tempFitness.begin(), tempFitness.end());// min and max turn it to int. Min beconme 0
	//if ((tempSum-Min)==0){
       // for(int j=0;j<MOA_popsize;j++){ 
			//MOA_pop[j].fitness+=rnd(0,1);
		//}
	//}
	//cout<<Min<<' ' ;
	if(gbest>Min)
		gbest=Min;
	//for(int j=0;j<MOA_popsize;j++){
		//if(Max!=Min)
			//MOA_pop[j].fitness=(MOA_pop[j].fitness - Min)/(Max-Min); //Normalization
		//else
			//MOA_pop[j].fitness=(MOA_pop[j].fitness - Min)/abs(rand()+Max-Min);
	//}
	//cout<<"  "<<MOA_pop[3].fitness;
	return;
}   
          

void Individualmagnetic::search(double ro,double alpha)   //Magnetic search
{   
	
	int gen=0;   
	int neighbor[4];
	double distance[4];
	int row=0;
	int column=0;
	vector<double> help(solution_size);
	vector<double> mass(MOA_popsize);
	vector<double> tempF(solution_size);
	
	
	while(gen<MOA_maxgen) 
	{   
		
		init();    
		results.at(gen)=gbest;

		for(int ii=0;ii<MOA_popsize;ii++)   
		{  
			row=floor(ii/L);;
			column=fmod(ii,L);
			mass.at(ii)=((ro*MOA_pop[ii].fitness)+alpha);
			neighbor[0]=ii-1; //left hand side neighbor
			neighbor[1]=ii+1;//right hand side neighbor
			neighbor[2]=ii+L; //downer neighbor
			neighbor[3]=ii-L; //upper neighbor
			if (column==0) neighbor[0]=ii+L-1;
			if (column==L-1) neighbor[1]=ii-L+1;
			if (row==L-1 && (L*L)==MOA_popsize) neighbor[2]=ii-((L-1)*L);
			if (row==L && (L*L)!=MOA_popsize) neighbor[2]=ii-(L*L); 
			if (row==0 && (L*L)==MOA_popsize) neighbor[3]=ii+((L-1)*L);
			if (row==0 && (L*L)!=MOA_popsize) neighbor[3]=ii+(L*L);

			if (neighbor[0]>MOA_popsize) neighbor[0]=ii-1;
			if (neighbor[1]>MOA_popsize) neighbor[1]=ii-1;
			if (neighbor[2]>MOA_popsize) neighbor[2]=ii-1;
			if (neighbor[3]>MOA_popsize) neighbor[3]=ii-1;

			for(int j=0;j<solution_size;j++)
				tempF.at(j)=0;
			for(int j=0;j<4;j++){
				double sumMean=0;
				for(int k=0;k<solution_size;k++){
					help.at(k)=(MOA_pop[neighbor[j]].indi[k]-MOA_pop[ii].indi[k]);
					sumMean+=abs(help.at(k));
				}
				distance[j]=(sumMean/solution_size)/R;
				for(int k=0;k<solution_size;k++)
					tempF[k]+=(help[k]*MOA_pop[neighbor[j]].fitness)/(distance[j]+0.0001);


			}
			for(int k=0;k<solution_size;k++){
				MOA_pop[ii].force[k]=(tempF[k]*rnd(0,1));
				MOA_pop[ii].v[k]=(MOA_pop[ii].force[k]/mass[ii]);
				MOA_pop[ii].indi[k]=MOA_pop[ii].indi[k]+MOA_pop[ii].v[k];
				MOA_pop[ii].indi[k]=(MOA_pop[ii].indi[k]>-R/2 ?MOA_pop[ii].indi[k]:-R/2); 
				MOA_pop[ii].indi[k]=(MOA_pop[ii].indi[k]<R/2 ?MOA_pop[ii].indi[k]:R/2); 
			}
		} 
		gen++;         
	}   
	return;
}   
//------------------------------------------------------------------------------------
//                 Main
//------------------------------------------------------------------------------------

int main (array<System::String ^> ^args) {

	int roSize=5;
	int alphaSize=5;
	double roArray[]={0.1,0.2,0.5,0.7,0.9};
	double alphaArray[]={0.1,0.2,0.5,0.7,0.9};
	vector<double> ro(roArray,roArray+roSize);
	vector<double> alpha(alphaArray,alphaArray+alphaSize);
	double meanResult=0;
    readDataset();
    srand((unsigned)time(NULL));
	filebuf buffer;
    ostream output(&buffer);
    istream input(&buffer);
	filebuf meanBuffer;
    ostream meanOutput(&meanBuffer);
	istream meanInput(&meanBuffer);
    buffer.open ("results.dat", ios::in | ios::out | ios::trunc);
	meanBuffer.open ("meanResults.dat", ios::in | ios::out | ios::trunc);
	for(int k=0;k<roSize;k++)
		for(int w=0;w<alphaSize;w++){
			meanResult=0;
			for (int i = 0; i < RUNTIMES; i++)     
			{   
				cout<<"Run :"<<i+1<<endl;
				Individualmagnetic p; 
				p.search(ro[k],alpha[w]);   
				meanResult+=p.results.at(MOA_maxgen-1);
				output << "Run number:"<<i<<"  "<< "ro:"<<ro[k]<<"  "<< "alpha:"<<alpha[w]<<" = "<<p.gbest<< endl ;
				output << "Results:"<<i<<"  ";

				for (int j=0;j<MOA_maxgen;j++)
					output <<"  "<<p.results.at(j);  // This Part can be commented 
				output <<endl;
			}
			meanOutput <<"ro:"<<ro[k]<<"  "<< "alpha:"<<alpha[w]<<" = "<<(meanResult/RUNTIMES)<< endl ;
		}
    system("pause");  
	return 0;
}