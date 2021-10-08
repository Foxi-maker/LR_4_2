#pragma once
#include "Header.h"

class TwoDEq
{
	//Space grid
	//X2 - Horizontal
	//X1 - Vertical 
	double HorizontalEnd = 1., VerticalEnd = 1.;
	int NumHorizontalStep = 5;
	int NumVerticalStep = 5;
	double HorizontalStep = HorizontalEnd / NumHorizontalStep;				//h_2
	double VerticalStep = VerticalEnd / NumVerticalStep;					//h_1

	//Time grid
	double TimeEnd = 1.;
	int NumTimeStep = 100;
	double TimeStep = TimeEnd / NumTimeStep;				//tau

	//Projecting the right side of equation
	std::vector<std::vector<double>> ProjRightSide;

	//Temperature distribution in space
	std::vector<std::vector<double>> U;

	static std::ofstream stream;
public:

	TwoDEq();
	
	void InitialCond(std::function<double(double)>, std::function<double(double)>, std::function<double(double)>, std::function<double(double)>, std::string);

	void RightSide(std::function<double(double,double)>);

	//Tridiagonal Matrix Algorithm
	void TridigAlg();

	~TwoDEq();
};

