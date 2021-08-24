#include <iostream>
#include <vector>
#include "SLAE.h"

int main()
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
	double TimeEnd = 10.;
	int NumTimeStep = 5;
	double TimeStep = TimeEnd / NumTimeStep;				//tau

	//Projecting the right side of equation
	//std::vector<std::vector<double>> ProjRightSide;

	//Temperature distribution in space
	std::vector<std::vector<double>> U(NumVerticalStep, std::vector<double>(NumHorizontalStep));

	//Initial conditions
	for (int ExternIter = 0; ExternIter < NumHorizontalStep; ExternIter++)
	{
		U[0][ExternIter] = 1.;
		U.back()[ExternIter] = 1.;
	}
	for (int ExternIter = 0; ExternIter < NumVerticalStep; ExternIter++)
	{
		U[ExternIter][0] = 1.;
		U[ExternIter].back() = 1.;
	}

	for (int RowIter = 1; RowIter < NumVerticalStep - 1; RowIter++)
		for (int ColIter = 1; ColIter < NumHorizontalStep - 1; ColIter++)
			U[RowIter][ColIter] = 0.;

	std::cout << "U_0\n";
	for (const auto& UExtern : U)
	{
		for (const auto& UInter : UExtern)
			std::cout << UInter << " ";
		std::cout << "\n";
	}

	std::vector<double> alpha, betta;

	for (int timeIter = 0; timeIter < 2 * NumTimeStep; timeIter++)
	{
		//F components
		std::vector<std::vector<double>> F(NumVerticalStep - 2, std::vector<double>(NumHorizontalStep - 2));
		for (int RowIter = 0; RowIter < NumVerticalStep - 2; RowIter++)
			for (int ColIter = 0; ColIter < NumHorizontalStep - 2; ColIter++)
				F[RowIter][ColIter] = -(2. / TimeStep * U[RowIter + 1][ColIter + 1] + 1. / (HorizontalStep * HorizontalStep) * (U[RowIter + 1][ColIter + 2] - 2. * U[RowIter + 1][ColIter + 1] + U[RowIter + 1][ColIter]));

		//std::cout << "F\n";
		//for (const auto& fExtern : F)
		//{
		//	for (const auto& fInter : fExtern)
		//		std::cout << fInter << " ";
		//	std::cout << "\n";
		//}

		double MainDiagonal = -2. * (1. / (VerticalStep * VerticalStep) + 1. / TimeStep),
			SideDiagonal = 1. / (VerticalStep * VerticalStep);

		for (int HorizIter = 0; HorizIter < NumHorizontalStep - 2; HorizIter++)
		{
			F[0][HorizIter] += SideDiagonal;
			F.back()[HorizIter] += SideDiagonal;
		}

		//Along the X1 direction (Vertical)
		for (int HorizIter = 0; HorizIter < NumHorizontalStep - 2; HorizIter++)
		{
			alpha.push_back(-SideDiagonal / MainDiagonal);
			betta.push_back(F[0][HorizIter] / MainDiagonal);
			for (int VertIter = 1; VertIter < NumVerticalStep - 3; VertIter++)
			{
				double denominator = MainDiagonal + SideDiagonal * alpha.back();
				alpha.push_back(-SideDiagonal / denominator);
				betta.push_back((-SideDiagonal * betta.back() + F[VertIter][HorizIter]) / denominator);
			}
			U[NumVerticalStep - 2][HorizIter + 1] = (-SideDiagonal * betta.back() + F[NumVerticalStep - 3][HorizIter]) / (MainDiagonal + SideDiagonal * alpha.back());

			for (int VertIter = NumVerticalStep - 3; VertIter > 0; VertIter--)
				U[VertIter][HorizIter + 1] = alpha[VertIter - 1] * U[VertIter + 1][HorizIter + 1] + betta[VertIter - 1];
		}

		//F_tilda components
		for (int RowIter = 0; RowIter < NumVerticalStep - 2; RowIter++)
			for (int ColIter = 0; ColIter < NumHorizontalStep - 2; ColIter++)
				F[RowIter][ColIter] = -(2. / TimeStep * U[RowIter + 1][ColIter + 1] + 1. / (VerticalStep * VerticalStep) * (U[RowIter + 2][ColIter + 1] - 2. * U[RowIter + 1][ColIter + 1] + U[RowIter][ColIter + 1]));

		//std::cout << "F_tilda\n";
		//for (const auto& fExtern : F)
		//{
		//	for (const auto& fInter : fExtern)
		//		std::cout << fInter << " ";
		//	std::cout << "\n";
		//}

		MainDiagonal = -2. * (1. / (HorizontalStep * HorizontalStep) + 1. / TimeStep);
		SideDiagonal = 1. / (HorizontalStep * HorizontalStep);

		for (int VertIter = 0; VertIter < NumVerticalStep - 2; VertIter++)
		{
			F[VertIter][0] += SideDiagonal;
			F[VertIter].back() += SideDiagonal;
		}

		//Along the X2 direction (Horizontal)
		for (int VertIter = 0; VertIter < NumVerticalStep - 2; VertIter++)
		{
			alpha.push_back(-SideDiagonal / MainDiagonal);
			betta.push_back(F[VertIter][0] / MainDiagonal);
			for (int HorizIter = 1; HorizIter < NumHorizontalStep - 3; HorizIter++)
			{
				double denominator = MainDiagonal + SideDiagonal * alpha.back();
				alpha.push_back(-SideDiagonal / denominator);
				betta.push_back((-SideDiagonal * betta.back() + F[VertIter][HorizIter]) / denominator);
			}
			U[VertIter + 1][NumHorizontalStep - 2] = (-SideDiagonal * betta.back() + F[VertIter][NumHorizontalStep - 3]) / (MainDiagonal + SideDiagonal * alpha.back());

			for (int HorizIter = NumHorizontalStep - 3; HorizIter > 0; HorizIter--)
				U[VertIter + 1][HorizIter] = alpha[HorizIter - 1] * U[VertIter + 1][HorizIter + 1] + betta[HorizIter - 1];
		}

		std::cout << "U\n";
		for (const auto& UExtern : U)
		{
			for (const auto& UInter : UExtern)
				std::cout << UInter << " ";
			std::cout << "\n";
		}

	}

	//Tridiagonal Matrix Algorithm
	//std::vector<std::vector<double>> A = { {3,2,0,0},{1,3,2,0},{0,1,3,2},{0,0,1,3} };
	//std::vector<double> b = { 1,2,3,4 };

	//std::vector<double> alpha, betta, solve;

	//alpha.push_back(-A[0][1] / A[0][0]);
	//betta.push_back(b[0] / A[0][0]);
	//for (int Iter = 1; Iter < (A.size() - 1); Iter++)
	//{
	//	double denominator = A[Iter][Iter] + A[Iter][Iter - 1] * alpha.back();
	//	alpha.push_back(-A[Iter][Iter + 1] / denominator);
	//	betta.push_back((-A[Iter][Iter - 1] * betta.back() + b[Iter]) / denominator);
	//}

	//solve.push_back((-A.back()[A.size() - 2] * betta.back()+b.back()) / (A.back().back() + A.back()[A.size() - 2] * alpha.back()));
	//for (int Iter = alpha.size()-1; Iter >= 0; Iter--)
	//	solve.push_back(alpha[Iter] * solve.back() + betta[Iter]);

	//for (const auto& s : solve)
	//	std::cout << s << " ";
}