#include "TwoDEq.h"

TwoDEq::TwoDEq()
{
	U.resize(NumVerticalStep);
	for (auto& u : U)
		u.resize(NumHorizontalStep);

	ProjRightSide.resize(NumVerticalStep);
	for (auto& pR : ProjRightSide)
		pR.resize(NumHorizontalStep);
}

void TwoDEq::InitialCond(std::function<double(double)> Gamma_1, std::function<double(double)> Gamma_2, std::function<double(double)> Gamma_3, std::function<double(double)> Gamma_4, std::string streamName)
{
	//Initial conditions
	for (int ExternIter = 0; ExternIter < NumHorizontalStep; ExternIter++)
	{
		U[0][ExternIter] = Gamma_3(ExternIter * HorizontalStep);
		U.back()[ExternIter] = Gamma_2(ExternIter * HorizontalStep);
	}
	for (int ExternIter = 0; ExternIter < NumVerticalStep; ExternIter++)
	{
		U[ExternIter][0] = Gamma_1(ExternIter * VerticalStep);
		U[ExternIter].back() = Gamma_4(ExternIter * VerticalStep);
	}

	for (int RowIter = 1; RowIter < NumVerticalStep - 1; RowIter++)
		for (int ColIter = 1; ColIter < NumHorizontalStep - 1; ColIter++)
			U[RowIter][ColIter] = 0.;

	stream.open("OUTPUT\\" + streamName + ".txt");
	if (!stream.is_open())
	{
		std::cout << "Error! File for writing was not open!";
		exit(1);
	}
}

void TwoDEq::RightSide(std::function<double(double, double)> RightS)
{
	for (int ExternIter = 0; ExternIter < NumVerticalStep - 2; ExternIter++)
		for (int InterIter = 0; InterIter < NumHorizontalStep - 2; InterIter++)
			ProjRightSide[ExternIter][InterIter] = RightS((ExternIter + 1) * VerticalStep, (InterIter + 1) * HorizontalStep);
}

void TwoDEq::TridigAlg()
{
	std::vector<double> alpha, betta;

	for (int timeIter = 0; timeIter < 2 * NumTimeStep; timeIter++)
	{
		//F components
		std::vector<std::vector<double>> F(NumVerticalStep - 2, std::vector<double>(NumHorizontalStep - 2));
		for (int RowIter = 0; RowIter < NumVerticalStep - 2; RowIter++)
			for (int ColIter = 0; ColIter < NumHorizontalStep - 2; ColIter++)
				F[RowIter][ColIter] = -(2. / TimeStep * U[RowIter + 1][ColIter + 1] + 1. / (HorizontalStep * HorizontalStep) * (U[RowIter + 1][ColIter + 2] - 2. * U[RowIter + 1][ColIter + 1] + U[RowIter + 1][ColIter])) + ProjRightSide[RowIter][ColIter];

		double MainDiagonal = -2. * (1. / (VerticalStep * VerticalStep) + 1. / TimeStep),
			SideDiagonal = 1. / (VerticalStep * VerticalStep);

		for (int HorizIter = 0; HorizIter < NumHorizontalStep - 2; HorizIter++)
		{
			F[0][HorizIter] -= SideDiagonal;
			F.back()[HorizIter] -= SideDiagonal;
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

		alpha.clear();
		betta.clear();

		//F_tilda components
		for (int RowIter = 0; RowIter < NumVerticalStep - 2; RowIter++)
			for (int ColIter = 0; ColIter < NumHorizontalStep - 2; ColIter++)
				F[RowIter][ColIter] = -(2. / TimeStep * U[RowIter + 1][ColIter + 1] + 1. / (VerticalStep * VerticalStep) * (U[RowIter + 2][ColIter + 1] - 2. * U[RowIter + 1][ColIter + 1] + U[RowIter][ColIter + 1])) + ProjRightSide[RowIter][ColIter];

		MainDiagonal = -2. * (1. / (HorizontalStep * HorizontalStep) + 1. / TimeStep);
		SideDiagonal = 1. / (HorizontalStep * HorizontalStep);

		for (int VertIter = 0; VertIter < NumVerticalStep - 2; VertIter++)
		{
			F[VertIter][0] -= SideDiagonal;
			F[VertIter].back() -= SideDiagonal;
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

		//std::cout << "U\n";
		//for (const auto& UExtern : U)
		//{
		//	for (const auto& UInter : UExtern)
		//		std::cout << UInter << " ";
		//	std::cout << "\n";
		//}
		alpha.clear();
		betta.clear();
	}

	std::cout << "U\n";
	for (const auto& UExtern : U)
	{
		for (const auto& UInter : UExtern)
			std::cout << UInter << " ";
		std::cout << "\n";
	}

	for (const auto& UExtern : U)
	{
		for (const auto& UInter : UExtern)
			stream << UInter << " ";
		stream << "\n";
	}
}

TwoDEq::~TwoDEq()
{
	stream.close();
}
