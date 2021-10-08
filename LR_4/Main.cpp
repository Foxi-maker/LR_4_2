#include "Header.h"

std::ofstream TwoDEq::stream;

int main()
{
	//TwoDEq ExFirst;
	//ExFirst.InitialCond(ExFirstGamma_1, ExFirstGamma_2, ExFirstGamma_3, ExFirstGamma_4, "ExFirst");
	//ExFirst.RightSide(ExFirstRight);
	//ExFirst.TridigAlg();

	TwoDEq ExSecond;
	ExSecond.InitialCond(ExSecondGamma_1, ExSecondGamma_2, ExSecondGamma_3, ExSecondGamma_4, "ExSecond");
	ExSecond.RightSide(ExSecondRight);
	ExSecond.TridigAlg();

	//TwoDEq ExThird;
	//ExThird.InitialCond(ExThirdGamma_1, ExThirdGamma_2, ExThirdGamma_3, ExThirdGamma_4, "ExThird");
	//ExThird.RightSide(ExThirdRight);
	//ExThird.TridigAlg();

	//TwoDEq ExNine;
	//ExNine.InitialCond(ExNineGamma_1, ExNineGamma_2, ExNineGamma_3, ExNineGamma_4, "ExNine");
	//ExNine.RightSide(ExNineRight);
	//ExNine.TridigAlg();
}