//#include "Potential.H"
#include "MetaMove.H"

MetaMove::MetaMove()
{
}

MetaMove::~MetaMove()
{
	pll.clear();
}

int 
MetaMove::setScoreImprovement(double aVal)
{
	scoreDelta=aVal;
	return 0;
}

int
MetaMove::setSrcMBScore(double aScore)
{
	mbscore=aScore;
	return 0;
}

int
MetaMove::setTargetMBScore(double aScore)
{
	targetMBScore=aScore;
	return 0;
}


int 
MetaMove::setSrcVertex(int vid)
{
	src=vid;
	return 0;
}


int 
MetaMove::setTargetVertex(int vid)
{	
	target=vid;
	return 0;
}

int 
MetaMove::setDestPot(Potential* apot)
{
	destPot=apot;
	return 0;
}


int 
MetaMove::getSrcVertex()
{
	return src;
}

int
MetaMove::getTargetVertex()
{
	return target;
}

double 
MetaMove::getSrcMBScore()
{
	return mbscore;
}

double
MetaMove::getTargetMBScore()
{
	return targetMBScore;
}


double 
MetaMove::getScoreImprovement()
{
	return scoreDelta;
}


Potential* 
MetaMove::getDestPot()
{
	return destPot;
}
