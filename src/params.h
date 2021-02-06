#pragma once

struct Params
{
    double bits;
    double FBB;
    double LPB;
    double M;
    double T;

    Params(double bits, double FBB, double LPB, double M, double T) :
        bits(bits), FBB(FBB), LPB(LPB), M(M), T(T)
    {
    }
};
