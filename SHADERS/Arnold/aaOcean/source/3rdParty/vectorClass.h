
class mVector
{
    private:
        // NO PRIVATE DATA
    public:
        double i, j, k;
        mVector();
        mVector(double,double);
        mVector(double,double,double);
        mVector operator + (mVector);
        mVector operator - (mVector);
        mVector operator = (mVector);
        mVector operator / (mVector);
        double  operator % (mVector);
        double  operator ^ (mVector);
		mVector operator * (mVector);
        mVector operator * (double);
        mVector operator / (double);
        double mod();
        mVector unit();

};

// INITIALIZE EMPTY mVector
mVector::mVector()
{
    i=0;
    j=0;
    k=0;
}

// INITIALIZE 2D mVector
mVector::mVector(double a, double b)
{
    i=a;
    j=b;
    k=0;
}

// INITIALIZE 3D mVector
mVector::mVector(double a, double b, double c)
{
    i=a;
    j=b;
    k=c;
}

// ADD TWO mVectorS TOGETHER
mVector mVector::operator + (mVector param)
{
    mVector c;
    c.i = i + param.i;
    c.j = j + param.j;
    c.k = k + param.k;
    return (c);
}

// SUBTRACT ONE mVector FROM ANOTHER
mVector mVector::operator - (mVector param)
{
    mVector c;
    c.i = i - param.i;
    c.j = j - param.j;
    c.k = k - param.k;
    return (c);
}

// SET A mVector TO ANOTHER mVector
mVector mVector::operator = (mVector param)
{
    i=param.i;
    j=param.j;
    k=param.k;
    return *this;
}



// CROSS PRODUCT
mVector mVector::operator / (mVector param)
{
    mVector c;
    c.i=j*param.k-k*param.j;
    c.j=k*param.i-i*param.k;
    c.k=i*param.j-j*param.i;
    return (c);
}

// DOT PRODUCT
double mVector::operator % (mVector param)
{
    return (i*param.i + j*param.j + k*param.k);
}

// ANGLE BETWEEN IN RADIANS
double mVector::operator ^ (mVector b)
{
    mVector a=*this;
    return (acos((a%b)/(a.mod()*b.mod())));
}

// MULTIPLY BY REAL
mVector mVector::operator * (double b)
{
    mVector c;
    c.i=i*b;
    c.j=j*b;
    c.k=k*b;
    return (c);
}

// DIVIDE BY REAL
mVector mVector::operator / (double b)
{
    mVector c;
    c.i=i/b;
    c.j=j/b;
    c.k=k/b;
    return (c);
}

// STRETCH mVectorS
mVector mVector::operator * (mVector b)
{
    mVector c;
    c.i=i*b.i;
    c.j=j*b.j;
    c.k=k*b.k;
    return (c);
}

// MODULUS
double mVector::mod()
{
    return sqrt(pow(i,2)+pow(j,2)+pow(k,2));
}

// UNIT mVector
mVector mVector::unit()
{
    mVector c;
    c.i=i/this->mod();
    c.j=j/this->mod();
    c.k=k/this->mod();
    return c;
}
