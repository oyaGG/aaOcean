// ContainArray.h
// Template classes for multidimensional arrays:
// Vector (1-D), Matrix (2-D), and Tensor (3-D) in namespace ContainArray
// Richard J. Wagner  v1.1  14 March 2005  wagnerr@umich.edu

// Copyright (c) 2005 Richard J. Wagner
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// It would be nice to send a note to the author at wagnerr@umich.edu if you
// find this software useful.

// Typical usage
// -------------
// 
// using namespace ContainArray;
// 
// Vector<int> vecint(3);             // vector with 3 rows
// Vector<int> vecone( 3, 1 );        // vector with 3 rows, filled with 1's
// Matrix<double> matdbl( 3, 2 );     // matrix with 3 rows and 2 columns
// Tensor<string> tenstr( 3, 4, 2 );  // tensor with 3 rows, 4 cols, 2 lays
// 
// vecint[1] = 7;              // assign element of vector
// matdbl[0][1] = 7.0;         // assign element of matrix
// tenstr[1][3][1] = "seven";  // assign element of tensor
// tenstr(1,3,1) = "eight";    // assign element of tensor again
// vecint.at(2) = 77;          // assign element with range checking
// matdbl = 77.7;              // assign all elements
// 
// vecint.swap(vecone);  // swap vectors of same type
// 
// std::cout << vecint << std::endl;  // write vector to output
// std::cin >> vecint;                // read vector from input
// 
// See file example.cpp for more examples.

#ifndef CONTAINARRAY_H
#define CONTAINARRAY_H

//#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace ContainArray {

/*----------------------------------------------------------------------------*/
/*                                   Vector                                   */
/*----------------------------------------------------------------------------*/

template<class T=int>
class Vector {
private:
	// Data
	size_t  mySize;
	T      *myArray;
	
public:
	// Constructors
	Vector( size_t rows = 0 );
	Vector( size_t rows, const T& value );
	Vector( const Vector<T> & orig );
	~Vector();
	
	// Inspectors
	size_t size() const { return mySize; }
	size_t rows() const { return size(); }
	
	// Manipulators
	Vector& operator=( const Vector<T> & v )
	{
		if( v.size() == size() )
		{
			for( size_t i = 0; i < size(); ++i ) (*this)[i] = v[i];
		}
		else
		{
			Vector<T> cp(v);
			(*this).swap(cp);
		}
		return *this;
	}
	Vector& operator=( const T& value )
	{
		for( size_t i = 0; i < size(); ++i ) (*this)[i] = value;
		return *this;
	}
	void resize( size_t rows, const T& value = T() )
	{
		if( rows == this->rows() ) return;
		Vector<T> v( rows, value );
		for( size_t i = 0; i < rows && i < this->rows(); ++i )
			v[i] = (*this)[i];
		swap(v);
	}
	void swap( Vector<T> & v )
	{
		std::swap( mySize, v.mySize );
		std::swap( myArray, v.myArray );
	}
	
	// Accessors
	T& operator[]( size_t i ) { return myArray[i]; }
	const T& operator[]( size_t i ) const { return myArray[i]; }
	T& operator()( size_t i ) { return (*this)[i]; }
	const T& operator()( size_t i ) const { return (*this)[i]; }
	T& at( size_t i )
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	const T& at( size_t i ) const
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	
	// Comparators
	bool operator!=( const Vector<T> & v ) const
	{
		if( size() != v.size() ) return true;
		for( size_t i = 0; i < size(); ++i )
			if( (*this)[i] != v[i] ) return true;
		return false;
	}
	bool operator==( const Vector<T> & v ) const
	{
		return !( *this != v );
	}
};


template<class T>
Vector<T>::Vector( size_t rows )
	: mySize(rows)
{
	myArray = new T[mySize];
}


template<class T>
Vector<T>::Vector( size_t rows, const T& value )
	: mySize(rows)
{
	myArray = new T[mySize];
	for( size_t i = 0; i < mySize; ++i ) myArray[i] = value;
}


template<class T>
Vector<T>::Vector( const Vector<T> & orig )
	: mySize(orig.mySize)
{
	myArray = new T[mySize];
	for( size_t i = 0; i < mySize; ++i ) myArray[i] = orig.myArray[i];
}


template<class T>
Vector<T>::~Vector()
{
	delete [] myArray;
}

/*----------------------------------------------------------------------------*/
/*                                   Matrix                                   */
/*----------------------------------------------------------------------------*/

template<class T=int>
class Matrix : public Vector< Vector<T> > {
public:
	// Constructors
	Matrix( size_t rows = 0, size_t cols = 1 );
	Matrix( size_t rows, size_t cols, const T& value );
	Matrix( const Matrix<T> & orig );
	
	// Inspectors
	size_t cols() const { return size() > 0 ? (*this)[0].size() : 0; }
	
	// Manipulators
	Matrix& operator=( const Matrix<T> & m )
	{
		if( m.size() == size() )
		{
			for( size_t i = 0; i < size(); ++i ) (*this)[i] = m[i];
		}
		else
		{
			Matrix<T> cp(m);
			(*this).swap(cp);
		}
		return *this;
	}
	Matrix& operator=( const T& value )
		{ for( size_t i = 0; i < size(); ++i ) (*this)[i] = value;
		  return *this; }
	void resize( size_t rows, size_t cols, const T& value = T() )
	{
		if( rows == this->rows() && cols == this->cols() ) return;
		Matrix<T> m( rows, cols, value );
		for( size_t i = 0; i < rows && i < this->rows(); ++i )
			for( size_t j = 0; j < cols && j < this->cols(); ++j )
				m[i][j] = (*this)[i][j];
		swap(m);
	}
	
	// Additional accessors
	Vector<T>& at( size_t i )
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	const Vector<T>& at( size_t i ) const
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	T& at( size_t i, size_t j )
		{ return at(i).at(j); }
	const T& at( size_t i, size_t j ) const
		{ return at(i).at(j); }
	Vector<T>& operator()( size_t i )
		{ return (*this)[i]; }
	const Vector<T>& operator()( size_t i ) const
		{ return (*this)[i]; }
	T& operator()( size_t i, size_t j )
		{ return (*this)[i][j]; }
	const T& operator()( size_t i, size_t j ) const
		{ return (*this)[i][j]; }
};


template<class T>
Matrix<T>::Matrix( size_t rows, size_t cols )
	: Vector< Vector<T> >( rows, Vector<T>(cols) ) {}


template<class T>
Matrix<T>::Matrix( size_t rows, size_t cols, const T& value )
	: Vector< Vector<T> >( rows, Vector<T>(cols,value) ) {}


template<class T>
Matrix<T>::Matrix( const Matrix<T> & orig )
	: Vector< Vector<T> >(orig) {}

/*----------------------------------------------------------------------------*/
/*                                   Tensor                                   */
/*----------------------------------------------------------------------------*/

template<class T=int>
class Tensor : public Vector< Matrix<T> > {
public:
	// Constructors
	Tensor( size_t rows = 0, size_t cols = 1, size_t lays = 1 );
	Tensor( size_t rows, size_t cols, size_t lays, const T& value );
	Tensor( const Tensor<T> & orig );
	
	// Inspectors
	size_t cols() const { return size() > 0 ? (*this)[0].size() : 0; }
	size_t lays() const { return size() > 0 ? (*this)[0].cols() : 0; }
	
	// Manipulators
	Tensor& operator=( const Tensor<T> & t )
	{
		if( t.size() == size() )
		{
			for( size_t i = 0; i < size(); ++i ) (*this)[i] = t[i];
		}
		else
		{
			Tensor<T> cp(t);
			(*this).swap(cp);
		}
		return *this;
	}
	Tensor& operator=( const T& value )
		{ for( size_t i = 0; i < size(); ++i ) (*this)[i] = value;
		  return *this; }
	void resize( size_t rows, size_t cols, size_t lays, const T& value = T() )
	{
		if( rows == this->rows() && cols == this->cols() && lays == this->lays() )
			return;
		Tensor<T> t( rows, cols, lays, value );
		for( size_t i = 0; i < rows && i < this->rows(); ++i )
			for( size_t j = 0; j < cols && j < this->cols(); ++j )
				for( size_t k = 0; k < lays && k < this->lays(); ++k )
					t[i][j][k] = (*this)[i][j][k];
		swap(t);
	}
	
	// Additional accessors
	Matrix<T>& at( size_t i )
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	const Matrix<T>& at( size_t i ) const
	{
		if( i >= size() ) throw std::out_of_range( "index out of range" );
		return (*this)[i];
	}
	T& at( size_t i, size_t j, size_t k )
		{ return at(i).at(j).at(k); }
	const T& at( size_t i, size_t j, size_t k ) const
		{ return at(i).at(j).at(k); }
	Matrix<T>& operator()( size_t i )
		{ return (*this)[i]; }
	const Matrix<T>& operator()( size_t i ) const
		{ return (*this)[i]; }
	T& operator()( size_t i, size_t j, size_t k )
		{ return (*this)[i][j][k]; }
	const T& operator()( size_t i, size_t j, size_t k ) const
		{ return (*this)[i][j][k]; }
};


template<class T>
Tensor<T>::Tensor( size_t rows, size_t cols, size_t lays )
	: Vector< Matrix<T> >( rows, Matrix<T>(cols,lays) ) {}


template<class T>
Tensor<T>::Tensor( size_t rows, size_t cols, size_t lays, const T& value )
	: Vector< Matrix<T> >( rows, Matrix<T>(cols,lays,value) ) {}


template<class T>
Tensor<T>::Tensor( const Tensor<T> & orig )
	: Vector< Matrix<T> >(orig) {}

}  // namespace ContainArray

/*----------------------------------------------------------------------------*/
/*                              Input and output                              */
/*----------------------------------------------------------------------------*/
//
//template<class T>
//std::ostream& operator<<( std::ostream& os, const ContainArray::Vector<T> & v )
//{
//	// Output a Vector to os
//	os << v.rows();
//	for( size_t i = 0; i < v.rows(); ++i ) os << '\t' << v[i];
//	return os;
//}
//
//
//template<class T>
//std::istream& operator>>( std::istream& is, ContainArray::Vector<T> & v )
//{
//	// Input a Vector from is
//	size_t rows;
//	is >> rows;
//	if( rows == v.rows() )
//	{
//		for( size_t i = 0; i < rows; ++i )
//			is >> v[i];
//		return is;
//	}
//	ContainArray::Vector<T> vin(rows);
//	for( size_t i = 0; i < rows; ++i )
//		is >> vin[i];
//	vin.swap(v);
//	return is;
//}
//
//
//template<class T>
//std::ostream& operator<<( std::ostream& os, const ContainArray::Matrix<T> & m )
//{
//	// Output a Matrix to os
//	size_t rows = m.rows();
//	size_t cols = m.cols();
//	os << rows << ' ' << cols << '\n';
//	for( size_t i = 0; i < rows; ++i )
//	{
//		for( size_t j = 0; j < cols; ++j )
//		{
//			os << m[i][j];
//			if( j != cols-1 ) os << '\t';
//		}
//		if( i != rows-1 ) os << '\n';
//	}
//	return os;
//}
//
//
//template<class T>
//std::istream& operator>>( std::istream& is, ContainArray::Matrix<T> & m )
//{
//	// Input a Matrix from is
//	size_t rows;
//	size_t cols;
//	is >> rows >> cols;
//	if( rows == m.rows() && cols == m.cols() )
//	{
//		for( size_t i = 0; i < rows; ++i )
//			for( size_t j = 0; j < cols; ++j )
//				is >> m[i][j];
//		return is;		
//	}
//	ContainArray::Matrix<T> min( rows, cols );
//	for( size_t i = 0; i < rows; ++i )
//		for( size_t j = 0; j < cols; ++j )
//			is >> min[i][j];
//	min.swap(m);
//	return is;
//}
//
//
//template<class T>
//std::ostream& operator<<( std::ostream& os, const ContainArray::Tensor<T> & t )
//{
//	// Output a Tensor to os
//	size_t rows = t.rows();
//	size_t cols = t.cols();
//	size_t lays = t.lays();
//	os << rows << ' ' << cols << ' ' << lays << '\n';
//	for( size_t k = 0; k < lays; ++k )
//	{
//		for( size_t i = 0; i < rows; ++i )
//		{
//			for( size_t kindent = 0; kindent < k; ++kindent ) os << "\t";
//			for( size_t j = 0; j < cols; ++j )
//			{
//				os << t[i][j][k];
//				if( j != cols-1 ) os << '\t';
//			}
//			if( i != rows-1 || k != lays-1 ) os << '\n';
//		}
//	}
//	return os;
//}
//
//
//template<class T>
//std::istream& operator>>( std::istream& is, ContainArray::Tensor<T> & t )
//{
//	// Input a Tensor from is
//	size_t rows;
//	size_t cols;
//	size_t lays;
//	is >> rows >> cols >> lays;
//	if( rows == t.rows() && cols == t.cols() && lays == t.lays() )
//	{
//		for( size_t k = 0; k < lays; ++k )
//			for( size_t i = 0; i < rows; ++i )
//				for( size_t j = 0; j < cols; ++j )
//					is >> t[i][j][k];
//		return is;
//	}
//	ContainArray::Tensor<T> tin( rows, cols, lays );
//	for( size_t k = 0; k < lays; ++k )
//		for( size_t i = 0; i < rows; ++i )
//			for( size_t j = 0; j < cols; ++j )
//				is >> tin[i][j][k];
//	tin.swap(t);
//	return is;
//}

#endif  // CONTAINARRAY_H

// Change log:
//
// v1.0  13 March 2005
//   + First release
// 
// v1.1  14 March 2005
//   + Added missing "#include <algorithm>"
//   + Removed debug statements
//   + Declared comparators as const
//   + Defined single-index at() and operator() for Matrix and Tensor
//   + Moved operator<< and operator>> outside of namespace ContainArray
