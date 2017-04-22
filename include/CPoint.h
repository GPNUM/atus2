
#ifndef __class_CPoint__
#define __class_CPoint__

#include <string>
#include <initializer_list>
#include <iostream>
#include <vector>

/**
  * Template class for n-dimensional objects, for example position
  * and momentum vectors.
  *
  * Multiple operators are overloaded such as addition, multiplication and so on.
  */
template <int dim>
class CPoint
{
public:
  /**
    *  Create CPoint initialised with 0
    */
  CPoint()
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] = 0;
  };

  /**
    *  Create CPoint initialised with doubles from a list
    */
  CPoint( std::initializer_list<double> list )
  {
    if ( list.size() != dim ) throw std::string("Error: initializer_list size mismatch\n");
    int i=0;
    for ( auto it : list )
    {
      this->m_data[i++] = it;
    }
  }

  /**
    * Overload = between CPoint and double. Set all elements of CPoint to the value of double
    */
  CPoint &operator=(const double s)
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] = s;
    return *this;
  }

  /** Set CPoint to Values of std::vector */
  CPoint &operator=(const std::vector<double> &vec)
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] = vec[i];
    return *this;
  }

  /**
    * Get value of element CPoint[index]
    */
  double &operator[](int index)
  {
    if ( index<0 || index>dim ) throw std::string("Error: Index out of bound\n");
    return m_data[index];
  }

  /**
    * Overload += operator for two CPoints
    */
  CPoint &operator+=(const CPoint &right)
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] += right.m_data[i];
    return *this;
  }

  /**
    * Overload -= operator for two CPoints
    */
  CPoint &operator-=(const CPoint &right)
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] -= right.m_data[i];
    return *this;
  }

  /**
    * Overload *= between CPoint and double
    */
  CPoint &operator*=( const double s )
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] *= s;
    return *this;
  }

  /**
    * Overload /= between CPoint and double
    */
  CPoint &operator/=( const double s )
  {
    for ( int i=0; i<dim; i++ )
      this->m_data[i] /= s;
    return *this;
  }

  /**
    * Overload + for two CPoints. Element-wise addition
    */
  CPoint operator+(CPoint &right)
  {
    CPoint retval;
    for ( int i=0; i<dim; i++ )
      retval[i] = this->m_data[i] + right[i];
    return retval;
  }

  /**
    * Overload - for two CPoints. Element-wise subtraction
    */
  CPoint operator-(CPoint &right)
  {
    CPoint retval;
    for ( int i=0; i<dim; i++ )
      retval[i] = this->m_data[i] - right[i];
    return retval;
  }

  /**
    * Vector multiplication between CPoints
    */
  CPoint scale( CPoint &b )
  {
    CPoint retval;
    for ( int i=0; i<dim; i++ )
      retval[i] = this->m_data[i] * b[i];
    return retval;
  }

  /**
    * Overload * for two CPoints. Element-wise multiplication
    */
  double operator*(CPoint &right)
  {
    double retval=0;
    for ( int i=0; i<dim; i++ )
      retval += this->m_data[i] * right[i];
    return retval;
  }

  bool operator==( const double s)
  {
    for ( int i=0; i<dim; i++ )
      if ( this->m_data[i] != s ) return false;
    return true;
  }

  /** Overload * for CPoint and Double */
  CPoint operator*(double right)
  {
    CPoint retval;
    for ( int i= 0; i<dim; i++ )
      retval = this->m_data[i] * right;
    return retval;
  }

  /**
    * Store elements of CPoint in std::ofstream variable os
    */
  void Dump( std::ostream &os )
  {
    //os << "(";
    for ( int i=0; i<dim; i++ )
    {
      //os << m_data[i] << ((i==dim-1)? ")" : ", ");
      os << m_data[i] << ((i==dim-1)? "" : ",");
    }
  }

protected:
  ///Elements of CPoint stored in an array of doubles of dimension dim
  double m_data[dim];
};

/** \relates CPoint
  * Overload operator << between ofstream and CPoint
  */
template<int dim>
std::ostream &operator<<(std::ostream &os, CPoint<dim> &obj)
{
  obj.Dump(os);
  return os;
}

#endif
