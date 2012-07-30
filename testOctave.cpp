
// g++ -O2 -I/usr/include/octave-3.4.2 testOctave.cpp -L/usr/lib/octave/3.4.2 -loctave

#include<iostream>
#include<octave/config.h>
#include<octave/Matrix.h>

using namespace std;

double calcError( Matrix* homography, Matrix& x, Matrix& u )
{
  Matrix err1 = ((*homography)*x - u );
  Matrix err2 = (((*homography).inverse())*u - x );
  double distance1 = 0.0;
  double distance2 = 0.0;
  for( int i=0; i<x.cols(); ++i ) {
    distance1 += sqrt( err1(0,i)*err1(0,i)+err1(1,i)*err1(1,i) );
    distance2 += sqrt( err2(0,i)*err2(0,i)+err2(1,i)*err2(1,i) );
  }
  return distance1 + distance2;
}

void calcHomographyFromPoints( Matrix* homography, Matrix x, Matrix u, int normal_flag=true )
{
  Matrix T1(3,3, 0);
  Matrix T2(3,3, 0);
  int npoint = x.cols();

  if( normal_flag ) {
    // ¬•ª–ˆ‚Ì•½‹Ï‚ðŽZo 
    int dim = x.rows();
    ColumnVector mean(dim);
    mean = (x.sum(1)/npoint).column(0);
    double x_ave1 = mean(0);
    double x_ave2 = mean(1);
    mean = (u.sum(1)/npoint).column(0);
    double u_ave1 = mean(0);
    double u_ave2 = mean(1);

    // •½‹Ï‹——£‚ðsqrt(2)‚É•â³‚·‚éŒW”ŽZo 
    double x_distance = 0.0;
    double u_distance = 0.0;
    for( int i=0; i<npoint; ++i ) {
      x_distance += sqrt( (x(0,i)-x_ave1)*(x(0,i)-x_ave1)+(x(1,i)-x_ave2)*(x(1,i)-x_ave2) );
      u_distance += sqrt( (u(0,i)-u_ave1)*(u(0,i)-u_ave1)+(u(1,i)-u_ave2)*(u(1,i)-u_ave2) );
    }
    T1(0,0) = T1(1,1) = npoint * sqrt(2.0) / x_distance;
    T2(0,0) = T2(1,1) = npoint * sqrt(2.0) / u_distance;

    // ¬•ª–ˆ‚Ì•½‹Ï‚ð0‚É•â³‚·‚éŒW”ŽZo 
    T1(0,2) = -T1(0,0) * x_ave1;
    T1(1,2) = -T1(0,0) * x_ave2;
    T1(2,2) = 1;
    T2(0,2) = -T2(0,0) * u_ave1;
    T2(1,2) = -T2(0,0) * u_ave2;
    T2(2,2) = 1;

    // ã‹L‚Ì³‹K‰»s—ñ‚Å•â³
    x = T1*x;
    u = T2*u;
    //cout << "T1" << endl << T1 << endl;
    //cout << "T2" << endl << T2 << endl;
    //cout << "x" << endl << x << endl;
    //cout << "u" << endl << u << endl;
  }

  Matrix A(npoint*2, 9);
  for( int i=0; i<npoint; ++i ) {
    A(2*i,0) = x(0,i);
    A(2*i,1) = x(1,i);
    A(2*i,2) = 1.0;
    A(2*i,3) = 0;
    A(2*i,4) = 0;
    A(2*i,5) = 0;
    A(2*i,6) = -u(0,i)*x(0,i);
    A(2*i,7) = -u(0,i)*x(1,i);
    A(2*i,8) = -u(0,i);

    A(2*i+1,0) = 0;
    A(2*i+1,1) = 0;
    A(2*i+1,2) = 0;
    A(2*i+1,3) = x(0,i);
    A(2*i+1,4) = x(1,i);
    A(2*i+1,5) = 1.0;
    A(2*i+1,6) = -u(1,i)*x(0,i);
    A(2*i+1,7) = -u(1,i)*x(1,i);
    A(2*i+1,8) = -u(1,i);
  }
  //cout << "Original Matrix" << endl << A << endl;

  SVD svd(A);
  //	cout << "Left Singular Matrix" << endl << svd.left_singular_matrix() << endl;
  //	cout << "Singular Values" << endl << svd.singular_values() << endl;
  //	cout << "Right Singular Matrix" << endl << svd.right_singular_matrix() << endl;

  // Homography s—ñ‚ðŽZo 
  ColumnVector V = svd.right_singular_matrix().column(8);
  for( int j=0; j<3; j++ ) {
    for( int i=0; i<3; i++ ) {
      (*homography)(j,i) = V(3*j+i) / V(8);
    }
  }

  if( normal_flag ) {
    // ã‹L‚Ì³‹K‰»s—ñ‚Ì‰e‹¿‚ðœ‹Ž 
    *homography = T2.inverse() * (*homography) * T1;
    *homography /= (*homography)(2,2);
  }
}

void calcHomographyFromLines( Matrix* homography, Matrix l, Matrix m, int normal_flag = true )
{
  Matrix T1(3,3, 0);
  Matrix T2(3,3, 0);
  int nline = l.cols();
  
  if( normal_flag ) {
    // ¬•ª–ˆ‚Ì•½‹Ï‚ðŽZo 
    int dim = l.rows();
    ColumnVector mean(dim);
    mean = (l.sum(1)/nline).column(0);
    double l_ave1 = mean(0);
    double l_ave2 = mean(1);
    mean = (m.sum(1)/nline).column(0);
    double m_ave1 = mean(0);
    double m_ave2 = mean(1);
    
    // •½‹Ï‹——£‚ðsqrt(2)‚É•â³‚·‚éŒW”ŽZo 
    double l_distance = 0.0;
    double m_distance = 0.0;
    for( int i=0; i<nline; ++i ) {
      l_distance += sqrt( (l(0,i)-l_ave1)*(l(0,i)-l_ave1)+(l(1,i)-l_ave2)*(l(1,i)-l_ave2) );
      m_distance += sqrt( (m(0,i)-m_ave1)*(m(0,i)-m_ave1)+(m(1,i)-m_ave2)*(m(1,i)-m_ave2) );
    }
    T1(0,0) = T1(1,1) = nline * sqrt(2.0) / l_distance;
    T2(0,0) = T2(1,1) = nline * sqrt(2.0) / m_distance;

    // ¬•ª–ˆ‚Ì•½‹Ï‚ð0‚É•â³‚·‚éŒW”ŽZo 
    T1(0,2) = -T1(0,0) * l_ave1;
    T1(1,2) = -T1(0,0) * l_ave2;
    T1(2,2) = 1;
    T2(0,2) = -T2(0,0) * m_ave1;
    T2(1,2) = -T2(0,0) * m_ave2;
    T2(2,2) = 1;

    // ã‹L‚Ì³‹K‰»s—ñ‚Å•â³ 
    for( int i=0; i<nline; ++i ) {
      l(0,i) = l(0,i)/(T1(0,0) - T1(0,2)*l(0,i) - T1(1,2)*l(1,i));
      l(1,i) = l(1,i)/(T1(0,0) - T1(0,2)*l(0,i) - T1(1,2)*l(1,i));
      m(0,i) = m(0,i)/(T2(0,0) - T2(0,2)*m(0,i) - T2(1,2)*m(1,i));
      m(1,i) = m(1,i)/(T2(0,0) - T2(0,2)*m(0,i) - T2(1,2)*m(1,i));
    }
    //cout << "T1" << endl << T1 << endl;
    //cout << "T2" << endl << T2 << endl;
    //cout << "l" << endl << l << endl;
    //cout << "m" << endl << m << endl;
  }

  Matrix A(nline*2, 9);
  for( int i=0; i<nline; ++i ) {
    A(2*i,0) = -m(0,i);
    A(2*i,1) = 0;
    A(2*i,2) = l(0,i)*m(0,i);
    A(2*i,3) = -m(1,i);
    A(2*i,4) = 0;
    A(2*i,5) = -l(0,i)*m(1,i);
    A(2*i,6) = -1.0;
    A(2*i,7) = 0;
    A(2*i,8) = l(0,i);

    A(2*i+1,0) = 0;
    A(2*i+1,1) = -m(0,i);
    A(2*i+1,2) = l(1,i)*m(0,i);
    A(2*i+1,3) = 0;
    A(2*i+1,4) = -m(1,i);
    A(2*i+1,5) = l(1,i)*m(1,i);
    A(2*i+1,6) = 0;
    A(2*i+1,7) = -1.0;
    A(2*i+1,8) = l(1,i);
  }
  //cout << "Original Matrix" << endl << A << endl;

  SVD svd(A);

  //	cout << "Left Singular Matrix" << endl << svd.left_singular_matrix() << endl;
  //	cout << "Singular Values" << endl << svd.singular_values() << endl;
  //	cout << "Right Singular Matrix" << endl << svd.right_singular_matrix() << endl;

  ColumnVector V = svd.right_singular_matrix().column(8);
  for( int j=0; j<3; j++ ) {
    for( int i=0; i<3; i++ ) {
      (*homography)(j,i) = V(3*j+i) / V(8);
    }
  }

  if( normal_flag ) {
    *homography = T2.inverse() * (*homography) * T1;
    *homography /= (*homography)(2,2);
  }
}

void calcHomographyFromPointsLines( Matrix* homography, Matrix x, Matrix u, Matrix l, Matrix m, int normal_flag=true )
{
  Matrix T1(3,3, 0);
  Matrix T2(3,3, 0);
  int npoint = x.cols();
  int nline = l.cols();

  if( normal_flag ) {
    // ¬•ª–ˆ‚Ì•½‹Ï‚ðŽZo
    int dim = x.rows();
    ColumnVector mean(dim);
    mean = (x.sum(1)/npoint).column(0);
    double x_ave1 = mean(0);
    double x_ave2 = mean(1);
    mean = (u.sum(1)/npoint).column(0);
    double u_ave1 = mean(0);
    double u_ave2 = mean(1);

    // •½‹Ï‹——£‚ðsqrt(2)‚É•â³‚·‚éŒW”ŽZo 
    double x_distance = 0.0;
    double u_distance = 0.0;
    for( int i=0; i<npoint; ++i ) {
      x_distance += sqrt( (x(0,i)-x_ave1)*(x(0,i)-x_ave1)+(x(1,i)-x_ave2)*(x(1,i)-x_ave2) );
      u_distance += sqrt( (u(0,i)-u_ave1)*(u(0,i)-u_ave1)+(u(1,i)-u_ave2)*(u(1,i)-u_ave2) );
    }
    T1(0,0) = T1(1,1) = npoint * sqrt(2.0) / x_distance;
    T2(0,0) = T2(1,1) = npoint * sqrt(2.0) / u_distance;

    // ¬•ª–ˆ‚Ì•½‹Ï‚ð0‚É•â³‚·‚éŒW”ŽZo 
    T1(0,2) = -T1(0,0) * x_ave1;
    T1(1,2) = -T1(0,0) * x_ave2;
    T1(2,2) = 1;
    T2(0,2) = -T2(0,0) * u_ave1;
    T2(1,2) = -T2(0,0) * u_ave2;
    T2(2,2) = 1;

    // ã‹L‚Ì³‹K‰»s—ñ‚Å•â³
    x = T1*x;
    u = T2*u;
    for( int i=0; i<nline; ++i ) {
      l(0,i) = l(0,i)/(T1(0,0) - T1(0,2)*l(0,i) - T1(1,2)*l(1,i));
      l(1,i) = l(1,i)/(T1(0,0) - T1(0,2)*l(0,i) - T1(1,2)*l(1,i));
      m(0,i) = m(0,i)/(T2(0,0) - T2(0,2)*m(0,i) - T2(1,2)*m(1,i));
      m(1,i) = m(1,i)/(T2(0,0) - T2(0,2)*m(0,i) - T2(1,2)*m(1,i));
    }
    //cout << "T1" << endl << T1 << endl;
    //cout << "T2" << endl << T2 << endl;
    //cout << "x" << endl << x << endl;
    //cout << "u" << endl << u << endl;
    //cout << "l" << endl << l << endl;
    //cout << "m" << endl << m << endl;
  }

  Matrix A((npoint+nline)*2, 9);
  for( int i=0; i<npoint; ++i ) {
    A(2*i,0) = x(0,i);
    A(2*i,1) = x(1,i);
    A(2*i,2) = 1.0;
    A(2*i,3) = 0;
    A(2*i,4) = 0;
    A(2*i,5) = 0;
    A(2*i,6) = -u(0,i)*x(0,i);
    A(2*i,7) = -u(0,i)*x(1,i);
    A(2*i,8) = -u(0,i);

    A(2*i+1,0) = 0;
    A(2*i+1,1) = 0;
    A(2*i+1,2) = 0;
    A(2*i+1,3) = x(0,i);
    A(2*i+1,4) = x(1,i);
    A(2*i+1,5) = 1.0;
    A(2*i+1,6) = -u(1,i)*x(0,i);
    A(2*i+1,7) = -u(1,i)*x(1,i);
    A(2*i+1,8) = -u(1,i);
  }
  for( int i=0; i<nline; ++i ) {
    A(2*i+2*npoint,0) = -m(0,i);
    A(2*i+2*npoint,1) = 0;
    A(2*i+2*npoint,2) = l(0,i)*m(0,i);
    A(2*i+2*npoint,3) = -m(1,i);
    A(2*i+2*npoint,4) = 0;
    A(2*i+2*npoint,5) = -l(0,i)*m(1,i);
    A(2*i+2*npoint,6) = -1.0;
    A(2*i+2*npoint,7) = 0;
    A(2*i+2*npoint,8) = l(0,i);

    A(2*i+2*npoint+1,0) = 0;
    A(2*i+2*npoint+1,1) = -m(0,i);
    A(2*i+2*npoint+1,2) = l(1,i)*m(0,i);
    A(2*i+2*npoint+1,3) = 0;
    A(2*i+2*npoint+1,4) = -m(1,i);
    A(2*i+2*npoint+1,5) = l(1,i)*m(1,i);
    A(2*i+2*npoint+1,6) = 0;
    A(2*i+2*npoint+1,7) = -1.0;
    A(2*i+2*npoint+1,8) = l(1,i);
  }
  //cout << "Original Matrix" << endl << A << endl;

  SVD svd(A);
  //	cout << "Left Singular Matrix" << endl << svd.left_singular_matrix() << endl;
  //	cout << "Singular Values" << endl << svd.singular_values() << endl;
  //	cout << "Right Singular Matrix" << endl << svd.right_singular_matrix() << endl;

  // Homography s—ñ‚ðŽZo 
  ColumnVector V = svd.right_singular_matrix().column(8);
  for( int j=0; j<3; j++ ) {
    for( int i=0; i<3; i++ ) {
      (*homography)(j,i) = V(3*j+i) / V(8);
    }
  }

  if( normal_flag ) {
    // ã‹L‚Ì³‹K‰»s—ñ‚Ì‰e‹¿‚ðœ‹Ž 
    *homography = T2.inverse() * (*homography) * T1;
    *homography /= (*homography)(2,2);
  }
}


int main( int argc, char* argv[] )
{
  srand(time(NULL));

  Matrix x(3,9, 1.0);
  //x(0,0)=1.0; x(0,1)=-1.0; x(0,2)=-1.0; x(0,3)= 1.0;
  //x(1,0)=1.0; x(1,1)= 1.0; x(1,2)=-1.0; x(1,3)=-1.0;
  x(0,0)=100; x(0,1)=110; x(0,2)=120; x(0,3)=100; x(0,4)=110; x(0,5)=120; x(0,6)=100; x(0,7)=110; x(0,8)=120;
  x(1,0)=100; x(1,1)=100; x(1,2)=100; x(1,3)=110; x(1,4)=110; x(1,5)=110; x(1,6)=120; x(1,7)=120; x(1,8)=120;

  Matrix u(3,9, 1.0);
  //u(0,0)=5.0+(double)rand()/RAND_MAX; u(0,1)= 1.0+(double)rand()/RAND_MAX; u(0,2)= 1.0+(double)rand()/RAND_MAX; u(0,3)= 5.0+(double)rand()/RAND_MAX;
  //u(1,0)=5.0+(double)rand()/RAND_MAX; u(1,1)= 5.0+(double)rand()/RAND_MAX; u(1,2)= 1.0+(double)rand()/RAND_MAX; u(1,3)= 1.0+(double)rand()/RAND_MAX;
  u(0,0)=186; u(0,1)=228; u(0,2)=281; u(0,3)=183; u(0,4)=228; u(0,5)=273; u(0,6)=180; u(0,7)=227; u(0,8)=275;
  u(1,0)= 81; u(1,1)= 81; u(1,2)= 82; u(1,3)=104; u(1,4)=104; u(1,5)=103; u(1,6)=129; u(1,7)=128; u(1,8)=128;

  Matrix l(3,4, 1.0);
  Matrix m(3,4, 1.0);

  for( int i=0; i<4; i++ ){
    int j = i>=3? 0: i+1;
    double d = x(0,j)*x(1,i)-x(0,i)*x(1,j);
    l(0,i) = d==0? 0:(x(1,j)-x(1,i))/d;
    l(1,i) = d==0? 0:(x(0,i)-x(0,j))/d;

    d = u(0,j)*u(1,i)-u(0,i)*u(1,j);
    m(0,i) = d==0? 0:(u(1,j)-u(1,i))/d;
    m(1,i) = d==0? 0:(u(0,i)-u(0,j))/d;
  }
  //cout << "original point" << endl << x << endl << "detected point" << endl << u << endl;

  Matrix* homography = new Matrix(3,3);
  
  calcHomographyFromPoints( homography, x, u, true );
  //cout << "Homography Matrix" << endl << *homography << endl;
  cout << calcError( homography, x, u );

  calcHomographyFromLines( homography, l, m, true );
  //cout << "Homography Matrix" << endl << *homography << endl;
  cout << ", " <<  calcError( homography, x, u );

  calcHomographyFromPointsLines( homography, x, u, l, m, true );
  //cout << "Homography Matrix" << endl << *homography << endl;
  cout << ", " << calcError( homography, x, u ) << endl;

  delete homography;
  return 0;
}
