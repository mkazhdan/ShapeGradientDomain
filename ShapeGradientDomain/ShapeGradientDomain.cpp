/*
Copyright (c) 2016, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifdef USE_EIGEN
// Enable this if you have an MKL-backed Eigen implementation
#undef EIGEN_USE_MKL_ALL
#endif // USE_EIGEN
// Enable this to force testing of array access
#undef ARRAY_DEBUG
#define FOR_RELEASE


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <Misha/cmdLineParser.h>
#include <Misha/Timer.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/LinearSolvers.h>
#include <Misha/FEM.h>

cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameter< float > ValueWeight( "vWeight" , 1e4f ) , GradientWeight( "gWeight" , 1.f ) , GradientScale( "gScale" , 1.f );
cmdLineParameter< float > CurvatureWeight( "kWeight" , 0.f );
cmdLineReadable UseColors( "useColors" ) , Verbose( "verbose" ) , Anisotropic( "aniso" );

cmdLineReadable* params[] = { &In , &Out , &ValueWeight , &GradientWeight , &GradientScale , &UseColors , &CurvatureWeight , &Anisotropic , &Verbose , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input geometry>\n" , In.name );
	printf( "\t[--%s <output geometry>]\n" , Out.name );
	printf( "\t[--%s <value weight>=%f]\n" , ValueWeight.name , ValueWeight.value );
#ifndef FOR_RELEASE
	printf( "\t[--%s <gradient weight>=%f]\n" , GradientWeight.name , GradientWeight.value );
#endif // FOR_RELEASE
	printf( "\t[--%s <gradient scale>=%f]\n" , GradientScale.name , GradientScale.value );
	printf( "\t[--%s <curvature weight>=%f]\n" , CurvatureWeight.name , CurvatureWeight.value );
	printf( "\t[--%s]\n" , UseColors.name );
	printf( "\t[--%s]\n" , Anisotropic.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< class Real > Real Clamp( Real v , Real min , Real max ){ return std::min< Real >( max , std::max< Real >( min , v ) ); }
template< class Real > Point3D< Real > ClampColor( Point3D< Real > c ){ return Point3D< Real >( Clamp( c[0] , (Real)0 , (Real)255 ) , Clamp( c[1] , (Real)0 , (Real)255 ) , Clamp( c[2] , (Real)0 , (Real)255 ) ); }

// Compute the SVD factorization of a _SYMMETRIC_ matrix M as M = R^t * D * R
template< class Real > void SVD( const SquareMatrix< Real , 2 >& M , SquareMatrix< Real , 2 >& D , SquareMatrix< Real , 2 >& R )
{
	R = D = SquareMatrix< Real , 2 >::Identity();
	// The characteristic polynomial is:
	//		P(x) = x^2 - tr(M) * x + det(M)
	// The roots are:
	//		x = ( tr(M) +/- sqrt( tr(M) * tr(M) - 4 * det(M) ) / 2
	Real tr = M.trace() , det = M.determinant();
	Real disc = tr*tr - 4 * det;
	if( disc<0 ) fprintf( stderr , "[WARNING] Negative discriminant set to zero: %g\n" , disc ) , disc = 0;
	disc = (Real)sqrt( disc );
	Real e1 = ( tr - disc ) / 2 , e2 = ( tr + disc ) / 2;
	D(0,0) = e1 , D(1,1) = e2;

	// Assuming that the columns of R are (x,y)^t and (-y,x)^t, we get:
	//		 x*x * e1 + y*y * e2 = M(0,0)
	//		 x*x * e2 + y*y * e1 = M(1,1)
	//		-x*y * e1 + x*y * e2 = M(0,1)
	// Using the first two equations we get the values of x*x and y*y.
	// We can then choose the sign minimizing the third equation.

	SquareMatrix< Real , 2 > temp;
	temp(0,0) = temp(1,1) = e1 , temp(1,0) = temp(0,1) = e2;
	// The determinant vanishes if e1 = e2 [Which is good because then there is multiplicity.]
	if( !temp.determinant() ) return;

	Point2D< Real > xx_yy = temp.inverse() * Point2D< Real >( M(0,0) , M(1,1) );
	if( xx_yy[0]<0 || xx_yy[1]<0 ) fprintf( stderr , "[ERROR] The square is negative: %g , %g >=0\n" , xx_yy[0] , xx_yy[1] ) , exit( 0 );
	Real x = (Real)sqrt( xx_yy[0] ) , y = (Real)sqrt( xx_yy[1] );
	if( fabs( -x*y * e1 + x*y * e2 - M(0,1) )>fabs( x*y * e1 - x*y * e2 - M(0,1) ) ) x = -x;
	R(0,0) = x , R(0,1) = y , R(1,0) = -y , R(1,1) = x;
}
template< class Real > SquareMatrix< Real , 2 > SecondFundamentalForm( const Point3D< Real > v[3] , const Point3D< Real > n[3] )
{
	SquareMatrix< Real , 2 > II;
	Point3D< Real > dv[] = { v[1]-v[0] , v[2]-v[0] } , dn[] = { n[1]-n[0] , n[2]-n[0] };
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) II( i , j ) = Point3D< Real >::Dot( dv[i] , dn[j] );
	II(1,0) = II (0,1) = ( II(0,1) + II(1,0) ) / (Real)2.;
	return II;
}
template< class Real > SquareMatrix< Real , 2 > SecondFundamentalForm( Point3D< Real > v0 , Point3D< Real > v1 , Point3D< Real > v2 , Point3D< Real > n0 , Point3D< Real > n1 , Point3D< Real > n2 )
{
	Point3D< Real > v[] = { v0 , v1 , v2 } , n[] = { n0 , n1 , n2 };
	return SecondFundamentalForm( v , n );
}

template< class Real >
int _main( void )
{
#if defined( USE_CHOLMOD )
	typedef CholmodSolver Solver;
#elif defined( USE_EIGEN )
	typedef EigenSolverCholeskyLLt< Real , typename SparseMatrix< Real , int >::RowIterator > Solver;
#else // !USE_CHOLMOD && !USE_EIGEN
#error "Uknown solver type"
#endif // USE_CHOLMOD || USE_EIGEN

	std::vector< TriangleIndex > triangles;
	std::vector< PlyOrientedColorVertex< float > > vertices;

	int file_type;
	bool hasNormals , hasColors;

	//////////////////////
	// Read in the data //
	{
		bool readFlags[ PlyOrientedColorVertex< float >::ReadComponents ];
		PlyReadTriangles( In.value , vertices , triangles , PlyOrientedColorVertex< float >::ReadProperties , readFlags , PlyOrientedColorVertex< float >::ReadComponents , file_type );

		hasNormals = ( readFlags[3] && readFlags[4] && readFlags[5] );
		hasColors = ( readFlags[6] && readFlags[7] && readFlags[8] ) || ( readFlags[9] && readFlags[10] && readFlags[11] );

		if( !hasNormals )
		{
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].normal = Point3D< float >();
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< float > n = Point3D< float >::CrossProduct( vertices[ triangles[i][1] ].point - vertices[ triangles[i][0] ].point , vertices[ triangles[i][2] ].point - vertices[ triangles[i][0] ].point );
				for( int j=0 ; j<3 ; j++ ) vertices[ triangles[i][j] ].normal += n;
			}
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].normal /= (float)Length( vertices[i].normal );
		}
		if( !hasColors ) for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].color = ( Point3D< float >( 1.f , 1.f , 1.f ) + vertices[i].normal ) / 2.f;
		if( Verbose.set ) printf( "Source Vertices / Triangles: %d / %d\n" , (int)vertices.size() , (int)triangles.size() );
	}
	// Read in the data //
	//////////////////////
	
	std::vector< Point3D< Real > > signal( vertices.size() );

	//////////////////
	// Set the data //
	if( UseColors.set ) for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].color );
	else                for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].point );
	// Set the data //
	//////////////////

	////////////////////
	// Smooth the signal
	{
		// Read in the mesh and normalize the scale
		FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
		mesh.setMetricFromEmbedding( GetPointer( vertices ) );
		Real area = mesh.area();
		mesh.makeUnitArea();

		////////////////////
		// Modify the metric
		if( CurvatureWeight.value>0 )
		{
			Timer t;
			Real s = (Real)( 1. / sqrt(area) );
#pragma omp parallel for
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ].point  ) * s , Point3D< Real >( vertices[ triangles[i][1] ].point  ) * s , Point3D< Real >( vertices[ triangles[i][2] ].point  ) * s };
				Point3D< Real > n[] = { Point3D< Real >( vertices[ triangles[i][0] ].normal )     , Point3D< Real >( vertices[ triangles[i][1] ].normal )     , Point3D< Real >( vertices[ triangles[i][2] ].normal )     };
				SquareMatrix< Real , 2 > II = SecondFundamentalForm( v , n );

				// The columns of A_inverse given an orthonormal frame with respect to g.
				SquareMatrix< Real , 2 > A = FEM::TensorRoot( mesh.g(i) ) , A_inverse = A.inverse() , D , R;
				// The matrix A_inverse^t * II * A_inverse gives the second fundamental form with respect to a basis that is orthonormal w.r.t. g
				SVD( A_inverse.transpose() * II * A_inverse , D , R );

				if( Anisotropic.set ) D(0,0) *= D(0,0) , D(1,1) *= D(1,1);
				else                  D(0,0) = D(1,1) = ( D(0,0) * D(0,0) + D(1,1) * D(1,1) ) / (Real)2.;
				mesh.g( i ) = A * R.transpose() * ( SquareMatrix< Real , 2 >::Identity() + D * CurvatureWeight.value ) * R * A;
			}
			if( Verbose.set ) printf( "\tUpdated metric: %.2f(s)\n" , t.elapsed() );
		}
		// Modify the metric
		////////////////////

		///////////////////////////////
		// Compute and solve the system
		{
			Timer t;
			SparseMatrix< Real , int > mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() * ValueWeight.value;
			SparseMatrix< Real , int > stiffness = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >() * GradientWeight.value;

			Solver solver( mass + stiffness );
			Pointer( Real ) x = AllocPointer< Real >( vertices.size() );
			Pointer( Real ) b = AllocPointer< Real >( vertices.size() );
			for( int c=0 ; c<3 ; c++ )
			{
				for( int i=0 ; i<vertices.size() ; i++ ) x[i] = signal[i][c];
				mass.Multiply( x , b );
				for( int i=0 ; i<vertices.size() ; i++ ) x[i] *= GradientScale.value;
				stiffness.Multiply( x , b , MULTIPLY_ADD );
				solver.solve( b , x );
				for( int i=0 ; i<vertices.size() ; i++ ) signal[i][c] = x[i];
			}
			if( Verbose.set ) printf( "\tSolved the system: %.2f(s)\n" , t.elapsed() );
		}
		// Compute and solve the system
		///////////////////////////////
	}
	// Smooth the signal
	////////////////////

	////////////////////////////////////////////
	// Copy the processedsignal to the vertices
	if( UseColors.set ){ hasColors = true ; for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].color = ClampColor( Point3D< float >( signal[i] ) ); }
	else                                    for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].point = Point3D< float >( signal[i] );
	// Copy the processedsignal to the vertices
	////////////////////////////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		if( hasColors && hasNormals ) PlyWriteTriangles( Out.value , vertices , triangles , PlyOrientedColorVertex< float >::WriteProperties , PlyOrientedColorVertex< float >::WriteComponents , file_type );
		else if( hasColors )
		{
			std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = vertices[i].point , _vertices[i].color = vertices[i].color;
			PlyWriteTriangles( Out.value , _vertices , triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , file_type );
		}
		else if( hasNormals )
		{
			std::vector< PlyOrientedVertex< float > > _vertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = vertices[i].point , _vertices[i].normal = vertices[i].normal;
			PlyWriteTriangles( Out.value , _vertices , triangles , PlyOrientedVertex< float >::WriteProperties , PlyOrientedVertex< float >::WriteComponents , file_type );
		}
		else
		{
			std::vector< PlyVertex< float > > _vertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = vertices[i].point;
			PlyWriteTriangles( Out.value , _vertices , triangles , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , file_type );
		}
	}
	// Output the result
	////////////////////
	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	return _main< double >();
}
