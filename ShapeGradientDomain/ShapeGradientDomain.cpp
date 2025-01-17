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

#include <Include/PreProcessor.h>


#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <mutex>
#include <Misha/CmdLineParser.h>
#include <Misha/Miscellany.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#include <Eigen/Sparse>
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // EIGEN_USE_MKL_ALL
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< float > ValueWeight( "vWeight" , 1.f ) , GradientWeight( "gWeight" , 1e-4f ) , GradientScale( "gScale" , 1.f );
Misha::CmdLineParameter< float > CurvatureWeight( "kWeight" , 0.f );
Misha::CmdLineReadable Anisotropic( "anisotropic" ) , UseColors( "useColors" ) , Verbose( "verbose" );
Misha::CmdLineReadable* params[] = { &In , &Out , &ValueWeight , &GradientWeight , &GradientScale , &Anisotropic , &UseColors , &CurvatureWeight , &Verbose , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input geometry>\n" , In.name.c_str() );
	printf( "\t[--%s <output geometry>]\n" , Out.name.c_str() );
	printf( "\t[--%s <value weight>=%f]\n" , ValueWeight.name.c_str() , ValueWeight.value );
	printf( "\t[--%s <gradient weight>=%f]\n" , GradientWeight.name.c_str() , GradientWeight.value );
	printf( "\t[--%s <gradient scale>=%f]\n" , GradientScale.name.c_str() , GradientScale.value );
	printf( "\t[--%s <curvature weight>=%f]\n" , CurvatureWeight.name.c_str() , CurvatureWeight.value );
	printf( "\t[--%s]\n" , Anisotropic.name.c_str() );
	printf( "\t[--%s]\n" , UseColors.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL

template< class Real > Real Clamp( Real v , Real min , Real max ){ return std::min< Real >( max , std::max< Real >( min , v ) ); }
template< class Real > Point3D< Real > ClampColor( Point3D< Real > c ){ return Point3D< Real >( Clamp( c[0] , (Real)0 , (Real)255 ) , Clamp( c[1] , (Real)0 , (Real)255 ) , Clamp( c[2] , (Real)0 , (Real)255 ) ); }

template< typename Real >
Eigen::SparseMatrix< Real > ToEigen( const SparseMatrix< Real , int > & M )
{
	Eigen::SparseMatrix< Real > _M;
	_M.resize( M.Rows() , M.Rows() );
	std::vector< Eigen::Triplet< Real > > triplets( M.Entries() );
	for( unsigned int i=0 , idx=0 ; i<M.Rows() ; i++ ) for( unsigned int j=0 ; j<M.rowSizes[i] ; j++ , idx++ )
		triplets[idx] = Eigen::Triplet< Real >( i , M[i][j].N , M[i][j].Value );
	_M.setFromTriplets( triplets.begin() , triplets.end() );
	return _M;
}

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
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
	using Vertex = typename Factory::VertexType;

	std::vector< TriangleIndex > triangles;
	std::vector< Vertex > vertices;

	int file_type;
	bool hasNormals , hasColors;

	//////////////////////
	// Read in the data //
	{
		Factory factory;
		bool *readFlags = new bool[ factory.plyReadNum() ];
		PLY::ReadTriangles( In.value , factory , vertices , triangles , readFlags , file_type );

		hasNormals = factory.template plyValidReadProperties<1>( readFlags );
		hasColors  = factory.template plyValidReadProperties<2>( readFlags );

		if( !hasNormals )
		{
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() = Point3D< float >();
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< float > n = Point3D< float >::CrossProduct( vertices[ triangles[i][1] ].template get<0>() - vertices[ triangles[i][0] ].template get<0>() , vertices[ triangles[i][2] ].template get<0>() - vertices[ triangles[i][0] ].template get<0>() );
				for( int j=0 ; j<3 ; j++ ) vertices[ triangles[i][j] ].template get<1>() += n;
			}
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() /= Point< float , 3 >::Length( vertices[i].template get<1>() );
		}
		if( !hasColors ) for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<2>() = ( Point3D< float >( 1.f , 1.f , 1.f ) + vertices[i].template get<1>() ) / 2.f;
		delete[] readFlags;
		if( Verbose.set ) std::cout << "Source Vertices / Triangles: " << vertices.size() << " / " << triangles.size() << std::endl;
	}
	// Read in the data //
	//////////////////////

	std::vector< Point3D< Real > > signal( vertices.size() );

	//////////////////
	// Set the data //
	if( UseColors.set ) for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get<2>() );
	else                for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get<0>() );
	// Set the data //
	//////////////////

	////////////////////
	// Smooth the signal
	{
		// Read in the mesh and normalize the scale
		FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
		mesh.template setMetricFromEmbedding< 3 >( [&]( unsigned int i ){ return vertices[i].template get<0>(); } );
		Real area = mesh.area();
		mesh.makeUnitArea();

		////////////////////
		// Modify the metric
		if( CurvatureWeight.value>0 )
		{
			Miscellany::Timer timer;
			Real s = (Real)( 1. / sqrt(area) );
			std::vector< Real > curvatureValues;


			std::mutex mut;
			ThreadPool::ParallelFor
			(
				0 , triangles.size() ,
				[&]( unsigned int , size_t i )
				{
					Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ].template get<0>() ) * s , Point3D< Real >( vertices[ triangles[i][1] ].template get<0>() ) * s , Point3D< Real >( vertices[ triangles[i][2] ].template get<0>() ) * s };
					Point3D< Real > n[] = { Point3D< Real >( vertices[ triangles[i][0] ].template get<1>() )     , Point3D< Real >( vertices[ triangles[i][1] ].template get<1>() )     , Point3D< Real >( vertices[ triangles[i][2] ].template get<1>() )     };
					SquareMatrix< Real , 2 > II = SecondFundamentalForm( v , n );

					// The columns of A_inverse given an orthonormal frame with respect to g.
					SquareMatrix< Real , 2 > A = FEM::TensorRoot( mesh.g(i) ) , A_inverse = A.inverse() , D , R;
					// The matrix A_inverse^t * II * A_inverse gives the second fundamental form with respect to a basis that is orthonormal w.r.t. g
					// The rows of R are the eigenvecctors.
					SVD( A_inverse.transpose() * II * A_inverse , D , R );
					D(0,0) *= D(0,0);
					D(1,1) *= D(1,1);
					if( !Anisotropic.set ) D(0,0) = D(1,1) = ( D(0,0) + D(1,1) ) / (Real)2.;
					D = SquareMatrix< Real , 2 >::Identity() + D * CurvatureWeight.value;
					// Multiplying by A gives the coefficients w.r.t. to an orthonormal basis.
					// Multiplying by R gives the coefficients w.r.t. to the principal curvature directions.
					mesh.g( i ) = A.transpose() * R.transpose() * D * R * A;
				}
			);

			if( Verbose.set ) std::cout << "\tUpdated metric: " << timer() << std::endl;
		}
		// Modify the metric
		////////////////////

		///////////////////////////////
		// Compute and solve the system
		{
			Miscellany::Timer timer;

			Eigen::SparseMatrix< Real > mass = ToEigen( mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() * ValueWeight.value );
			Eigen::SparseMatrix< Real > stiffness = ToEigen( mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >() * GradientWeight.value );

			Eigen::SparseMatrix< Real > M = mass + stiffness;
			double aTime = timer.elapsed();
			Solver solver;
			solver.analyzePattern( M );
			aTime = timer.elapsed() - aTime;
			double fTime = timer.elapsed();
			solver.factorize( M );
			fTime = timer.elapsed() - fTime;
			Eigen::VectorXd x;
			x.resize( vertices.size() );
			double sTime = 0;
			for( int c=0 ; c<3 ; c++ )
			{
				for( int i=0 ; i<vertices.size() ; i++ ) x[i] = signal[i][c];
				double t = timer.elapsed();
				Eigen::VectorXd b = mass * x + stiffness * x * GradientScale.value;
				x = solver.solve( b );
				sTime += timer.elapsed()-t;
				for( int i=0 ; i<vertices.size() ; i++ ) signal[i][c] = x[i];
			}
			if( Verbose.set ) std::cout << "\tSolved the system: " << timer() << ": " << aTime << " + " << fTime << " + " << sTime << std::endl;
		}
		// Compute and solve the system
		///////////////////////////////
	}
	// Smooth the signal
	////////////////////

	////////////////////////////////////////////
	// Copy the processed signal to the vertices
	if( UseColors.set ){ hasColors = true ; for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<2>() = ClampColor( Point3D< float >( signal[i] ) ); }
	else                                    for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<0>() = Point3D< float >( signal[i] );
	// Copy the processedsignal to the vertices
	////////////////////////////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		if( hasColors && hasNormals ) PLY::WriteTriangles( Out.value , Factory() , vertices , triangles , file_type );
		else if( hasColors )
		{
			using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i].template get<0>() , _vertices[i].template get<1>() = vertices[i].template get<2>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
		else if( hasNormals )
		{
			using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i].template get<0>() , _vertices[i].template get<1>() = vertices[i].template get<1>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
		else
		{
			using Factory = VertexFactory::PositionFactory< float , 3 >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = vertices[i].template get<0>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
	}
	// Output the result
	////////////////////
	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	return _main< double >();
}
