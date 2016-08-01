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

// Enable this to force testing of array access
#undef ARRAY_DEBUG
#define FOR_RELEASE

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <Misha/cmdLineParser.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/LinearSolvers.h>
#include <Misha/Timer.h>
#include <Misha/FEM.h>

cmdLineParameter< char* > In( "in" ) , Out( "out" ) , Sphere( "sphere" );
cmdLineParameter< int > Iters( "iters" , 1 );
cmdLineParameter< float > ValueWeight( "vWeight" , 1e4f ) , GradientWeight( "gWeight" , 1.f ) , GradientScale( "gScale" , 1.f );
cmdLineReadable Verbose( "verbose" );

cmdLineReadable* params[] = { &In , &Out , &ValueWeight , &GradientWeight , &Verbose , &Iters , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output mesh>]\n" , Out.name );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s <value weight>=%f]\n" , ValueWeight.name , ValueWeight.value );
#ifndef FOR_RELEASE
	printf( "\t[--%s <gradient/smoothing weight>=%f]\n" , GradientWeight.name , GradientWeight.value );
#endif // !FOR_RELEASE
	printf( "\t[--%s]\n" , Verbose.name );
}

template< class Real >
void WriteMesh( char* fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& normals )
{
	std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].color = Point3D< float >( normals[i] + Point3D< Real >( 1 , 1 , 1 ) ) * 128.f;
	PlyWriteTriangles( fileName , _vertices ,  triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}
template< class Real >
void WriteColorMesh( char* fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors )
{
	std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].color = Point3D< float >( colors[i] );
	PlyWriteTriangles( fileName , _vertices ,  triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}
template< class Real >
void WriteMesh( char* fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices )
{
	std::vector< PlyVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] );
	PlyWriteTriangles( fileName , _vertices ,  triangles , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
}


template< class Real >
void _main( void )
{
#if defined( USE_CHOLMOD )
	typedef CholmodSolver Solver;
#elif defined( USE_EIGEN )
	typedef EigenSolverCholeskyLLt< Real , typename SparseMatrix< Real , int >::RowIterator > Solver;
#else // !USE_CHOLMOD && !USE_EIGEN
#error "Uknown solver type"
#endif // USE_CHOLMOD || USE_EIGEN
	int file_type;
	std::vector< TriangleIndex > triangles;
	std::vector< Point3D< Real > > vertices , normals;

	//////////////////////
	// Read in the data //
	{
		std::vector< PlyOrientedVertex< float > > _vertices;
		bool vertexFlags[ PlyOrientedVertex< float >::ReadComponents ];
		PlyReadTriangles( In.value , _vertices ,  triangles , PlyOrientedVertex< float >::ReadProperties , vertexFlags , PlyOrientedVertex< float >::ReadComponents , file_type );
		vertices.resize( _vertices.size() ) , normals.resize( _vertices.size() );
		for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].point );
		if( vertexFlags[3] && vertexFlags[4] && vertexFlags[5] ) for( int i=0 ; i<_vertices.size() ; i++ ) normals[i] = Point3D< Real >( _vertices[i].normal );
		else
		{
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< Real > v[] = { Point3D< Real >( _vertices[ triangles[i][0] ].point ) , Point3D< Real >( _vertices[ triangles[i][1] ].point ) , Point3D< Real >( _vertices[ triangles[i][2] ].point ) };
				Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
				for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
			}
			for( int i=0 ; i<normals.size() ; i++ ) normals[i] /= (Real)Length( normals[i] );
		}
		if( Verbose.set ) printf( "Source Vertices / Triangles: %d / %d\n" , (int)vertices.size() , (int)triangles.size() );
	}
	// Read in the data //
	//////////////////////

	/////////////////////
	// Smooth the normals
	{
		FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
		mesh.setMetricFromEmbedding( GetPointer( vertices ) );
		mesh.makeUnitArea();
		SparseMatrix< Real , int > M , _M = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() , _S = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
		M.resize( 2*vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ )
		{
			M.SetRowSize( 2*i+ii , 2*_M.rowSizes[i] );
			for( int j=0 ; j<_M.rowSizes[i] ; j++ ) for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].N = _M[i][j].N*2+jj;
		}
		std::vector< Point3D< Real > > tangents( vertices.size()*2 );
		std::vector< Real > b( vertices.size()*2 ) , o( vertices.size()*2 );

		Solver solver( M , true );

		Real mWeight = (Real)ValueWeight.value , sWeight = (Real)GradientWeight.value;

		for( int iter=0 ; iter<Iters.value ; iter++ )
		{
			Timer t;

			// Set the tangent directions
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ )
			{
				Point3D< Real > v( 1 , 0 , 0 );
				if( fabs( Point3D< Real >::Dot( v , normals[i] ) )>0.99 ) v = Point3D< Real >( 0 , 1 , 0 );
				tangents[2*i+0] = Point3D< Real >::CrossProduct( normals[i] , v               ) ; tangents[2*i+0] /= Length( tangents[2*i+0] );
				tangents[2*i+1] = Point3D< Real >::CrossProduct( normals[i] , tangents[2*i+0] ) ; tangents[2*i+1] /= Length( tangents[2*i+1] );
			}

			// Solve for the tangent offsets minimizing the dirichlet energy:
			// E( o1 , o2 ) = || \sum o[i] * T[i] ||^2 + e * || \nabla( \sum n[i] + o[i] * T[i] ) ||^2
			//              = o^t * T^t * M * T * o + e * [ o^t * T^t * S * T * o + 2 * o^t * T^t * S * n + n^t * S * n ]
			// \nabla E = 0:
			// 0 = T^t * ( M + e * S ) * T * o + e * T^t * S * n
			{
				Timer t;
#pragma omp parallel for 
				for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ ) 
				{
					b[2*i+ii] = 0;
					for( int j=0 ; j<_M.rowSizes[i] ; j++ )
					{
						for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].Value = ( _M[i][j].Value*mWeight + _S[i][j].Value * sWeight ) * Point3D< Real >::Dot( tangents[2*i+ii] , tangents[ 2*_M[i][j].N+jj ] );
						b[2*i+ii] -= _S[i][j].Value * Point3D< Real >::Dot( normals[ _S[i][j].N ] , tangents[2*i+ii] ) * sWeight;
					}
				}
				if( Verbose.set ) printf( "\tSet system matrix[%d]: %.2f(s)\n" , iter+1 , t.elapsed() );
			}
			{
				Timer t;
				solver.update( M );
				solver.solve( GetPointer( b ) , GetPointer( o ) );
#pragma omp parallel for
				for( int i=0 ; i<vertices.size() ; i++ ) normals[i] += tangents[2*i+0] * o[2*i+0] + tangents[2*i+1] * o[2*i+1] , normals[i] /= Length( normals[i] );
				if( Verbose.set ) printf( "\tSolved system[%d]: %.2f(s)\n" , iter+1 , t.elapsed() );
			}
		}
	}
	// Smooth the normals
	/////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		std::vector< PlyOrientedVertex< float > > _vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].normal = Point3D< float >( normals[i] );
		PlyWriteTriangles( Out.value , _vertices , triangles , PlyOrientedVertex< float >::WriteProperties , PlyOrientedVertex< float >::WriteComponents , file_type );
	}
	// Output the result
	////////////////////
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	_main< double >();
	return EXIT_SUCCESS;
}
