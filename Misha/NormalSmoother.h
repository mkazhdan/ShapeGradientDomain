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

#ifndef NORMAL_SMOOTHER
#define NORMAL_SMOOTHER

#include <Eigen/Sparse>
#include "SimplexMesh.h"
#include "Miscellany.h"
#include "MultiThreading.h"

namespace NormalSmoother
{
	template< unsigned int K , typename Index , typename EigenSolver=Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > , unsigned int Degree=1 >
	void Smooth( const std::vector< Point< double , K+1 > > &vertices , std::vector< Point< double , K+1 > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , unsigned int iters , double timeStep , bool unitMeasure=true );

	template< unsigned int K , typename EigenSolver=Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > >
	void Smooth( const Eigen::SparseMatrix< double > &mass , const Eigen::SparseMatrix< double > &stiffness , std::vector< Point< double , K+1 > > &normals , unsigned int iters , double timeStep );

	template< unsigned int K >
	Eigen::SparseMatrix< double > Prolongation( const std::vector< Point< double , K+1 > > &normals );

	template< unsigned int K >
	Eigen::SparseMatrix< double > DirectSum( const Eigen::SparseMatrix< double > &M );

	////////////////////
	// Implementation //
	////////////////////

	template< unsigned int K >
	struct _TangentFrame
	{
		static const unsigned int Dim = K+1;

		Point< double , Dim > t[K];

		const Point< double , Dim > &operator[]( unsigned int idx ) const { return t[idx]; }

		void set( Point< double , Dim > n )
		{
			Point< double , Dim > frame[Dim];
			{
				frame[0] = n;

				for( unsigned int d=1 ; d<Dim-1 ; d++ )
				{
					frame[d] = RandomSpherePoint< double , Dim >();
					while( true )
					{
						for( unsigned int dd=0 ; dd<d ; dd++ ) frame[d] -= Point< double , Dim >::Dot( frame[dd] , frame[d] ) * frame[dd];
						if( frame[d].squareNorm()>1e-10 )
						{
							frame[d] /= sqrt( frame[d].squareNorm() );
							break;
						}
					}
				}
				frame[Dim-1] = Point< double , Dim >::CrossProduct( frame );
			}
			for( unsigned int k=0 ; k<K ; k++ ) t[k] = frame[k+1];
		}
	};

	template< unsigned int K >
	Eigen::SparseMatrix< double > Prolongation( const std::vector< _TangentFrame< K > > &frames )
	{
		static const unsigned int Dim = K+1;

		Eigen::SparseMatrix< double > _P;
		_P.resize( frames.size()*Dim , frames.size()*K );
		_P.reserve( Eigen::VectorXi::Constant( frames.size()*K , Dim ) );
		ThreadPool::ParallelFor
		(
			0 , frames.size() ,
			[&]( unsigned int , size_t f ){ for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _P.insert( f*Dim+d , f*K+k ) = frames[f][k][d]; }
		);
		return _P;
	}

	template< unsigned int K >
	Eigen::SparseMatrix< double > Prolongation( const std::vector< Point< double , K+1 > > &normals )
	{
		static const unsigned int Dim = K+1;

		std::vector< _TangentFrame< K > > frames( normals.size() );
		ThreadPool::ParallelFor( 0 , normals.size() , [&]( unsigned int , size_t i ){ frames[i].set( normals[i] ); } );

		Eigen::SparseMatrix< double > _P;
		_P.resize( frames.size()*Dim , frames.size()*K );
		_P.reserve( Eigen::VectorXi::Constant( frames.size()*K , Dim ) );
		ThreadPool::ParallelFor
		(
			0 , frames.size() ,
			[&]( unsigned int , size_t f ){ for( unsigned int k=0 ; k<K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _P.insert( f*Dim+d , f*K+k ) = frames[f][k][d]; }
		);
		return _P;
	}

	template< unsigned int K >
	Eigen::SparseMatrix< double > DirectSum( const Eigen::SparseMatrix< double > &M )
	{
		static const unsigned int Dim = K+1;

		Eigen::SparseMatrix< double > _M;
		_M.resize( M.rows() * Dim , M.cols() * Dim );
		Eigen::VectorXi sizes( M.rows() * Dim );

		ThreadPool::ParallelFor
		(
			0 , M.outerSize() ,
			[&]( unsigned int , size_t k )
			{
				unsigned int rowSize = 0;
				for( Eigen::SparseMatrix<double>::InnerIterator it(M,k) ; it ; ++it ) rowSize++;
				for( unsigned int d=0 ; d<Dim ; d++ ) sizes[k*Dim+d] = rowSize;
			}
		);
		_M.reserve( sizes );

		ThreadPool::ParallelFor
		(
			0 , M.outerSize() ,
			[&]( unsigned int , size_t k )
			{
				for( Eigen::SparseMatrix<double>::InnerIterator it(M,k) ; it ; ++it )
					for( unsigned int d=0 ; d<Dim ; d++ ) _M.insert( it.row()*Dim+d , it.col()*Dim+d ) = it.value();
			}
		);

		return _M;
	}


	template< unsigned int K , typename Index , typename EigenSolver , unsigned int Degree >
	void Smooth( const std::vector< Point< double , K+1 > > &vertices , std::vector< Point< double , K+1 > > &normals , const std::vector< SimplexIndex< K , Index > > &simplices , unsigned int iters , double timeStep , bool unitMeasure )
	{
		static const unsigned int Dim = K+1;

		SimplexMesh< K , Degree > sMesh = SimplexMesh< K , Degree >::template Init< Dim , Index >( simplices , [&]( unsigned int v ){ return vertices[v]; } );
		if( unitMeasure ) sMesh.makeUnitVolume();
		Smooth< K >( sMesh.mass() , sMesh.stiffness() , normals , iters , timeStep );
	}

	template< unsigned int K , typename EigenSolver >
	void Smooth( const Eigen::SparseMatrix< double > &mass , const Eigen::SparseMatrix< double > &stiffness , std::vector< Point< double , K+1 > > &normals , unsigned int iters , double timeStep )
	{
		static const unsigned int Dim = K+1;

		Eigen::SparseMatrix< double > A , M = DirectSum< K >( mass ) , S = DirectSum< K >( stiffness );
		A = M + S * timeStep;

		EigenSolver solver;
		for( unsigned int iter=0 ; iter<iters ; iter++ )
		{
			Eigen::SparseMatrix< double > P = Prolongation< K >( normals );
			Eigen::SparseMatrix< double > Pt = P.transpose();


			// E( o ) = M( n + P*o - n ) + t * S( n + P*o )
			//        = o^t * P^t * ( M + t*S ) * P * o + 2 * t * o^t * P^t * S * n + ...
			// =>
			//  P^t * ( M + t*S ) * P * o = - t * P^t * S * n
			if( !iter ) solver.compute( Pt * A * P );
			else        solver.factorize( Pt * A * P );


			Eigen::VectorXd n( normals.size()*Dim );
			for( unsigned int i=0 ; i<normals.size() ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) n[i*Dim+d] = normals[i][d];
			Eigen::VectorXd b = - Pt * S * n * timeStep;
			n += P * solver.solve( b );

			ThreadPool::ParallelFor
			(
				0 , normals.size() ,
				[&]( unsigned int , size_t i )
				{
					for( unsigned int d=0 ; d<Dim ; d++ ) normals[i][d] = n[i*Dim+d];
					normals[i] /= Point< double , Dim >::Length( normals[i] );
				}
			);
		}
	}
}
#endif // NORMAL_SMOOTHER