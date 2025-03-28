/*
Copyright (c) 2025, Michael Kazhdan
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

#ifndef CURVATURE_METRIC_INCLUDED
#define CURVATURE_METRIC_INCLUDED

#include <Eigen/Dense>
#include <Misha/MultiThreading.h>
#include <Misha/Geometry.h>
#include <Misha/FEM.h>

namespace MishaK
{
	namespace CurvatureMetric
	{
		template< typename Real , unsigned int K >
		SquareMatrix< Real , K > SecondFundamentalForm( Simplex< Real , K+1 , K > s , const Point< Real , K+1 > n[K+1] , bool symmetric=true );

		template< typename Real , unsigned int K >
		Matrix< Real , K , K+1 > PrincipalCurvatures( Simplex< Real , K+1 , K > s , const Point< Real , K+1 > n[K+1] , Real curvatureValues[K] );

		template< typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 > ( unsigned int ) */ , typename NormalFunctor /* = std::function< Point< Real , 3 > ( unsigned int ) */ , typename CurvatureFunctor /* = std::function< Point< Real , 2 > ( Point< Real , 2 > ) > */ >
		void SetCurvatureMetric( FEM::RiemannianMesh< Real > &mesh , VertexFunctor && V , NormalFunctor && N , CurvatureFunctor && F );

		////////////////////
		// Implementation //
		////////////////////

		template< typename Real , unsigned int K >
		SquareMatrix< Real , K > SecondFundamentalForm( Simplex< Real , K+1 , K > s , const Point< Real , K+1 > n[K+1] , bool symmetric )
		{
			SquareMatrix< Real , K > II;
			Point< Real , K+1 > dv[K] , dn[K];
			for( unsigned int k=0 ; k<K ; k++ ) dv[k] = s[k+1]-s[0] , dn[k] = n[k+1]-n[0];
			for( int i=0 ; i<K ; i++ ) for( int j=0 ; j<K ; j++ ) II( i , j ) = Point< Real , K+1 >::Dot( dv[i] , dn[j] );
			if( symmetric ) return ( II + II.transpose() ) / (Real)2;
			else            return II;
		}

		template< typename Real , unsigned int K >
		Matrix< Real , K , K+1 > PrincipalCurvatures( Simplex< Real , K+1 , K > s , const Point< Real , K+1 > n[K+1] , Real curvatureValues[K] )
		{
			auto ToEigen = []( const SquareMatrix< Real , K > & M )
				{
					Eigen::Matrix< Real , K , K > _M;
					for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) _M(i,j) = M(j,i);
					return _M;
				};

			auto FromEigen = []( const Eigen::Matrix< Real , K , K > & M )
				{
					SquareMatrix< Real , K > _M;
					for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) _M(i,j) = M(j,i);
					return _M;
				};

			Point< Real , K+1 > d[K];
			for( unsigned int k=0 ; k<K ; k++ ) d[k] = s[k+1]-s[0];

			Matrix< Real , K , K+1 > D;
			for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<=K ; j++ ) D(i,j) = d[i][j];

			SquareMatrix< Real , K > II = SecondFundamentalForm< Real , K >( s , n , true );

			Eigen::GeneralizedSelfAdjointEigenSolver< Eigen::MatrixXd > ges( ToEigen( II ) , ToEigen( D.transpose() * D ) );
			for( unsigned int k=0 ; k<K ; k++ ) curvatureValues[k] = ges.eigenvalues()[k];

			return D * FromEigen( ges.eigenvectors() );
		}

		template< typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 > ( unsigned int ) */ , typename NormalFunctor /* = std::function< Point< Real , 3 > ( unsigned int ) */ , typename PrincipalCurvatureFunctor /* = std::function< Point< Real , 2 > ( unsigned int , Point< Real , 2 > ) > */ >
		void SetCurvatureMetric( FEM::RiemannianMesh< Real > &mesh , VertexFunctor && V , NormalFunctor && N , PrincipalCurvatureFunctor && PCF )
		{
			static const unsigned int K = 2;
			static const unsigned int Dim = K+1;

			static_assert( std::is_convertible_v<             VertexFunctor , std::function< Point< Real , Dim > ( unsigned int )                   > > , "[ERROR] VertexFunctor is poorly formed" );
			static_assert( std::is_convertible_v<             NormalFunctor , std::function< Point< Real , Dim > ( unsigned int )                   > > , "[ERROR] NormalFunctor is poorly formed" );
			static_assert( std::is_convertible_v< PrincipalCurvatureFunctor , std::function< Point< Real , K > ( unsigned int , Point< Real , K > ) > > , "[ERROR] PrincipalCurvatureFunctor is poorly formed" );

			auto ToEigen = []( const SquareMatrix< Real , K > & M )
				{
					Eigen::Matrix< Real , K , K > _M;
					for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) _M(i,j) = M(j,i);
					return _M;
				};

			auto FromEigen = []( const Eigen::Matrix< Real , K , K > & M )
				{
					SquareMatrix< Real , K > _M;
					for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=0 ; j<K ; j++ ) _M(i,j) = M(j,i);
					return _M;
				};

			ConstPointer( SimplexIndex< K , int > ) simplices = mesh.triangles();
			unsigned int simplexNum = (unsigned int)mesh.tCount();
			ThreadPool::ParallelFor
			(
				0 , simplexNum ,
				[&]( unsigned int , size_t i )
				{
					Simplex< Real , Dim , K > s;
					Point< Real , Dim > n[K+1];
					for( unsigned int k=0 ; k<=K ; k++ ) s[k] = V( simplices[i][k] ) , n[k] = N( simplices[i][k] );
					SquareMatrix< Real , K > II = SecondFundamentalForm< Real , K >( s , n );
					Eigen::GeneralizedSelfAdjointEigenSolver< Eigen::MatrixXd > ges( ToEigen( II ) , ToEigen( mesh.g(i) ) );
					Eigen::Matrix< Real , K  , K > D( Eigen::DiagonalMatrix< Real , K >( ges.eigenvalues() ) );
					Point< Real , K > pCurvatures;
					for( unsigned int k=0 ; k<K ; k++ ) pCurvatures[k] = D(k,k);
					pCurvatures = PCF( i , pCurvatures );
					for( unsigned int k=0 ; k<K ; k++ ) D(k,k) = pCurvatures[k];
					mesh.g(i) = mesh.g(i) * FromEigen( Eigen::Matrix< Real , K , K >( ges.eigenvectors() * D * ges.eigenvectors().inverse() ) );
				}
			);
		}
	}
}
#endif // CURVATURE_METRIC_INCLUDED