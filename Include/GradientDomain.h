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

#ifndef GRADIENT_DOMAIN_INCLUDED
#define GRADIENT_DOMAIN_INCLUDED

#include <Eigen/Dense>
#include <Misha/MultiThreading.h>
#include <Misha/Geometry.h>
#include <Misha/FEM.h>

namespace GradientDomain
{
	template< typename EigenSolver , typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ , typename NormalFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ >
	std::vector< Point< Real , 3 > > FitToNormals( const FEM::RiemannianMesh< Real > &mesh , Real vWeight , Real nWeight , VertexFunctor && V , NormalFunctor && N );

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High );

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High );

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Variants of the above with a solver for which symbolic factorization has already been performed
	template< typename EigenSolver , typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ , typename NormalFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ >
	std::vector< Point< Real , 3 > > FitToNormals( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real vWeight , Real nWeight , VertexFunctor && V , NormalFunctor && N );

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High );

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High );
	// Variants of the above with a solver for which symbolic factorization has already been performed
	//////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////
	// Implementation //
	////////////////////

	template< typename T >
	struct Array{ static const unsigned int Dim=0; };

	template<>
	struct Array< float >
	{
		using T = float;
		using Real = float;
		static const unsigned int Dim=1;
		static const Real & Entry( const T &v , unsigned int ){ return v; }
		static       Real & Entry(       T &v , unsigned int ){ return v; }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
	};

	template<>
	struct Array< double >
	{
		using T = double;
		using Real = double;
		static const unsigned int Dim=1;
		static const Real & Entry( const T &v , unsigned int ){ return v; }
		static       Real & Entry(       T &v , unsigned int ){ return v; }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
	};

	template< typename _T , unsigned int _Dim , typename _Real >
	struct Array< Point< _T , _Dim , _Real > >
	{
		using T = Point< _T , _Dim , _Real >;
		using Real = _Real;
		static const unsigned int Dim = _Dim * Array< _T >::Dim;
		static const Real & Entry( const T &v , unsigned int d ){ return Array< _T >::Entry( v[ d%_Dim ] , d/_Dim ); }
		static       Real & Entry(       T &v , unsigned int d ){ return Array< _T >::Entry( v[ d%_Dim ] , d/_Dim ); }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return Array< _T >::template Entry< D/_Dim >( v[ D%_Dim ] ); }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return Array< _T >::template Entry< D/_Dim >( v[ D%_Dim ] ); }
	};

	template< unsigned int D , typename T , typename Real , typename EigenSolver , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	void _Solve( const EigenSolver &solver , const Eigen::SparseMatrix< Real > &lowM , const Eigen::SparseMatrix< Real > &highM , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High , std::vector< T > &out  )
	{
		if constexpr( D < Array< T >::Dim )
		{
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > l , h;
			l.resize( lowM.cols() ) , h.resize( highM.cols() );
			for( unsigned int i=0 ; i< lowM.cols() ; i++ ) l[i] = Array< T >::template Entry< D >(  Low(i) );
			for( unsigned int i=0 ; i<highM.cols() ; i++ ) h[i] = Array< T >::template Entry< D >( High(i) );

#ifdef EIGEN_USE_MKL_ALL
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > b = lowM * l + highM * h;
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > x = solver.solve( b );
#else // !EIGEN_USE_MKL_ALL
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > x = solver.solve( lowM * l + highM * h );
#endif // EIGEN_USE_MKL_ALL

			for( unsigned int i=0 ; i<x.size() ; i++ ) Array< T >::template Entry< D >( out[i] ) = x[i];
			_Solve< D+1 >( solver , lowM , highM , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , out );
		}
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > _ProcessVertexVertex( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High , bool preAnalyzed )
	{
		static_assert( Array< T >::Dim , "[ERROR] T is not of array type" );
		static_assert( std::is_convertible_v< LowFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] LowFrequencyVertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< HighFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] HighFrequencyVertexFunctor is poorly formed" );

		Eigen::SparseMatrix< Real > mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY , true >() * lowWeight;
		Eigen::SparseMatrix< Real > stiffness = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY , true >() * highWeight;

		if( preAnalyzed ) solver.factorize( mass + stiffness );
		else              solver.compute  ( mass + stiffness );
		if( solver.info()!=Eigen::Success ) ERROR_OUT( "Failed to factorize matrix" );

		std::vector< T > out( mesh.vCount() );
		_Solve< 0 >( solver , mass , stiffness , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , out );
		return out;
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High )
	{
		EigenSolver solver;
		return _ProcessVertexVertex< EigenSolver , T , Real >( solver , mesh , lowWeight , highWeight , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , false );
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High )
	{
		return _ProcessVertexVertex< EigenSolver , T , Real >( solver , mesh , lowWeight , highWeight , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , true );
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > _ProcessVertexEdge( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High , bool preAnalyzed )
	{
		static_assert( Array< T >::Dim , "[ERROR] T is not of array type" );
		static_assert( std::is_convertible_v< LowFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] LowFrequencyVertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< HighFrequencyEdgeFunctor , std::function< T ( unsigned int ) > > , "[ERROR] HighFrequencyEdgeFunctor is poorly formed" );

		Eigen::SparseMatrix< Real > mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY , true >() * lowWeight;
		Eigen::SparseMatrix< Real > d = mesh.template dMatrix< FEM::BASIS_0_WHITNEY , FEM::BASIS_1_WHITNEY , true >();
		Eigen::SparseMatrix< Real > divergence = d.transpose() * mesh.template massMatrix< FEM::BASIS_1_WHITNEY , true >( true ) * highWeight;

		if( preAnalyzed ) solver.factorize( mass + divergence * d );
		else              solver.compute  ( mass + divergence * d );
		if( solver.info()!=Eigen::Success ) ERROR_OUT( "Failed to factorize matrix" );

		std::vector< T > out( mesh.vCount() );
		_Solve< 0 >( solver , mass , divergence , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyEdgeFunctor >( High ) , out );
		return out;
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High )
	{
		EigenSolver solver;
		return _ProcessVertexEdge< EigenSolver , T , Real >( solver , mesh , lowWeight , highWeight , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyEdgeFunctor >( High ) , false );
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High )
	{
		return _ProcessVertexEdge< EigenSolver , T , Real >( solver , mesh , lowWeight , highWeight , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyEdgeFunctor >( High ) , true );
	}

	template< typename EigenSolver , typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ , typename NormalFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ >
	std::vector< Point< Real , 3 > > FitToNormals( const FEM::RiemannianMesh< Real > &mesh , Real vWeight , Real nWeight , VertexFunctor && V , NormalFunctor && N )
	{
		static_assert( std::is_convertible_v< VertexFunctor , std::function< Point< Real , 3 > ( unsigned int ) > > , "[ERROR] VertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< NormalFunctor , std::function< Point< Real , 3 > ( unsigned int ) > > , "[ERROR] NormalFunctor is poorly formed" );

		auto LowFrequencyVertexValues = [&]( unsigned int v ){ return V(v); };

		// High frequencies given projection of edge offset onto the normal plane
		auto HighFrequencyEdgeValues = [&]( unsigned int e )
			{
				int v1 , v2;
				mesh.edgeVertices( e , v1 , v2 );
				Point< Real , 3 > d = V(v2) - V(v1);
				Point< Real , 3 > n = N(v2) + N(v1);
				return d - Point< Real , 3 >::Dot( d , n ) / Point< Real , 3 >::SquareNorm( n ) * n;
			};

		return ProcessVertexEdge< EigenSolver , Point3D< Real > , Real >( mesh , vWeight , nWeight , LowFrequencyVertexValues , HighFrequencyEdgeValues );
	}

	template< typename EigenSolver , typename Real , typename VertexFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ , typename NormalFunctor /* = std::function< Point< Real , 3 >( unsigned int ) > */ >
	std::vector< Point< Real , 3 > > FitToNormals( EigenSolver & solver , const FEM::RiemannianMesh< Real > &mesh , Real vWeight , Real nWeight , VertexFunctor && V , NormalFunctor && N )
	{
		static_assert( std::is_convertible_v< VertexFunctor , std::function< Point< Real , 3 > ( unsigned int ) > > , "[ERROR] VertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< NormalFunctor , std::function< Point< Real , 3 > ( unsigned int ) > > , "[ERROR] NormalFunctor is poorly formed" );

		auto LowFrequencyVertexValues = [&]( unsigned int v ){ return V(v); };

		// High frequencies given projection of edge offset onto the normal plane
		auto HighFrequencyEdgeValues = [&]( unsigned int e )
			{
				int v1 , v2;
				mesh.edgeVertices( e , v1 , v2 );
				Point< Real , 3 > d = V(v2) - V(v1);
				Point< Real , 3 > n = N(v2) + N(v1);
				return d - Point< Real , 3 >::Dot( d , n ) / Point< Real , 3 >::SquareNorm( n ) * n;
			};

		return ProcessVertexEdge< EigenSolver , Point3D< Real > , Real >( solver , mesh , vWeight , nWeight , LowFrequencyVertexValues , HighFrequencyEdgeValues );
	}

};
#endif // GRADIENT_DOMAIN_INCLUDED