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

#ifndef HISTOGRAMS_INCLUDED
#define HISTOGRAMS_INCLUDED
#include <stdio.h>
#include <Misha/Array.h>

namespace Histogram
{
	template< class Real > struct Sample{ Real weight , value; };
	template< class Real > inline void Print( ConstPointer( Real ) histogram , int bins , int units , Real min , Real max );
	template< class Real > inline void Print( const std::vector<         Real   >& values , int bins , int units , Real xMin , Real xMax );
	template< class Real > inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , Real xMin , Real xMax );
	template< class Real > inline void Print( const std::vector< Real >& values , int bins , int units , bool symmetric=false );
	template< class Real > inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , bool symmetric=false );

	template< class Real >
	inline void Print( ConstPointer( Real ) histogram , int bins , int units , Real xMin , Real xMax )
	{
		Real yMax = 0;
		for( int i=0 ; i<bins ; i++ ) yMax = std::max< Real >( yMax , histogram[i] );
		printf( "[ %f , %f ] x [ 0 , %f ]\n" , xMin , xMax , yMax );
		for( int i=0 ; i<units ; i++ ) 
		{
			Real cutOff = (Real) ( yMax / units * ( units - i ) );
			for( int j=0 ; j<bins ; j++ )
			{
				if( histogram[j]>=cutOff ) printf( " " );
				else                       printf( "*" );
			}
			printf( "\n" );
		}
		for( int j=0 ; j<bins ; j++ )
		{
			Real start = xMin + ( xMax - xMin ) / bins * j , end = xMin + ( xMax - xMin ) / bins * ( j+1 );
			if( start<0 && end>=0 ) printf( "0" );
			else                    printf( "-" );
		}
		printf( "\n" );
	}
	template< class Real >
	inline void Print( const std::vector< Real >& values , int bins , int units , Real xMin , Real xMax )
	{
		Pointer( Real ) histogram = AllocPointer< Real >( bins );
		memset( histogram , 0 , sizeof(Real)* bins );
		for( int i=0 ; i<values.size() ; i++ )
		{
			Real b = ( values[i] - xMin ) / ( xMax - xMin ) * ( bins-1 );
			int ib1 = (int)floor( b ) , ib2 = ib1+1;
			Real d1 = (Real)( 1. - (b-ib1) ) , d2 = (Real)( b - ib1 );
			if( ib1>=0 && ib1<bins ) histogram[ib1] += d1;
			if( ib2>=0 && ib2<bins ) histogram[ib2] += d2;
		}
		for( int i=0 ; i<bins ; i++ ) histogram[i] /= values.size();
		Print( histogram , bins , units , xMin , xMax );
		FreePointer( histogram );
	}
	template< class Real >
	inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , Real xMin , Real xMax )
	{
		Real* histogram = AllocPointer< Real >( bins );
		memset( histogram , 0 , sizeof(Real)* bins );
		Real weightSum = 0;
		for( int i=0 ; i<values.size() ; i++ )
		{
			Real b = ( values[i].value - xMin ) / ( xMax - xMin ) * ( bins-1 );
			int ib1 = (int)floor( b ) , ib2 = ib1+1;
			Real d1 = (Real)( 1. - (b-ib1) ) * values[i].weight , d2 = (Real)( b - ib1 ) * values[i].weight;
			weightSum += values[i].weight;
			if( ib1>=0 && ib1<bins ) histogram[ib1] += d1;
			if( ib2>=0 && ib2<bins ) histogram[ib2] += d2;
		}
		for( int i=0 ; i<bins ; i++ ) histogram[i] /= weightSum;
		Print( histogram , bins , units , xMin , xMax );
		FreePointer( histogram );
	}
	template< class Real >
	inline void Print( const std::vector< Real >& values , int bins , int units , bool symmetric )
	{
		Real xMin , xMax;
		xMin = xMax = values[0];
		for( int i=0 ; i<values.size() ; i++ ) xMin = std::min< Real >( xMin , values[i] ) , xMax = std::max< Real >( xMax , values[i] );
		if( symmetric ) xMin = std::min< Real >( xMin , -xMax ) , xMax = std::max< Real >( xMax , -xMin );
		Print( values , bins , units , xMin , xMax );
	}
	template< class Real >
	inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , bool symmetric )
	{
		Real xMin , xMax;
		xMin = xMax = values[0];
		for( int i=0 ; i<values.size() ; i++ ) xMin = std::min< Real >( xMin , values[i].value ) , xMax = std::max< Real >( xMax , values[i].value );
		if( symmetric ) xMin = std::min< Real >( xMin , -xMax ) , xMax = std::max< Real >( xMax , -xMin );
		Print( values , bins , units , xMin , xMax );
	}
}

#endif // HISTOGRAMS_INCLUDED