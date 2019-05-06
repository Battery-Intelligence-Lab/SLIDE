/*
 * ReadCSVfiles.cpp
 *
 * groups functions for reading csv files into arrays and matrices
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "ReadCSVfiles.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

using namespace std;

void loadCSV_1col(string name, int nin, double x[]){
	/*
	 * read data from a CSV file with one column
	 *
	 * IN
	 * name 	the name of the file
	 * nin		the length of the file
	 *
	 * OUT
	 * x 		array in which the data will be put
	 *
	 * THROWS
	 * 2 		could not open the specified file
	 */
	ifstream inputfilex;
	inputfilex.open(name);

	// check if we could open the file
	if ( (inputfilex.rdstate() & std::ifstream::failbit ) != 0 ){
		cerr << "Error in ReadCSVfiles::loadCSV_1col. File " << name << " could not be opened\n";
		throw 2;
	}

	// variables
	int j = 0;							// row-index
	string ss;							// string in which the value will be read
	const char* te;						// char array to which the string will be converted
	double nn;							// double to which the char array will be converted

	// loop through the rows of the file
	while(!inputfilex.eof()) {			// read while there are still lines to read
		inputfilex>>ss;					// read data into a string
		te = ss.c_str();				// convert to a c string
		nn = strtod(te, NULL);			// convert to a double
		x[j] = nn;						// store in the array
		j++;							// increase the row index
		if (j == nin)					// reached the length of the file
			break;
	}
}

void loadCSV_2col(string name, int nin, double x[], double y[]){
	/*
	 * Reads data from a CSV file with 2 columns
	 *
	 * IN
	 * name 	the name of the file
	 * nin		the length of the file
	 *
	 * OUT
	 * x 		array in which the data from the first column will be put
	 * y 		array in which the data from the second column will be put
	 *
	 * THROWS
	 * 2 		could not open the specified file
	 */

	ifstream in;
	in.open(name);

	// check if we could open the file
	if ( (in.rdstate() & std::ifstream::failbit ) != 0 ){
		cerr << "Error in ReadCSVfiles::loadCSV_2col. File " << name << " could not be opened\n";
		throw 2;
	}

	int j = 0;
	double nn;
	string ss;
	const char* te;

	while(!in.eof()) {
		getline( in, ss, ',' );			// read until first , (i.e. the first column)
		te = ss.c_str();
		nn = strtod(te, NULL);
		x[j] = nn;

		in>>ss;							// read the rest of the line (i.e. the second column)
		te = ss.c_str();
		nn = strtod(te, NULL);
		y[j] = nn;

		j++;
		if (j == nin)
			break;
	}
}

void loadCSV_2colMatrix(string name, int nin, double x[][2]){
	/*
	 * Reads a csv file with 2 columns and stores the result in one matrix-variable
	 *
	 * IN
	 * name 	the name of the file
	 * nin		the length of the file
	 *
	 * OUT
	 * x 		matrix in which the data will be put
	 * 			x[0][0] has the data from the first column
	 * 			x[0][1] has the data from the second column
	 *
	 * THROWS
	 * 2 		could not open the specified file
	 */

	// Call the function to read the data in two arrays
	double x1[nin], x2[nin];
	try{
		loadCSV_2col(name, nin, x1, x2);
	}
	catch(int e){
		cout<<"Error in ReadCSVfiles::loadCSV_2colMatrix when reading the file: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// Combine both arrays into the matrix
	for(int i=0;i<nin;i++){
		x[i][0] = x1[i];
		x[i][1] = x2[i];
	}
}


void loadCSV_mat(string name, int length, int width, double x[][5*nch]){
	/*
	 * Reads a matrix CSV file with multiple columns.
	 * C++ doesn't allow a variable column number, so the matrix must have 5*nch columns.
	 * Only the first number of columns are filled in, the rest is left as they were
	 * If the file has more columns, an error is thrown.
	 * If you want to read a file with more column, you have to increase the number of columns of the matrix provided (x).
	 *
	 * IN
	 * name 	the name of the file
	 * length	the length of the file [number of rows]
	 * width	the width of the file [number of columns]
	 *
	 * OUT
	 * x 		matrix in which the data from the file will be put.
	 * 			Variable number of rows, 5*nch columns
	 *
	 * THROWS
	 * 2 		could not open the file
	 * 3 		the data file has too many columns
	 */

	ifstream in;
	in.open(name);

	// check if we could open the file
	if ( (in.rdstate() & std::ifstream::failbit ) != 0 ){
		cerr << "Error in ReadCSVfiles::loadCSV_mat. File " << name << " could not be opened\n";
		throw 2;
	}

	// we must provide a bound for the second dimension in the function declaration.
	// Check that this bound is large enough
	if (width >= 5*nch){
		cerr<<"Error in ReadCSVfiles::loadCSV_mat. The input data file has "<<width<<" columns, while it should have fewer than "<<5*nch<<endl;
		throw 3;
	}

	double nn;
	string ss;
	const char* te;

	for(int i=0;i<length;i++){

		// read all but the last column
		for(int j=0;j<width-1;j++){			// loop to read the first until the one but last column
			getline( in, ss, ',' );			// read until next , (i.e. read the next column)
			te = ss.c_str();
			nn = strtod(te, NULL);
			x[i][j] = nn;
		}

		// read the last column
		in>>ss;								// read the rest of the line, including the /n
		te = ss.c_str();
		nn = strtod(te, NULL);
		x[i][width-1] = nn;
	}
	in.close();
}

