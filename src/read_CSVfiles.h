/*
 * ReadCSVfiles.h
 *
 * groups functions for reading csv files into arrays and matrices
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <sstream>
#include <cstring>

#include "state.hpp"
#include "slide_aux.hpp"

template <typename Tpath, typename T, size_t ROW, size_t COL>
void loadCSV_mat(const Tpath &name, slide::Matrix<T, ROW, COL> &x)
{
	// 	 * Reads a matrix CSV file with multiple columns.
	// 	 * IN
	// 	 * name 	the name of the file
	// 	 *
	// 	 * OUT
	// 	 * x 		matrix in which the data from the file will be put.
	// 	 *
	// 	 * THROWS
	// 	 * 2 		could not open the file
	// 	 */
	std::ifstream in(name, std::ios_base::in);

	if (!in.good()) // check if we could open the file
	{
		std::cerr << "Error in ReadCSVfiles::loadCSV_mat. File " << name << " could not be opened\n";
		throw 2;
	}

	char c = '.';

	while (in.peek() <= 32) // Ignore byte order mark (BOM) at the beginning
		in >> c;

	for (size_t i = 0; i < ROW; i++)
	{
		// read all but the last column
		for (size_t j = 0; j < COL - 1; j++)
		{
			in >> x[i][j] >> c; // Read number and comma
		}
		in >> x[i][COL - 1]; // Read only number since there is no comma.
	}

	in.close();
}

template <typename Tpath, typename T, size_t ROW>
void loadCSV_1col(const Tpath &name, std::array<T, ROW> &x)
{
	// read data from a CSV file with one column
	slide::Matrix<T, ROW, 1> temp;
	loadCSV_mat(name, temp);
	std::memcpy(&x, &temp, sizeof temp);
}

template <typename Tpath, typename Tx, typename Ty>
void loadCSV_2col(const Tpath &name, Tx &x, Ty &y)
{
	/*
	 * Reads data from a CSV file with 2 columns
	 *
	 * IN
	 * name 	the name of the file	 *
	 * 
	 * OUT
	 * x 		array in which the data from the first column will be put
	 * y 		array in which the data from the second column will be put
	 *
	 * THROWS
	 * 2 		could not open the specified file
	 */

	std::ifstream in(name, std::ios_base::in);

	if (!in.good()) // check if we could open the file
	{
		std::cerr << "Error in ReadCSVfiles::loadCSV_2col. File " << name << " could not be opened.\n";
		throw 2;
	}

	int j = 0;
	char c = '.';

	while (!std::isalnum(in.peek())) // Ignore byte order mark (BOM) at the beginning
		in >> c;

	if constexpr (std::is_same<std::vector<double>, Tx>::value)
	{
		x.clear(); // Sometimes pre-allocated vectors are passed; therefore, cleared to be able to use push_back.
		y.clear();

		double x_i, y_i;
		while (in >> x_i >> c >> y_i) // Read file.
		{
			x.push_back(x_i);
			y.push_back(y_i);
		}
	}
	else
	{									// It must be a std::array, then just read without clear.
		while (in >> x[j] >> c >> y[j]) // Read file.
			j++;
	}
}
