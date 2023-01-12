/*
 * ReadCSVfiles.hpp
 *
 * groups functions for reading csv files into arrays and matrices
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "../../types/matrix.hpp"
#include "../util_debug.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <vector>
#include <map>
#include <span>
#include <cstring> //!< -> for memcpy

namespace slide {
template <typename Tpath>
std::string getFileContents(const Tpath &name)
{

  std::cerr << "NOT IMPLEMENTED YET!\n";
  throw 1234;

  //!< For more info see: https://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html

  std::ifstream in(name, std::ios::in | std::ios::binary);
  if (!in.good()) //!< check if we could open the file
  {
    std::cerr << "Error in getFileContents. File " << name << " could not be opened.\n";
    throw 2;
  }

  std::ostringstream fileContents;
  fileContents << in.rdbuf();
  in.close();
  return (fileContents.str());
}

inline auto getFile(std::string name)
{
  std::cerr << "NOT IMPLEMENTED YET!\n";
  throw 1234;
  static std::map<std::string, std::string> fileMap;
}

template <typename Tpath, typename T, size_t ROW, size_t COL>
void loadCSV_mat(const Tpath &name, Matrix<T, ROW, COL> &x)
{
  //!< 	 * Reads a matrix CSV file with multiple columns.
  //!< 	 * IN
  //!< 	 * name 	the name of the file
  //!< 	 *
  //!< 	 * OUT
  //!< 	 * x 		matrix in which the data from the file will be put.
  //!< 	 *
  //!< 	 * THROWS
  //!< 	 * 2 		could not open the file
  //!< 	 */
  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) //!< check if we could open the file
  {
    std::cerr << "Error in ReadCSVfiles::loadCSV_mat. File " << name << " could not be opened\n";
    throw 2;
  }

  char c = '.';

  while (in.peek() <= 32) //!< Ignore byte order mark (BOM) at the beginning
    in >> c;

  for (size_t i = 0; i < ROW; i++) {
    //!< read all but the last column
    for (size_t j = 0; j < COL - 1; j++) {
      in >> x[i][j] >> c; //!< Read number and comma
    }
    in >> x[i][COL - 1]; //!< Read only number since there is no comma.
  }

  in.close();
}

template <typename Tpath, typename T, size_t ROW>
void loadCSV_1col(const Tpath &name, std::array<T, ROW> &x)
{
  //!< read data from a CSV file with one column
  slide::Matrix<T, ROW, 1> temp;
  loadCSV_mat(name, temp);
  std::memcpy(&x, &temp, sizeof temp);
}

inline void ignoreBOM(std::ifstream &in)
{
  char c = '.';

  do {
    in >> c;
  } while (c < 45); //!< Ignore byte order mark (BOM) at the beginning

  in.putback(c);
}

template <typename Tpath, typename Tx, typename Ty>
void loadCSV_2col(const Tpath &name, Tx &x, Ty &y, int n = 0)
{
  /*
   * Reads data from a CSV file with 2 columns
   *
   * IN
   * name 	the name of the file
   * n 		the number of rows to read (if n==0, read all), it is used to read a portion of a *.csv file.
   *
   * OUT
   * x 		array in which the data from the first column will be put
   * y 		array in which the data from the second column will be put
   *
   * THROWS
   * 2 		could not open the specified file
   */

  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) //!< check if we could open the file
  {
    std::cerr << "Error in ReadCSVfiles::loadCSV_2col. File " << name
              << " could not be opened.\n";
    throw 2;
  }

  ignoreBOM(in);

  int j = 0;
  char c;
  if constexpr (std::is_same<std::vector<double>, Tx>::value) {
    x.clear(); //!< Sometimes pre-allocated vectors are passed; therefore, cleared to be able to use push_back.
    y.clear();

    double x_i, y_i;
    while ((n == 0 || j < n) && (in >> x_i >> c >> y_i)) //!< Read file.
    {
      x.push_back(x_i);
      y.push_back(y_i);
      j++;
    }
  } else {                                                 //!< It must be a std::array, then just read without clear.
    while ((n == 0 || j < n) && (in >> x[j] >> c >> y[j])) //!< Read file.
      j++;
  }
}

template <typename Tpath, typename Tx>
void loadCSV_Ncol(const Tpath &name, DynamicMatrix<Tx> &x, int n = 0)
{
  /*
   * Reads data from a CSV file with 2 columns
   *
   * IN
   * name 	the name of the file
   * n 		the number of rows to read (if n==0, read all), it is used to read a portion of a *.csv file.
   *
   * OUT
   * x 		array in which the data from the first column will be put
   * y 		array in which the data from the second column will be put
   *
   * THROWS
   * 2 		could not open the specified file
   */

  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) //!< check if we could open the file
  {
    std::cerr << "Error in ReadCSVfiles::loadCSV_2col. File " << name
              << " could not be opened.\n";
    throw 2;
  }

  ignoreBOM(in);

  x.data.clear(); //!< Sometimes pre-allocated vectors are passed; therefore, cleared to be able to use push_back.

  std::string line;

  int n_rows{ 0 };
  while ((n == 0 || n_rows < n) && std::getline(in, line)) //!< Read file.
  {
    n_rows++;
    std::istringstream in_line(line);
    double x_i;
    char c;
    while (in_line >> x_i) {
      x.data.push_back(x_i);
      in_line >> c;
    }
  }

  int n_cols = x.data.size() / n_rows;

  x.reshape(n_cols, n_rows);

  x.data.shrink_to_fit();
}

struct XYplain
{
  std::vector<double> x_vec, y_vec;
};

//!< struct Data
template <typename Tpath>
void loadCSV_2col(const Tpath &name, std::span<double> &x, std::span<double> &y, int n = 0)
{
  /*
   * Reads data from a CSV file with 2 columns
   *
   * IN
   * name 	the name of the file
   * n 		the number of rows to read (if n==0, read all), it is used to read a portion of a *.csv file.
   *
   * OUT
   * x 		array in which the data from the first column will be put
   * y 		array in which the data from the second column will be put
   *
   * THROWS
   * 2 		could not open the specified file
   */

  static std::map<std::string, XYplain> XYdataMap;

  auto name_str = name.string();

  auto fm = XYdataMap.find(name_str);

  if (fm == XYdataMap.end()) {
    XYplain xyp{};
    loadCSV_2col(name, xyp.x_vec, xyp.y_vec); //!< #TODO -> for some reason it does not take XYdataMap[name] directly.

    XYdataMap[name_str] = std::move(xyp);
  }

  if (n == 0) {
    x = std::span<double>(XYdataMap[name_str].x_vec);
    y = std::span<double>(XYdataMap[name_str].y_vec);
  } else {
    x = std::span<double>(&XYdataMap[name_str].x_vec[0], &XYdataMap[name_str].x_vec[0] + n); // #TODO temporary solution using pointers instead of iterators.
    y = std::span<double>(&XYdataMap[name_str].y_vec[0], &XYdataMap[name_str].y_vec[0] + n);
  }
}

} // namespace slide