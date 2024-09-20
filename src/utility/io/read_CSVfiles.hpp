/**
 * @file ReadCSVfiles.hpp
 * @brief Functions for reading CSV files into arrays and matrices.
 *
 * This file contains functions for reading CSV files into arrays and matrices.
 *
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 09 Nov 2022
 */

#pragma once

#include "io_util.hpp"
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

namespace slide::io {

/**
 * @brief Get the contents of a file as a string.
 *
 * @tparam Tpath The type of the file path.
 * @param name The path to the file.
 * @return The contents of the file as a string.
 * @throws std::runtime_error if the file could not be opened.
 */
template <typename Tpath>
std::string getFileContents(const Tpath &name)
{
  //!< For more info see: https://insanecoding.blogspot.com/2011/11/how-to-read-in-file-in-c.html
  std::ifstream in(name, std::ios::in | std::ios::binary);
  if (!in.good()) {
    throw std::runtime_error("Error in getFileContents. File " + std::string(name) + " could not be opened.");
  }

  std::ostringstream fileContents;
  fileContents << in.rdbuf();
  in.close();
  return fileContents.str();
}

/**
 * @brief Get the contents of a file.
 *
 * @param name The name of the file.
 * @return The contents of the file.
 * @throws std::runtime_error if the function is not implemented.
 */
inline auto getFile(std::string name)
{
  throw std::runtime_error("getFile function is not implemented yet!");
  static std::map<std::string, std::string> fileMap;
  // #TODO: Implement the function
}

/**
 * @brief Load data from a CSV file into a matrix.
 *
 * @tparam Tpath The type of the file path.
 * @tparam T The type of the matrix elements.
 * @tparam ROW The number of rows in the matrix.
 * @tparam COL The number of columns in the matrix.
 * @param name The path to the CSV file.
 * @param x The matrix to store the loaded data.
 * @throws std::runtime_error if the file could not be opened.
 */
template <typename Tpath, typename T, size_t ROW, size_t COL>
void loadCSV_mat(const Tpath &name, Matrix<T, ROW, COL> &x)
{
  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) {
    throw std::runtime_error("Error in ReadCSVfiles::loadCSV_mat. File " + std::string(name) + " could not be opened.");
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

/**
 * @brief Load data from a CSV file with one column into an array.
 *
 * @tparam Tpath The type of the file path.
 * @tparam T The type of the array elements.
 * @tparam ROW The number of rows in the array.
 * @param name The path to the CSV file.
 * @param x The array to store the loaded data.
 */
template <typename Tpath, typename T, size_t ROW>
void loadCSV_1col(const Tpath &name, std::array<T, ROW> &x)
{
  slide::Matrix<T, ROW, 1> temp;
  loadCSV_mat(name, temp);
  std::memcpy(&x, &temp, sizeof temp);
}

/**
 * @brief Load data from a CSV file with two columns into separate arrays.
 *
 * @tparam Tpath The type of the file path.
 * @tparam Tx The type of the first array elements.
 * @tparam Ty The type of the second array elements.
 * @param name The path to the CSV file.
 * @param x The array to store the data from the first column.
 * @param y The array to store the data from the second column.
 * @param n The number of rows to read (0 to read all rows).
 * @throws std::runtime_error if the file could not be opened.
 */
template <typename Tpath, typename Tx, typename Ty>
void loadCSV_2col(const Tpath &name, Tx &x, Ty &y, int n = 0)
{
  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) {
    throw std::runtime_error("Error in ReadCSVfiles::loadCSV_2col. File " + std::string(name) + " could not be opened.");
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

/**
 * @brief Load data from a CSV file with multiple columns into a dynamic matrix.
 *
 * @tparam Tpath The type of the file path.
 * @tparam Tx The type of the matrix elements.
 * @param name The path to the CSV file.
 * @param x The dynamic matrix to store the loaded data.
 * @param n The number of rows to read (0 to read all rows).
 * @throws std::runtime_error if the file could not be opened.
 */
template <typename Tpath, typename Tx>
void loadCSV_Ncol(const Tpath &name, DynamicMatrix<Tx> &x, int n = 0)
{
  std::ifstream in(name, std::ios_base::in);

  if (!in.good()) {
    throw std::runtime_error("Error in ReadCSVfiles::loadCSV_2col. File " + std::string(name) + " could not be opened.");
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

/**
 * @brief Structure to store plain X and Y data.
 */
struct XYplain
{
  std::vector<double> x_vec, y_vec;
};

/**
 * @brief Load data from a CSV file with two columns into spans.
 *
 * @tparam Tpath The type of the file path.
 * @param name The path to the CSV file.
 * @param x The span to store the data from the first column.
 * @param y The span to store the data from the second column.
 * @param n The number of rows to read (0 to read all rows).
 * @throws std::runtime_error if the file could not be opened.
 */
template <typename Tpath>
void loadCSV_2col(const Tpath &name, std::span<double> &x, std::span<double> &y, int n = 0)
{
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

} // namespace slide::io