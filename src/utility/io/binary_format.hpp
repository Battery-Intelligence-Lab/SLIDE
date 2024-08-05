/**
 * @file binary_format.hpp
 * @brief binary writing/reading functions, defines `.slide` format
 * @author Volkan Kumtepeli
 * @date 05 Aug 2024
 */

#pragma once

#include <bit>
#include <fstream>
#include <string>
#include <type_traits>
#include <fstream>
#include <vector>
#include <cstring>

namespace slide::io {

/**
 * @brief Function to check the type of ValueType and return a character indicating the data type.
 *
 * @tparam T_vec The vector type.
 * @return char 'u' for unsigned integer, 'i' for signed integer, 'f' single and double precision float.
 */
template <typename ValueType>
constexpr char getDataType()
{
  if constexpr (std::is_unsigned_v<ValueType>) {
    return 'u';
  } else if constexpr (std::is_signed_v<ValueType> && std::is_integral_v<ValueType>) {
    return 'i';
  } else if constexpr (std::is_same_v<ValueType, float> || std::is_same_v<ValueType, double>) {
    return 'f';
  } else {
    return 'x'; // Unknown type
  }
}

template <typename T_name, typename T_vec, typename T_header>
void binary_writer(T_name &&name, T_vec &data, T_header &header,
                   std::ios_base::openmode mode = std::ios::out)
{
  //!< Create metadata:
  constexpr char version = 0; // Versioning for compatability.
  constexpr char isLittleEndian = static_cast<char>(std::endian::native == std::endian::little);
  const char dType = getDataType<typename T_vec::value_type>();
  const char dSize = static_cast<char>(sizeof(typename T_vec::value_type));

  std::ofstream out{ name, mode | std::ios::binary }; //!< Open file

  out << version << isLittleEndian << dType << dSize; ///!< Write metadata

  for (size_t i = 0; i < header.size(); i++) {
    if (i != 0)
      out << ',';
    if constexpr (std::is_same_v<typename T_header::value_type, std::string>) {
      out.write(header[i].c_str(), header[i].size());
    } else if constexpr (std::is_same_v<typename T_header::value_type,
                                        const char *>) {
      out.write(header[i], std::strlen(header[i]));
    }
  }

  out << '\0'; // So we know when header is over.
  out.write(reinterpret_cast<const char *>(data.data()),
            data.size() * sizeof(typename T_vec::value_type));
  out.close();
}
} // namespace slide::io
