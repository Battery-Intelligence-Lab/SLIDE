/**
 * @file macros.hpp
 * @brief Defines macros, especially used for external libraries.
 * @author Volkan Kumtepeli
 * @date 26 Jul 2024
 */

// Default of this is __forceinline but we use __inline to make compilation faster at the moment.
// Probably better to make it __forceinline in final build. #TODO Make this SLIDE_DEVELOPER_MODE dependent macro.
#if (defined _MSC_VER) || (defined __INTEL_COMPILER)
#define EIGEN_STRONG_INLINE __inline //
#else
#define EIGEN_STRONG_INLINE inline
#endif