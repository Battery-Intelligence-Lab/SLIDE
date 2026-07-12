/**
 * @file PackSolverValidation.cpp
 * @brief Strict-IEEE public scalar validation for the fast pack solver.
 */

#include "PackSolver.hpp"

namespace slide::core {

slide::Status PackSolver::solve(real_t applied_current,
                                PackSolveMode mode,
                                real_t current_tolerance,
                                int max_iterations)
{
  if (!is_finite(applied_current) || !is_finite(current_tolerance))
    return slide::Status::Invalid_parameters;
  return solveImpl(applied_current,
                   mode,
                   current_tolerance,
                   max_iterations,
                   true);
}

slide::Status PackSolver::setRelaxationGain(real_t alpha)
{
  if (!configured_ || !is_finite(alpha) || !(alpha > 0.0 && alpha <= 1.0))
    return slide::Status::Invalid_parameters;
  relaxation_alpha_ = alpha;
  return slide::Status::Success;
}

} // namespace slide::core
