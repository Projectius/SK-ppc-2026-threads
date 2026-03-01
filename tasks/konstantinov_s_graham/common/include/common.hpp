#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace konstantinov_s_graham {

using InType = std::pair<std::vector<double>, std::vector<double>>;
using OutType = std::vector<std::pair<double, double>>;
using TestType = std::tuple<InType, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace konstantinov_s_graham
