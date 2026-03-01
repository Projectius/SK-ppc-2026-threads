#pragma once

#include "konstantinov_s_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace konstantinov_s_graham {

class KonstantinovAGrahamSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KonstantinovAGrahamSEQ(const InType &in);
  static constexpr double K_EPS = 1e-10;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  void remove_duplicates(std::vector<double> &xs, std::vector<double> &ys);
  size_t find_anchor_index(const std::vector<double> &xs, const std::vector<double> &ys);
  double dist2(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j);
  double cross_val(const std::vector<double> &xs, const std::vector<double> &ys, size_t i, size_t j, size_t k);
  std::vector<size_t> collect_and_sort_indices(const std::vector<double> &xs, const std::vector<double> &ys,
                                               size_t anchor_idx);
  bool all_collinear_with_anchor(const std::vector<double> &xs, const std::vector<double> &ys, size_t anchor_idx,
                                 const std::vector<size_t> &sorted_idxs);
std::vector<std::pair<double,double>> build_hull_from_sorted(const std::vector<double> &xs,
                                                             const std::vector<double> &ys,
                                                             size_t anchor_idx,
                                                             const std::vector<size_t> &sorted_idxs);                                 
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace konstantinov_s_graham
