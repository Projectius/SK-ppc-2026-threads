#include <gtest/gtest.h>

//#include "konstantinov_s_graham/all/include/ops_all.hpp"
#include "konstantinov_s_graham/common/include/common.hpp"
//#include "konstantinov_s_graham/omp/include/ops_omp.hpp"
#include "konstantinov_s_graham/seq/include/ops_seq.hpp"
//#include "konstantinov_s_graham/stl/include/ops_stl.hpp"
//#include "konstantinov_s_graham/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace konstantinov_s_graham {

class KonstantinovSRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  InType test_input_;
  OutType test_expected_output_;

  void SetUp() override {
    //input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data==test_expected_output_;
  }

  InType GetTestInputData() final {
    return test_input_;
  }
};

TEST_P(KonstantinovSRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType,  KonstantinovAGrahamSEQ
                                >(PPC_SETTINGS_konstantinov_s_graham);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KonstantinovSRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KonstantinovSRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace konstantinov_s_graham
